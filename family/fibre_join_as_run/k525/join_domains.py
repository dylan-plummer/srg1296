"""Exact three-bucket join through signed Walsh templates and bitsets."""
from unsigned_domains import *
from cayley_model import save_graph


def bitsets(column):
    n=(len(column)+7)//8;buffers={}
    for i,code in enumerate(column):
        key=int(code)
        if key not in buffers:buffers[key]=bytearray(n)
        buffers[key][i>>3]|=1<<(i&7)
    return {key:int.from_bytes(buf,'little') for key,buf in buffers.items()}


def indices(mask):
    while mask:
        bit=mask&-mask;yield bit.bit_length()-1;mask^=bit


def reconstruct(selected,arrays,D,U,buckets):
    P=np.array(D['points'],dtype=np.int16);V=list(itertools.product(range(3),repeat=4));ix={v:i for i,v in enumerate(V)}
    T=np.zeros((16,81),dtype=np.uint8)
    for p in P:
        if (int(p@p)+int(p.sum())**2)%3==2:T[0,ix[tuple(p)]]=T[0,ix[tuple(-p%3)]]=1
    g=U['action']
    for b,option in enumerate(selected):
        size=buckets[b]['cross_size'];local=int(arrays[b]['local_indices'][option]);mask=int(np.load(HERE/f'local{size}.npz')['masks'][local]);d=U['difference_cycles'][b][0]
        for _ in range(5):
            T[d,0]=int(size==27)
            for j,p in enumerate(P):
                if mask>>j&1:T[d,ix[tuple(p)]]=T[d,ix[tuple(-p%3)]]=1
            mask=sum(1<<g['projective'][j] for j in range(40) if mask>>j&1);d=g['binary'][d]
    assert T.sum()==525
    N=P@P.T%3==0;fine=T[:,0,None]+T[:,[ix[tuple(p)] for p in P]]@np.where(N,2,-1).T
    H=np.array([[(-1)**((a&b).bit_count()%2) for b in range(16)] for a in range(16)],dtype=np.int16)
    assert set(map(int,(H@fine).ravel()))<={21,-15}
    quotient=H@T.sum(1);assert int(quotient[0])==525 and set(map(int,quotient[1:]))<={21,-15}
    return T


def main():
    started=time.monotonic();D=json.loads((HERE/'source.json').read_text());U=json.loads((HERE/'unsigned_domains.json').read_text());S=json.loads((HERE/'signed_domains.json').read_text())
    arrays=[np.load(HERE/f'signed_bucket{b}.npz') for b in range(3)];A,B,C=[a['codes'] for a in arrays]
    assert len(C)<=min(len(A),len(B))
    abit=[bitsets(A[:,f]) for f in range(8)];bbit=[bitsets(B[:,f]) for f in range(8)]
    allowed=[];left_by_cap=[];right_by_cap=[]
    for f,templates in enumerate(S['frequency_templates']):
        pairs=collections.defaultdict(set);left=collections.defaultdict(set);right=collections.defaultdict(set)
        for a,b,c in templates:pairs[a,c].add(b);left[c].add(a);right[c].add(b)
        table={key:__import__('functools').reduce(int.__or__,(bbit[f].get(b,0) for b in vals),0) for key,vals in pairs.items()}
        allowed.append(table)
        left_by_cap.append({key:__import__('functools').reduce(int.__or__,(abit[f].get(a,0) for a in vals),0) for key,vals in left.items()})
        right_by_cap.append({key:__import__('functools').reduce(int.__or__,(bbit[f].get(b,0) for b in vals),0) for key,vals in right.items()})
    local36=np.load(HERE/'local36.npz')['coefficients'];positive=np.array(D['diagonal_fourier'])==3
    energyA=(local36[arrays[0]['local_indices']][:,positive].astype(np.int16)**2).sum(1)
    energyB=(local36[arrays[1]['local_indices']][:,positive].astype(np.int16)**2).sum(1)
    ebit=bitsets(energyB);local27=np.load(HERE/'local27.npz')['coefficients']
    energyC=(local27[arrays[2]['local_indices']][:,positive].astype(np.int16)**2).sum(1)
    tested=0;surviving_left=0;found=None;allA=(1<<len(A))-1;allB=(1<<len(B))-1;stats=[]
    print(json.dumps(dict(stage='indexed',domain_sizes=[len(A),len(B),len(C)],seconds=time.monotonic()-started)),flush=True)
    for ci,cap in enumerate(C):
        liveA=allA;liveB=allB
        for f in range(8):liveA &=left_by_cap[f].get(int(cap[f]),0);liveB &=right_by_cap[f].get(int(cap[f]),0)
        surviving_left+=liveA.bit_count();stats.append([ci,liveA.bit_count(),liveB.bit_count()])
        if not liveA or not liveB:continue
        # Start with the most selective character constraint for this cap.
        order=sorted(range(1,8),key=lambda f:right_by_cap[f].get(int(cap[f]),0).bit_count())
        for ai in indices(liveA):
            tested+=1;alive=liveB & ebit.get(20-int(energyC[ci])-int(energyA[ai]),0)
            if not alive:continue
            for f in order:
                alive &=allowed[f].get((int(A[ai,f]),int(cap[f])),0)
                if not alive:break
            if alive:found=[ai,next(indices(alive)),ci];break
        if found:break
        if ci%24==0:print(json.dumps(dict(cap_processed=ci+1,left_pairs_tested=tested,seconds=time.monotonic()-started)),flush=True)
    out=dict(scope=S['scope'],domain_sizes=[len(A),len(B),len(C)],complete_join=found is None,
             pairs_tested=tested,left_candidates_after_cap_filter=surviving_left,
             found_bucket_choices=found,cap_stats=stats,seconds=time.monotonic()-started)
    if found:
        for f,templates in enumerate(S['frequency_templates']):assert (int(A[found[0],f]),int(B[found[1],f]),int(C[found[2],f])) in set(map(tuple,templates))
        T=reconstruct(found,arrays,D,U,S['buckets']);out['verification']=save_graph(T,np.array(D['points']))
    (HERE/'join_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='cap_stats'}),flush=True)


if __name__=='__main__':main()
