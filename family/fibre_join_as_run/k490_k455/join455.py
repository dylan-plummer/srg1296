"""Bounded, exact three-domain join for degree 455 with five matching channels."""
from audit import *
import functools


def bitsets(column):
    length=(len(column)+7)//8;buffers={}
    for i,x in enumerate(column):
        key=int(x)
        if key not in buffers:buffers[key]=bytearray(length)
        buffers[key][i//8]|=1<<(i%8)
    return {key:int.from_bytes(value,'little') for key,value in buffers.items()}


def indices(mask):
    while mask:
        bit=mask&-mask;yield bit.bit_length()-1;mask^=bit


def prepare(fibre_mask=124):
    D=data();P=D['P'];cert=json.loads((SOURCE/'result/construction.json').read_text())
    F=np.array(cert['ternary_matrix'],dtype=np.int16);lookup={tuple(p):i for i,p in enumerate(P)};perm=[]
    for p in P:
        v=F.T@p%3;v=v*next(int(x) for x in v if x)%3;perm.append(lookup[tuple(v)])
    frequency=cycles(perm);fc=D['fine_cycles'];z=np.zeros(41,dtype=np.int16)
    for i,orbit in enumerate(fc):
        if fibre_mask>>i&1:z[np.array(orbit)+1]=1
    e=(D['character'][:41,:41]@z)[1:];H=D['H'];patterns={}
    for k in [4,8]:
        patterns[k]=[]
        for S in itertools.combinations(range(16),k):
            if functools.reduce(int.__xor__,S,0)==0:patterns[k].append(H[S,:].sum(0)[1:]//4)
        patterns[k]=np.array(patterns[k],dtype=np.int16)
    assert len(patterns[4])==140 and len(patterns[8])==870
    templates=[]
    for cycle in frequency:
        pats=patterns[8 if e[cycle[0]]==max(e) else 4]
        codes=np.column_stack([(pats[:,np.array(dc)-1]+2)@np.array([5**i for i in range(5)]) for dc in D['binary_cycles']])
        templates.append(codes)
    return D,frequency,e,templates,z


def save_found(D,frequency,e,arrays,choices,z,case,parameters=(1296,455,166,156),fibre_mask=124):
    cert=json.loads((SOURCE/'result/construction.json').read_text());P=D['P'];lookup={tuple(p):i for i,p in enumerate(P)}
    F=np.array(cert['ternary_matrix'],dtype=np.int16);forward=[]
    for p in P:
        v=F@p%3;v=v*next(int(x) for x in v if x)%3;forward.append(lookup[tuple(v)])
    zz=np.zeros((16,41),dtype=np.int16);zz[0]=z;seeds=[]
    for b,choice in enumerate(choices):
        mask=int(arrays[b]['masks'][choice]);d=D['binary_cycles'][b][0];initial=mask
        origin=int(mask.bit_count()==13)
        seeds.append(dict(binary_label=d,projective_mask=hex(mask),include_zero=bool(origin),
            projective_points=[list(map(int,p)) for i,p in enumerate(P) if mask>>i&1]))
        for _ in range(5):
            zz[d,0]=origin
            for p in range(40):zz[d,p+1]=mask>>p&1
            mask=sum(1<<forward[p] for p in range(40) if mask>>p&1);d=D['binary'][d]
        assert mask==initial
    v,k,lam,mu=parameters;gap=math.isqrt((lam-mu)**2+4*(k-mu));r,s=((lam-mu+gap)//2,(lam-mu-gap)//2)
    ev=D['character']@zz.ravel();assert ev[0]==k and set(map(int,ev[1:]))=={r,s}
    fine=D['fine'];fi=D['fi'];connections=[]
    for b in range(16):
        if zz[b,0]:connections.append((b,(0,0,0,0)))
        for p in range(40):
            if zz[b,p+1]:connections.extend([(b,tuple(map(int,P[p]))),(b,tuple(int(-x)%3 for x in P[p]))])
    assert len(connections)==k
    # An independent direct difference count, without spectral assumptions.
    connset=set(connections);counts=collections.Counter((b^c,tuple((x-y)%3 for x,y in zip(v,w))) for b,v in connections for c,w in connections)
    for b in range(16):
        for v in fine:assert counts[b,v]==(k if b==0 and not any(v) else (lam if (b,v) in connset else mu))
    multiplicities=collections.Counter()
    for i,value in enumerate(ev):multiplicities[int(value)]+=1 if i%41==0 else 2
    assert sum(multiplicities.values())==1296
    result=dict(parameters=list(parameters),case=case,fibre_mask=fibre_mask,
        binary_matrix=cert['binary_matrix'],ternary_matrix=cert['ternary_matrix'],
        fibre_projective_points=[list(map(int,p)) for p in P[z[1:]!=0]],seeds=seeds,
        all_1296_PDS_difference_counts_verified=True,character_spectrum=dict(multiplicities))
    (HERE/f'FOUND{k}_certificate.json').write_text(json.dumps(result,indent=2)+'\n')
    (HERE/f'FOUND{k}_connections.txt').write_text(''.join(str(b)+' '+' '.join(map(str,v))+'\n' for b,v in sorted(connections)))
    return result


def main(seconds=30):
    start=time.monotonic();D,frequency,e,templates,z=prepare();out=[]
    for case in range(3):
        source=np.load(HERE/f'matching455_{case}.npz');arrays=[]
        for b in range(3):
            C=source['cap_coefficients' if b==2 else 'zero_coefficients'];masks=source['cap_masks' if b==2 else 'zero_masks']
            codes=np.column_stack([(C[:,cycle]+2)@np.array([5**i for i in range(5)]) for cycle in frequency])
            keep=np.ones(len(C),dtype=bool)
            for f,ts in enumerate(templates):keep &=np.isin(codes[:,f],ts[:,b])
            arrays.append(dict(masks=masks[keep],codes=codes[keep],energy=(C[keep][:,e==5].astype(np.int16)**2).sum(1)))
        A,B,C=[a['codes'] for a in arrays]
        print(json.dumps(dict(case=case,domains=list(map(len,[A,B,C])),seconds=time.monotonic()-start)),flush=True)
        abit=[bitsets(A[:,f]) for f in range(8)];bbit=[bitsets(B[:,f]) for f in range(8)]
        allowed=[];left_by_cap=[];right_by_cap=[]
        for f,ts in enumerate(templates):
            pair=collections.defaultdict(set);left=collections.defaultdict(set);right=collections.defaultdict(set)
            for aa,bb,cc in ts:pair[int(aa),int(cc)].add(int(bb));left[int(cc)].add(int(aa));right[int(cc)].add(int(bb))
            def union(values,index):return functools.reduce(int.__or__,(index.get(v,0) for v in values),0)
            allowed.append({key:union(vals,bbit[f]) for key,vals in pair.items()})
            left_by_cap.append({key:union(vals,abit[f]) for key,vals in left.items()})
            right_by_cap.append({key:union(vals,bbit[f]) for key,vals in right.items()})
        energy_bit=bitsets(arrays[1]['energy']);tested=0;found=None;complete=True
        for ci,cap in enumerate(C):
            liveA=(1<<len(A))-1;liveB=(1<<len(B))-1
            for f in range(8):liveA &=left_by_cap[f].get(int(cap[f]),0);liveB &=right_by_cap[f].get(int(cap[f]),0)
            if not liveA or not liveB:continue
            order=sorted(range(8),key=lambda f:right_by_cap[f].get(int(cap[f]),0).bit_count())
            for ai in indices(liveA):
                tested+=1;alive=liveB &energy_bit.get(12-int(arrays[0]['energy'][ai])-int(arrays[2]['energy'][ci]),0)
                for f in order:
                    alive &=allowed[f].get((int(A[ai,f]),int(cap[f])),0)
                    if not alive:break
                if alive:found=[ai,next(indices(alive)),ci];break
            if found:break
            if time.monotonic()-start>seconds:complete=False;break
        record=dict(case=case,domain_sizes=list(map(len,[A,B,C])),pairs_tested=tested,complete_join=complete and found is None,found=found,
            seconds_since_start=time.monotonic()-start)
        if found:record['certificate']=save_found(D,frequency,e,arrays,found,z,case)
        out.append(record);print(json.dumps(record),flush=True)
        if found or not complete:break
    (HERE/'join455_results.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
