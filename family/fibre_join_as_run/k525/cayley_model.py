"""Exact 120-choice Cayley realization forced by the local uniqueness check.

Keep the fiber, quotient, C5 action and five matching channels. Impose
integer cross Fourier values and every global character equation directly.
The optional Walsh tables strengthen the same spectral constraints.
"""
from unsigned_domains import *
from ortools.sat.python import cp_model
import argparse


def setup():
    D=json.loads((HERE/'source.json').read_text());P=np.array(D['points'],dtype=np.int16)
    N=P@P.T%3==0;actions=D['group_actions'];units={(b,p) for b in range(16) for p in range(40)}|{(b,-1) for b in range(1,16)}
    oo=[];where={}
    while units:
        seed=min(units);orbit=sorted({(g['binary'][seed[0]],-1 if seed[1]<0 else g['projective'][seed[1]]) for g in actions})
        assert set(orbit)<=units;units-=set(orbit)
        for u in orbit:where[u]=len(oo)
        oo.append(orbit)
    assert len(oo)==131
    U=json.loads((HERE/'unsigned_domains.json').read_text())
    return D,P,N,oo,where,U


def build(tables=True,selfdual=False):
    D,P,N,oo,where,U=setup();md=cp_model.CpModel();x=[md.NewBoolVar(f'x{i}') for i in range(len(oo))]
    cap=set(D['binary_cap']);fixed={}
    for i,orbit in enumerate(oo):
        if orbit[0][0]==0:
            vals={int((int(P[p]@P[p])+int(P[p].sum())**2)%3==2) for b,p in orbit};assert len(vals)==1
            fixed[i]=vals.pop()
        elif orbit[0][1]<0:
            vals={int(b in cap) for b,p in orbit};assert len(vals)==1;fixed[i]=vals.pop()
        if i in fixed:md.Add(x[i]==fixed[i])
    # Exact principal degree.
    weights=[sum(1 if p<0 else 2 for b,p in orbit) for orbit in oo]
    md.Add(sum(w*v for w,v in zip(weights,x))==525)
    C={};square={}
    for bucket,dc in enumerate(U['difference_cycles']):
        b=dc[0];size=27 if b in cap else 36;origin=x[where[b,-1]]
        md.Add(origin+2*sum(x[where[b,p]] for p in range(40))==size)
        for q in range(40):
            bound=1 if D['diagonal_fourier'][q]==-6 else 2
            c=md.NewIntVar(-bound,bound,f'c{bucket}_{q}');C[bucket,q]=c
            md.Add(origin+sum((2 if N[q,p] else -1)*x[where[b,p]] for p in range(40))==9*c)
            t=md.NewIntVar(0,4,'');md.AddMultiplicationEquality(t,[c,c]);square[bucket,q]=t
            if q in D['matching_projective_channels']:
                target=2 if size==27 and D['binary_cap'][D['matching_projective_channels'].index(q)]==b else 0
                md.Add(c==target)
    # All full-group characters, merged only if their coefficient vectors
    # are identical. Pairing +/- ternary characters is exact for these sets.
    profiles={}
    for a in range(16):
        for q in range(-1,40):
            if a==0 and q==-1:continue
            row=[]
            for orbit in oo:
                row.append(sum((-1)**((a&b).bit_count()%2)*(1 if p<0 else (2 if q<0 or N[q,p] else -1)) for b,p in orbit))
            profiles.setdefault(tuple(row),(a,q))
    eigen={}
    for row in profiles:
        eig=md.NewBoolVar('');eigen[row]=eig;md.Add(sum(c*v for c,v in zip(row,x))==-15+36*eig)
    # The positive eigenspace has multiplicity 525. These explicit rank
    # counts help propagation beyond the individual character equations.
    profile_for={}
    for a in range(16):
        for q in range(40):
            row=tuple(sum((-1)**((a&b).bit_count()%2)*(1 if p<0 else (2 if N[q,p] else -1)) for b,p in orbit) for orbit in oo)
            profile_for[a,q]=eigen[row]
    for fo in U['frequency_cycles']:
        q=fo[0];md.Add(sum(profile_for[a,q] for a in range(16))==(4 if D['diagonal_fourier'][q]==-6 else 8))
    if selfdual:
        G2=(np.eye(4,dtype=np.int16)+np.ones((4,4),dtype=np.int16))%2
        G3=(np.eye(4,dtype=np.int16)+np.ones((4,4),dtype=np.int16))%3
        lookup={tuple(p):i for i,p in enumerate(P)}
        for i,orbit in enumerate(oo):
            b,p=orbit[0];bits=np.array([(b>>(3-j))&1 for j in range(4)],dtype=np.int16);a=0
            for digit in G2@bits%2:a=2*a+int(digit)
            if p<0:q=-1
            else:
                v=G3@P[p]%3;v=v*int(next(c for c in v if c))%3;q=lookup[tuple(v)]
            row=tuple(sum((-1)**((a&bb).bit_count()%2)*(1 if pp<0 else (2 if q<0 or N[q,pp] else -1)) for bb,pp in block) for block in oo)
            md.Add(x[i]==eigen[row])
    inv=np.argsort(U['action']['dual'])
    for fo in U['frequency_cycles']:
        p=fo[0];idx=[]
        for _ in range(5):idx.append(int(p));p=inv[p]
        md.Add(sum(square[b,q] for b in range(3) for q in idx)==(3 if D['diagonal_fourier'][fo[0]]==-6 else 4))
    if tables:
        H=np.array([[(-1)**((a&b).bit_count()%2) for b in range(16)] for a in range(16)],dtype=np.int16)
        patterns={4:[],8:[]}
        for k in [4,8]:
            for S in itertools.combinations(range(16),k):
                xx=0
                for a in S:xx^=a
                if not xx:patterns[k].append((H[1:,S].sum(1)//4).tolist())
        for fo in U['frequency_cycles']:
            if fo[0] in D['matching_projective_channels']:continue
            refs={}
            for b,dc in enumerate(U['difference_cycles']):
                q=fo[0]
                for d in dc:refs[d]=C[b,int(q)];q=inv[q]
            md.AddAllowedAssignments([refs[d] for d in range(1,16)],patterns[4 if D['diagonal_fourier'][fo[0]]==-6 else 8])
    assert not md.Validate(),md.Validate()
    # Independent coefficient controls: complete multipartite and the
    # prescribed fiber character table, evaluated without a solver.
    vector=np.array([int(b!=0) for orbit in oo for b,p in [orbit[0]]])
    for row,(a,q) in profiles.items():assert int(np.array(row)@vector)==(0 if q>=0 else -81)
    fibre=np.array([fixed.get(i,0) if oo[i][0][0]==0 else 0 for i in range(len(oo))])
    for row,(a,q) in profiles.items():assert int(np.array(row)@fibre)==(30 if q<0 else D['diagonal_fourier'][q])
    return md,x,(D,P,N,oo,where,U),dict(orbit_variables=131,unfixed_connection_choices=131-len(fixed),
                                         distinct_nonprincipal_character_profiles=len(profiles),walsh_tables=tables,selfdual=selfdual,controls_passed=True)


def save_graph(T,P):
    V=list(itertools.product(range(3),repeat=4));index={v:i for i,v in enumerate(V)}
    assert T.shape==(16,81) and int(T.sum())==525 and not T[0,0]
    negative=[index[tuple(-x%3 for x in v)] for v in V];assert np.array_equal(T,T[:,negative])
    rows=[];A=np.zeros((1296,1296),dtype=np.uint8)
    for b in range(16):
        for v in V:
            i=81*b+index[v]
            for c in range(16):
                for w in V:
                    delta=index[tuple((y-x)%3 for x,y in zip(v,w))]
                    A[i,81*c+index[w]]=T[b^c,delta]
            assert int(A[i].sum())==525
            rows.append(sum(int(x)<<j for j,x in enumerate(A[i])))
    assert np.array_equal(A,A.T) and not np.any(np.diag(A))
    for i in range(1296):
        for j in range(i):assert (rows[i]&rows[j]).bit_count()==(216 if A[i,j] else 210)
    np.savetxt(HERE/'FOUND_srg1296_525_216_210.txt',A,fmt='%d')
    np.savetxt(HERE/'FOUND_connections.txt',T,fmt='%d')
    return dict(parameters=[1296,525,216,210],all_degrees_verified=True,all_common_neighbour_pairs_verified=1296*1295//2)


def main(seconds,tables,selfdual=False):
    started=time.monotonic();md,x,data,out=build(tables,selfdual);out['build_seconds']=time.monotonic()-started
    print(json.dumps(out),flush=True)
    sv=cp_model.CpSolver();sv.parameters.num_workers=1;sv.parameters.max_time_in_seconds=seconds
    start=time.monotonic();status=sv.Solve(md)
    out.update(status=sv.StatusName(status),solver_seconds=time.monotonic()-start,branches=sv.NumBranches(),conflicts=sv.NumConflicts())
    if status in [cp_model.FEASIBLE,cp_model.OPTIMAL]:
        D,P,N,oo,where,U=data;V=list(itertools.product(range(3),repeat=4));ix={v:i for i,v in enumerate(V)};T=np.zeros((16,81),dtype=np.uint8)
        for i,orbit in enumerate(oo):
            if not sv.Value(x[i]):continue
            for b,p in orbit:
                if p<0:T[b,0]=1
                else:T[b,ix[tuple(P[p])]]=T[b,ix[tuple(-P[p]%3)]]=1
        out['verification']=save_graph(T,P)
    name='cayley_'+('selfdual_' if selfdual else '')+('tables' if tables else 'linear')+'.json'
    (HERE/name).write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--seconds',type=float,default=30);ap.add_argument('--no-tables',action='store_true')
    ap.add_argument('--selfdual',action='store_true');a=ap.parse_args();main(a.seconds,not a.no_tables,a.selfdual)
