"""Exact finite transfer screens for the newly constructed degree-525 PDS."""
from pathlib import Path
import collections,csv,importlib.util,itertools,json,math,time
import numpy as np
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
SOURCE=ROOT/'work/construct525_20261116'


def cycles(perm):
    remaining=set(range(len(perm)));out=[]
    while remaining:
        x=min(remaining);orbit=[]
        while x not in orbit:orbit.append(x);x=perm[x]
        assert set(orbit)<=remaining;remaining.difference_update(orbit);out.append(orbit)
    return out


def data():
    spec=importlib.util.spec_from_file_location('standalone525',SOURCE/'result/construct.py')
    builder=importlib.util.module_from_spec(spec);spec.loader.exec_module(builder)
    certificate=json.loads((SOURCE/'result/construction.json').read_text())
    rows,connections=builder.build(certificate)
    source=json.loads((SOURCE/'source.json').read_text());P=np.array(source['points'],dtype=np.int16)
    points={tuple(p):i for i,p in enumerate(P)}
    def normalized(v):
        return tuple(x*next(a for a in v if a)%3 for x in v)
    F=certificate['ternary_matrix'];B=certificate['binary_matrix']
    perm=[points[normalized(builder.apply(F,p,3))] for p in P]
    fine_cycles=cycles(perm);assert len(fine_cycles)==8
    binary=[]
    for b in range(16):
        bits=tuple((b>>(3-j))&1 for j in range(4))
        binary.append(builder.number(builder.apply(B,bits,2),2))
    bc=[c for c in cycles(binary) if c!=[0]]
    N=P@P.T%3==0
    H=np.array([[(-1)**((a&b).bit_count()%2) for b in range(16)] for a in range(16)],dtype=np.int16)
    ternary=np.column_stack([np.ones(41,dtype=np.int16),np.vstack([np.full((1,40),2),np.where(N,2,-1)])])
    # Character indices: 41*a + q+1, where q=-1 is zero ternary character.
    character=np.kron(H,ternary)
    T=np.array([[(rows[0]>>(81*b+t))&1 for t in range(81)] for b in range(16)],dtype=np.int16)
    fine=list(itertools.product(range(3),repeat=4));fi={v:i for i,v in enumerate(fine)}
    z=np.zeros((16,41),dtype=np.int16);z[:,0]=T[:,0]
    for j,p in enumerate(P):z[:,j+1]=T[:,fi[tuple(p)]]
    assert int((character@z.ravel())[0])==525
    assert set(map(int,(character@z.ravel())[1:]))=={21,-15}
    return dict(P=P,N=N,H=H,character=character,z=z,T=T,fine=fine,fi=fi,
                fine_cycles=fine_cycles,binary_cycles=bc,binary=binary,fine_perm=perm)


def subspaces(p):
    vectors=list(itertools.product(range(p),repeat=4));seen=set()
    for u,v in itertools.combinations(vectors[1:],2):
        span=frozenset(tuple((a*x+b*y)%p for x,y in zip(u,v)) for a,b in itertools.product(range(p),repeat=2))
        if len(span)==p*p:seen.add(span)
    return sorted(seen,key=lambda s:sorted(s))


def main():
    start=time.monotonic();D=data();C=D['character'];P=D['P'];z=D['z'];fine_cycles=D['fine_cycles']
    fibre_vectors=[]
    for mask in range(256):
        zz=np.zeros(41,dtype=np.int16)
        for j,orbit in enumerate(fine_cycles):
            if mask>>j&1:zz[np.array(orbit)+1]=1
        fibre_vectors.append(zz)
    fibre_vectors=np.array(fibre_vectors);fibre_sums=fibre_vectors@C[:41,:41].T
    # Dual binary character orbits use B^{-T}, computed from the character table.
    H=D['H'];dual=[]
    for a in range(16):
        matches=[b for b in range(16) if np.array_equal(H[a,D['binary']],H[b])]
        assert len(matches)==1;dual.append(matches[0])
    dual_cycles=[c for c in cycles(dual) if c!=[0]]
    records=[]
    table=list(csv.DictReader((ROOT/'work/sweep_20261115/brouwer_live.csv').open()))
    for row in table:
        if int(row['v'])!=1296 or row['status']!='?' or int(row['k'])>647:continue
        k,l,mu=[int(row[x]) for x in ['k','lam','mu']];a=l-mu;disc=a*a+4*(k-mu);gap=math.isqrt(disc)
        if gap*gap!=disc:continue
        r,s=(a+gap)//2,(a-gap)//2
        trace=[(d,j) for j in [0,5,10,15] for d in range(0,81,10) if 16*d==k+15*s+gap*j]
        candidates=[]
        for d,j in trace:
            for mask,sums in enumerate(fibre_sums):
                if sums[0]!=d:continue
                ranks=16*(sums[1:]-s)
                if np.any(ranks%gap) or np.any(ranks<0) or np.any(ranks>16*gap):continue
                quotients=[]
                for bits in range(8):
                    ev=np.full(16,s,dtype=np.int16);ev[0]=k
                    for t,orbit in enumerate(dual_cycles):
                        if bits>>t&1:ev[orbit]=r
                    sizes=H@ev
                    if np.any(sizes%16):continue
                    sizes=sizes//16
                    if sizes[0]!=d or np.any(sizes<0) or np.any(sizes>81):continue
                    quotients.append(dict(positive_dual_orbits=bits,cross_sizes=sizes.tolist()))
                if quotients:candidates.append(dict(fibre_mask=mask,fibre_degree=d,
                    fibre_spectrum=dict(collections.Counter(map(int,sums[1:]))),quotients=quotients))
        records.append(dict(parameters=[1296,k,l,mu],restricted_eigenvalues=[r,s],trace_options=trace,
                            admissible_fibre_count=len(candidates),fibres=candidates))
    # Natural fusions distinguish each binary orbit, membership in D, origin,
    # and the three nonzero quadratic levels in the zero binary fiber.
    labels={};columns=[]
    binary_class={b:i for i,cycle in enumerate(D['binary_cycles']) for b in cycle}
    for b in range(16):
        for j in range(41):
            if b==j==0:continue
            if b==0:
                v=P[j-1];key=('within',int((v@v+int(v.sum())**2)%3))
            else:key=('cross',binary_class[b],j==0,int(z[b,j]))
            if key not in labels:labels[key]=len(labels);columns.append(np.zeros(656,dtype=np.int16))
            columns[labels[key]][41*b+j]=1
    parts=np.column_stack(columns);spectra=C@parts;nparts=len(columns)
    selectors=((np.arange(1<<nparts)[:,None]>>np.arange(nparts))&1).astype(np.int16)
    sums=selectors@spectra.T;fusions=[]
    for mask,ev in enumerate(sums):
        values=set(map(int,ev[1:]))
        if len(values)!=2:continue
        r,s=sorted(values,reverse=True);k=int(ev[0]);mu=k+r*s;l=mu+r+s
        if mu<=0 or l<0 or not 0<k<1295:continue
        params=[1296,k,l,mu]
        status=next((row['status'] for row in table if [int(row[x]) for x in ['v','k','lam','mu']]==params),'absent')
        fusions.append(dict(mask=mask,parameters=params,table_status=status))
    S2=subspaces(2);S3=subspaces(3);assert len(S2)==35 and len(S3)==130
    inside=[];outside=[]
    for i,A in enumerate(S2):
        for j,B in enumerate(S3):
            hits=sum(int(D['T'][sum(x<<(3-t) for t,x in enumerate(a)),D['fi'][b]]) for a in A for b in B)
            if hits==35:inside.append([i,j])
            if hits==0:outside.append([i,j])
    out=dict(scope='Exact C5-and-inversion invariant Cayley constructions on F2^4 x F3^4; not all SRGs',
        fibre_orbits=fine_cycles,binary_difference_orbits=D['binary_cycles'],binary_dual_orbits=dual_cycles,
        open_row_screens=records,natural_fusion_parts=nparts,natural_fusions_checked=1<<nparts,
        strongly_regular_fusions=fusions,order36_subgroups_checked=4550,
        order36_subgroups_inside_D=inside,order36_subgroups_disjoint_from_D=outside,seconds=time.monotonic()-start)
    (HERE/'audit.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ['open_row_screens','strongly_regular_fusions']},indent=2))
    print('Surviving rows:',[(r['parameters'],r['admissible_fibre_count']) for r in records if r['admissible_fibre_count']])
    print('SRG fusions:',fusions)


if __name__=='__main__':main()
