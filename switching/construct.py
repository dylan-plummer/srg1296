"""Develop and independently check a switch of the 16-point binary fibres."""
from collections import Counter
from itertools import combinations
from hashlib import sha256
from pathlib import Path
import argparse,json
import numpy as np
from scan import ROOT,X,Y,H,signs,bilinear,data

def build_cut(path,side,w):
    param,masks=data(path)
    blocks=[]
    for m in masks.values():
        row=np.array([(m>>j)&1 for j in X])
        switched=(row@H*signs[w])@H
        blocks.append((row,switched))
    diffs=[[Y.index(tuple((a-b)%3 for a,b in zip(z,y))) for z in Y] for y in Y]
    A=np.empty((1296,1296),dtype=np.uint8);xor=np.bitwise_xor.outer(X,X)
    for i in range(81):
        for j in range(81):
            a,b=blocks[diffs[i][j]]
            if side[i]!=side[j]:
                assert np.all((b==0)|(b==16));a=b//16
            A[16*i:16*i+16,16*j:16*j+16]=a[xor]
    return param,A

def construct(path,u,w):
    return build_cut(path,[bilinear(u,y)==0 for y in Y],w)

def check(A,param):
    v,k,lam,mu=param
    assert A.shape==(v,v) and (A==A.T).all() and not np.any(np.diag(A))
    masks=[int.from_bytes(r,'little') for r in np.packbits(A,axis=1,bitorder='little')]
    assert all(row.bit_count()==k for row in masks)
    observed=Counter()
    for j,r in enumerate(masks):
        for i in range(j):observed[(r>>i&1,(r&masks[i]).bit_count())]+=1
    assert observed=={(1,lam):v*k//2,(0,mu):v*(v-k-1)//2}
    return {str(c):n for c,n in sorted(observed.items())}

def invariant(A):
    import igraph as ig
    i,j=np.triu(A,1).nonzero();g=ig.Graph(n=len(A),edges=list(zip(i.tolist(),j.tolist())))
    gens=g.automorphism_group();p=list(range(len(A)))
    def root(i):
        while p[i]!=i:
            p[i]=p[p[i]];i=p[i]
        return i
    for perm in gens:
        for i,j in enumerate(perm):p[root(i)]=root(j)
    sizes=sorted(Counter(root(i) for i in range(len(A))).values())
    return {'automorphism_order':g.count_automorphisms(),'orbit_sizes':sizes}

def graph6(A):
    n=len(A);i,j=np.tril_indices(n,-1)
    data=np.concatenate([A[i,j],np.zeros((-len(i))%6,dtype=np.uint8)]).reshape(-1,6)
    encoded=data@np.array([32,16,8,4,2,1],dtype=np.uint8)
    return bytes([126,63+(n>>12),63+((n>>6)&63),63+(n&63)])+(encoded+63).tobytes()+b'\n'

def main():
    ap=argparse.ArgumentParser();ap.add_argument('degree',type=int);ap.add_argument('limit',type=int);args=ap.parse_args()
    runs=[json.loads(l) for l in (ROOT/'scan.jsonl').read_text().splitlines()]
    saved=False;count=0
    for r in runs:
        if r['parameters'][1]!=args.degree:continue
        path=ROOT.parent/'connsets'/r['file']
        for u,ws in r['examples']:
            for w in ws:
                param,A=construct(path,u,w)
                summary=dict(source=path.name,normal=u,gauge_index=w,negative_frequencies=np.flatnonzero(signs[w]<0).tolist(),
                             parameters=param,pair_census=check(A,param),**invariant(A))
                print(json.dumps(summary),flush=True)
                if len(summary['orbit_sizes'])>1:
                    out=ROOT/f'switched_{args.degree}.g6';b=graph6(A);out.write_bytes(b)
                    summary.update(sha256=sha256(b).hexdigest(),graph=out.name)
                    (ROOT/f'switched_{args.degree}.json').write_text(json.dumps(summary,indent=2)+'\n')
                    return
                count+=1
                if count>=args.limit:return

if __name__=='__main__':main()
