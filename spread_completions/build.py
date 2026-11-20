"""Trade two coset resolutions inside selected translates of their sum."""
import argparse
from hashlib import sha256
import json
from itertools import combinations
import numpy as np
from survey import HERE, ROOT, ADD3, kernel, meet, sources, subgroups
from construct import check, graph6

V=np.arange(1296)
NEG3=np.array([np.flatnonzero(a==0)[0] for a in ADD3])
def coset(points, offset):
    points=np.array(points,dtype=int)
    return np.sort(81*((points//81)^(offset//81))+ADD3[points%81,offset%81])
def elements(u):
    return np.array([81*a+b for a in u[0] for b in u[1]])
def partition(group,universe):
    left=set(universe);out=[]
    while left:
        block=coset(group,min(left))
        left.difference_update(block)
        out.append(block)
    return out
def graph(d):
    return d[(V//81)[:,None]^(V//81),ADD3[NEG3[V%81,None],V%81]]
def resolutions(u,v):
    a,b=elements(u),elements(v)
    w=np.unique(np.concatenate([coset(a,t) for t in b]))
    return [(partition(a,s),partition(b,s)) for s in partition(w,V)]
def complete(a,choices,cut):
    a=a.copy()
    for i,blocks in enumerate(choices):
        for block in blocks[(cut>>i)&1]:
            sub=np.ix_(block,block)
            assert not a[sub].any()
            a[sub]=1
    np.fill_diagonal(a,0)
    return a
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('source')
    ap.add_argument('intersection',type=int)
    ap.add_argument('cut',type=lambda x:int(x,0))
    ap.add_argument('out')
    args=ap.parse_args()
    p=next(p for p in sources() if p.name==args.source)
    d=kernel(p);us=subgroups(d,0)
    u,v=next((u,v) for u,v in combinations(us,2) if meet(u,v)==args.intersection)
    choices=resolutions(u,v)
    assert 0<args.cut<2**len(choices)-1
    a=complete(graph(d),choices,args.cut)
    m=int(d.sum())//35+1
    param=(1296,m*35,36-2+(m-1)*(m-2),m*(m-1))
    print(json.dumps(dict(source=p.name,subgroups=[u,v],intersection=args.intersection,cosets=len(choices),
                          cut=args.cut,parameters=param,pair_counts=check(a,param))),flush=True)
    out=HERE/args.out
    out.write_bytes(graph6(a))
    print(sha256(out.read_bytes()).hexdigest(),out.name)
if __name__=='__main__':main()
