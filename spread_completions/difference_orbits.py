"""Classify maximum cocliques by the earliest group difference, modulo affine symmetries."""
import json
import igraph as ig
import numpy as np
from build import HERE,V,ADD3,NEG3,coset,elements,graph,kernel,sources,subgroups

def add(v,w):
    return 81*((v//81)^(w//81))+ADD3[v%81,w%81]
def negative(v):
    return (v//81)*81+NEG3[v%81]
def affine_generators(a):
    i,j=np.triu(a,1).nonzero()
    host=ig.Graph(n=1296,edges=list(zip(i.tolist(),j.tolist())))
    gens=[]
    for raw in host.automorphism_group(color=[1]+[0]*1295,sh='fl'):
        g=np.array(raw)
        if all(np.array_equal(g[add(V,t)],add(g,g[t])) for t in (81,162,324,648,1,3,9,27)):
            gens.append(g)
    return gens+[negative(V)]
def orbits(d,gens):
    left=set(np.flatnonzero(1-d.ravel()))-{0}
    classes=[]
    while left:
        seen=set();todo=[min(left)]
        while todo:
            x=todo.pop()
            if x not in seen:
                seen.add(x)
                todo.extend(int(g[x]) for g in gens)
        assert seen<=left
        left-=seen;classes.append(seen)
    return sorted(classes,key=lambda s:(-len(s),min(s)))
def audit(p):
    d=kernel(p);a=graph(d);gens=affine_generators(a);classes=orbits(d,gens)
    active=(1-d).ravel();active[0]=0;seeds=set();hits=[]
    for orbit in classes:
        x=min(orbit)
        # All differences of a coclique whose first orbit is this one
        # still belong to active. Translation and a linear automorphism
        # put a point-pair with that difference at 0 and x.
        h=np.flatnonzero(active&active[add(V,negative(x))])
        if len(h)>=34:
            b=active[add(negative(h)[:,None],h)]
            i,j=np.triu(b,1).nonzero()
            local=ig.Graph(n=len(h),edges=list(zip(i.tolist(),j.tolist())))
            cliques=local.cliques(min=34,max=34)
            hits.append(len(cliques))
            for c in cliques:
                seeds.add(tuple(sorted([0,x]+[int(h[i]) for i in c])))
        else:
            hits.append(0)
        active[list(orbit)]=0
    all_cliques=set();queue=list(seeds)
    while queue:
        c=queue.pop()
        if c in all_cliques:continue
        all_cliques.add(c)
        arr=np.array(c)
        assert len(c)==36 and c[0]==0 and not a[np.ix_(arr,arr)].any()
        queue.extend(tuple(sorted(map(int,g[arr]))) for g in gens)
        queue.extend(tuple(map(int,coset(arr,negative(t)))) for t in arr)
    planes={tuple(elements(z)) for z in subgroups(d,0)}
    assert planes<=all_cliques
    return dict(source=p.name,difference_orbits=len(classes),local_hits=hits,through_zero=len(all_cliques),
                subgroups=len(planes),total=36*len(all_cliques),extra=len(all_cliques-planes))
if __name__=='__main__':
    for p in sources():
        if '560' in p.name:print(json.dumps(audit(p)),flush=True)
