"""Enumerate every 36-coclique through zero, using no subgroup shape assumption."""
import json
from itertools import combinations
from collections import Counter
import igraph as ig
import numpy as np
from build import HERE, V, coset, elements, graph, kernel, sources, subgroups

def enumerate_zero(a):
    complement=1-a
    np.fill_diagonal(complement,0)
    out=[]
    for x in np.flatnonzero(complement[0]):
        # Uniquely assigns each clique to its smallest vertex after zero.
        h=np.flatnonzero(complement[0]&complement[x]&(V>x))
        if len(h)<34: continue
        i,j=np.triu(complement[np.ix_(h,h)],1).nonzero()
        g=ig.Graph(n=len(h),edges=list(zip(i.tolist(),j.tolist())))
        for cl in g.cliques(min=34,max=34):
            block=(0,int(x),*sorted(int(h[y]) for y in cl))
            assert len(set(block))==36 and not a[np.ix_(block,block)].any()
            out.append(block)
    assert len(out)==len(set(out))
    return out
if __name__=='__main__':
    for p in sources():
        if '560' not in p.name:continue
        d=kernel(p);a=graph(d)
        candidates=enumerate_zero(a)
        planes={tuple(elements(z)) for z in subgroups(d,0)}
        assert planes.issubset(candidates)
        print(json.dumps(dict(source=p.name,through_zero=len(candidates),subgroups=len(planes),
                         total=36*len(candidates),non_subgroup_examples=[c for c in candidates if c not in planes][:2])),flush=True)
        if candidates:
            (HERE/(p.stem+'_cocliques.json')).write_text(json.dumps(candidates)+'\n')
