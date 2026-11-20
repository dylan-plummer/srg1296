"""Compare graph6 files by canonical forms; independently count common neighbours."""
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import sys
import igraph as ig
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'fibre_switch_20261120'))
from verify_graph6 import read_graph,local_cliques
from construct import graph6

def canon(g):
    return g.permute_vertices(g.canonical_permutation(sh='fl'))
def encoding(g):
    edges=np.array(g.get_edgelist(),dtype=int)
    a=np.zeros((g.vcount(),g.vcount()),dtype=np.uint8)
    i,j=edges.T
    a[i,j]=a[j,i]=1
    return graph6(a)
def analyze(p):
    rows,digest=read_graph(p)
    n=len(rows)
    size=Counter(r.bit_count() for r in rows)
    assert len(size)==1
    adjacent=Counter()
    separate=Counter()
    edges=[]
    for j,r in enumerate(rows):
        for i in range(j):
            if (r>>i)&1:
                adjacent[(r&rows[i]).bit_count()]+=1
                edges.append((i,j))
            else:
                separate[(r&rows[i]).bit_count()]+=1
    assert len(adjacent)==len(separate)==1
    g=ig.Graph(n=n,edges=edges)
    canonical=encoding(canon(g))
    parent=list(range(n))
    def find(i):
        while parent[i]!=i:
            parent[i]=parent[parent[i]]
            i=parent[i]
        return i
    for perm in g.automorphism_group(sh='fl'):
        for i,t in enumerate(perm):
            parent[find(i)]=find(t)
    orbits=Counter(find(i) for i in range(n))
    loc=Counter()
    for rep,sz in orbits.items():
        loc[local_cliques(rows,rep)]+=sz
    return dict(file=p.name,sha256=digest,parameters=[n,next(iter(size)),next(iter(adjacent)),next(iter(separate))],
                bliss_sha256=sha256(canonical).hexdigest(),automorphisms=g.count_automorphisms(sh='fl'),
                orbits=sorted(orbits.values()),K4_per_vertex=sorted(loc.items()))
if __name__=='__main__':
    control=ig.Graph(n=7,edges=[(0,1),(0,2),(1,2),(2,3),(3,4),(3,5),(4,6)])
    assert encoding(canon(control))==encoding(canon(control.permute_vertices([2,5,0,4,1,6,3])))
    for name in sys.argv[1:]:
        print(json.dumps(analyze(Path(name))),flush=True)
