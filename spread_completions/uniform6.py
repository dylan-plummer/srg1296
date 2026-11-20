"""Compare the five parallel resolutions of the other extensible degree-560 host."""
from build import HERE, V, check, elements, graph, graph6, kernel, partition, sources, subgroups
from compare import canon, encoding
from hashlib import sha256
import json
import igraph as ig
import numpy as np
if __name__=='__main__':
    source=next(p for p in sources() if p.name=='srg1296_560_244_240_class6.txt')
    d=kernel(source);u=subgroups(d,0);classes={}
    for z in u:
        a=graph(d)
        for c in partition(elements(z),V):
            a[np.ix_(c,c)]=1
        np.fill_diagonal(a,0)
        check(a,(1296,595,274,272))
        i,j=np.triu(a,1).nonzero()
        g=ig.Graph(n=len(a),edges=list(zip(i.tolist(),j.tolist())))
        h=sha256(encoding(canon(g))).hexdigest()
        if h not in classes:
            name='uniform6_%d.g6'%(1+len(classes))
            classes[h]=name
            (HERE/name).write_bytes(graph6(a))
        print(json.dumps(dict(subgroup=z,canonical_sha256=h,file=classes[h])),flush=True)
