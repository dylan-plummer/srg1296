"""Construct the positive-character dual of each supplied graph, in dot-product coordinates."""
from hashlib import sha256
import numpy as np
from survey import coefficients
from build import HERE, graph, graph6, kernel, sources, check
for path in sources():
    if path.suffix!='.g6': continue
    d=kernel(path).astype(int)
    b=np.array([coefficients(i,2) for i in range(16)])
    t=np.array([coefficients(i,3) for i in range(81)])
    character2=1-2*((b@b.T)%2)
    t_phase=(t@t.T)%3
    assert not np.any(character2@d@((t_phase==1).astype(int)-(t_phase==2)))
    ev=character2@d@(3*(t_phase==0)-1)
    assert sorted(np.unique(ev))==[-34,38,1190]
    dual=(ev==38).astype(np.uint8)
    assert dual.sum()==595
    g=graph(dual)
    check(g,(1296,595,274,272))
    out=HERE/('dual_'+path.stem[-1]+'.g6')
    out.write_bytes(graph6(g))
    print(out.name,sha256(out.read_bytes()).hexdigest())
