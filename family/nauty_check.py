#!/usr/bin/env python3
"""nauty (pynauty) confirmation for Table 6 and Proposition 2.5: |Aut| of every graph in ../connsets, pairwise
non-isomorphism within each parameter set, and which duals Cay(G, D*) are isomorphic to Cay(G, D).
usage: python3 nauty_check.py > nauty_check.out   (needs numpy, pynauty)"""
import os, glob, collections
import numpy as np, pynauty
from fields import *
from duality import dual_set
HERE = os.path.dirname(os.path.abspath(__file__)); CONN = os.path.join(HERE, '..', 'connsets')
F2 = [tuple((t >> s) & 1 for s in range(4)) for t in range(16)]; F3 = [tuple((t // 3 ** s) % 3 for s in range(4)) for t in range(81)]
i2 = {a: i for i, a in enumerate(F2)}; i3 = {b: i for i, b in enumerate(F3)}
d2 = np.array([[i2[tuple((p - q) % 2 for p, q in zip(x, y))] for y in F2] for x in F2])
d3 = np.array([[i3[tuple((p - q) % 3 for p, q in zip(x, y))] for y in F3] for x in F3])
ea, eb = np.divmod(np.arange(1296), 81)
def graph(D):
    M = np.zeros((16, 81), bool)
    for a, b in D: M[i2[a], i3[b]] = True
    A = M[d2[ea][:, ea], d3[eb][:, eb]]
    return pynauty.Graph(1296, adjacency_dict={v: list(np.flatnonzero(A[v])) for v in range(1296)})
roman = {560: 'I', 592: 'II', 595: 'III', 525: 'IV', 490: 'V', 455: 'VI', 420: 'VII', 385: 'VIII'}
certs = collections.defaultdict(dict)
for fn in sorted(glob.glob(os.path.join(CONN, 'srg1296_*_class*.txt')), key=lambda f: (-int(f.split('_')[1]), int(f.split('class')[1][:-4]))):
    (v, k, lam, mu), D = read_connset(fn); name = f'{roman[k]}.{int(fn.split("class")[1][:-4])}'
    g = graph(D); gen, order, _, _, _ = pynauty.autgrp(g); aut = int(round(order))
    c = pynauty.certificate(g); cd = pynauty.certificate(graph(dual_set(D, k)))
    certs[k][name] = (c, cd)
    print(f'{name:6s} |Aut| = {aut}', flush=True)
for k, d in certs.items():
    names = list(d); distinct = len({d[n][0] for n in names}) == len(names)
    selfdual = [n for n in names if d[n][0] == d[n][1]]
    dual_of = {n: [m for m in names if d[m][0] == d[n][1]] for n in names}
    print(f'k={k}: {len(names)} graphs, pairwise non-isomorphic: {distinct}; Cay(D*) = Cay(D) for {selfdual};'
          f' dual classes: { {n: v for n, v in dual_of.items() if n not in selfdual} }', flush=True)
