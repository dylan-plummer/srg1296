#!/usr/bin/env python3
"""Graphs VII.2, VII.3 contain a subgroup U of order 36 (a spread of Delsarte cliques, omega = 36).  An OA(36,12) graph
has every edge in a Delsarte clique; here we exhibit edges {0, x} lying in no 36-clique, i.e. with omega(N(0) & N(x)) < 34.
usage: python3 oa_test.py > oa_test.out   (needs numpy, igraph)"""
import sys, numpy as np, igraph as ig
from fields import *
F2 = [tuple((t >> s) & 1 for s in range(4)) for t in range(16)]; F3 = [tuple((t // 3 ** s) % 3 for s in range(4)) for t in range(81)]
i2 = {a: i for i, a in enumerate(F2)}; i3 = {b: i for i, b in enumerate(F3)}
d2 = np.array([[i2[tuple((p - q) % 2 for p, q in zip(x, y))] for y in F2] for x in F2]); d3 = np.array([[i3[tuple((p - q) % 3 for p, q in zip(x, y))] for y in F3] for x in F3])
ea, eb = np.divmod(np.arange(1296), 81)
for c in (2, 3):
    _, D = read_connset(f'../connsets/srg1296_420_144_132_class{c}.txt')
    M = np.zeros((16, 81), bool)
    for a, b in D: M[i2[a], i3[b]] = True
    A = M[d2[ea][:, ea], d3[eb][:, eb]]
    nb0 = np.flatnonzero(A[0]); res = []
    for x in nb0[:60]:                                         # first 60 neighbors of 0 (all orbits are not needed)
        cn = np.flatnonzero(A[0] & A[x]); g = ig.Graph.Adjacency(A[np.ix_(cn, cn)].tolist(), mode='undirected')
        res.append(g.clique_number() + 2)
    print(f'VII.{c}: clique numbers of edges {{0,x}} for 60 neighbours x: min {min(res)}, max {max(res)}; '
          f'edges in no Delsarte clique: {sum(r < 36 for r in res)} of 60 -> not an OA graph: {min(res) < 36}', flush=True)
