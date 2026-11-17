#!/usr/bin/env python3
"""Graphs VIII.1, VIII.2 (LS_11) each have one free subgroup U of order 36 with D u U* an LS_12 PDS (Proposition 6.2).
Find U, build E = D u U*, verify E by difference counting, and compare Cay(G,E) with VII.1 and between the two (nauty).
usage: python3 spread_link.py > spread_link.out   (needs numpy, pynauty)"""
import itertools, collections, sys, os
sys.argv = [sys.argv[0]]
__file__ = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'nauty_check.py')
exec(open(__file__).read().split('roman =')[0])
from fields import *
F2 = [tuple((t >> s) & 1 for s in range(4)) for t in range(16)]; F3 = [tuple((t // 3 ** s) % 3 for s in range(4)) for t in range(81)]
def span(gens, p):
    S = {(0,) * 4}
    for g in gens: S = {tuple((x + c * y) % p for x, y in zip(s, g)) for s in S for c in range(p)}
    return frozenset(S)
sub2 = {span([a, b], 2) for a, b in itertools.combinations(F2[1:], 2)}; sub2 = [s for s in sub2 if len(s) == 4]
sub3 = {span([a, b], 3) for a, b in itertools.combinations(F3[1:], 2)}; sub3 = [s for s in sub3 if len(s) == 9]
assert len(sub2) == 35 and len(sub3) == 130
def pds_check(E, k, lam, mu):
    cnt = collections.Counter((tuple((p - q) % 2 for p, q in zip(x[0], y[0])), tuple((p - q) % 3 for p, q in zip(x[1], y[1]))) for x in E for y in E)
    return all(cnt[(a, b)] == (k if (a == Z4 and b == Z4) else lam if (a, b) in E else mu) for a in F2 for b in F3)
_, D7 = read_connset('../connsets/srg1296_420_144_132_class1.txt'); c7 = pynauty.certificate(graph(D7))
Es = []
for cls in (1, 2):
    _, D = read_connset(f'../connsets/srg1296_385_124_110_class{cls}.txt')
    free = [(U2, U3) for U2 in sub2 for U3 in sub3 if not any((a, b) in D for a in U2 for b in U3)]
    print(f'VIII.{cls}: free subgroups of order 36: {len(free)}')
    for U2, U3 in free:
        E = D | {(a, b) for a in U2 for b in U3 if (a, b) != (Z4, Z4)}
        ok = len(E) == 420 and pds_check(E, 420, 144, 132)
        c = pynauty.certificate(graph(E)); aut = pynauty.autgrp(graph(E))[1]
        print(f'  E = D u U*: |E| = {len(E)}, LS_12 PDS (difference count): {ok}, |Aut| = {aut:.0f}, isomorphic to VII.1: {c == c7}')
        Es.append((cls, E, c))
        write_connset(f'/mnt/jinstore/JinLab02/dmp131/paut_runs/E_from_VIII{cls}.txt', (1296, 420, 144, 132), E, f'D(VIII.{cls}) u U*')
print('E from VIII.1 isomorphic to E from VIII.2:', Es[0][2] == Es[1][2])
