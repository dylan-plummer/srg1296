#!/usr/bin/env python3
"""For every connection set: subgroups U of order 36 with U* inside D (cliques; deleting U gives LS_{m-1} by Kovaľ's
lemma) and free subgroups (U n D empty).  Plain Python 3.  usage: python3 spread_inside.py > spread_inside.out"""
import glob, os, itertools, collections
from fields import *
F2 = [tuple((t >> s) & 1 for s in range(4)) for t in range(16)]; F3 = [tuple((t // 3 ** s) % 3 for s in range(4)) for t in range(81)]
def span(gens, p):
    S = {(0,) * 4}
    for g in gens: S = {tuple((x + c * y) % p for x, y in zip(s, g)) for s in S for c in range(p)}
    return frozenset(S)
sub2 = [s for s in {span([a, b], 2) for a, b in itertools.combinations(F2[1:], 2)} if len(s) == 4]
sub3 = [s for s in {span([a, b], 3) for a, b in itertools.combinations(F3[1:], 2)} if len(s) == 9]
roman = {560: 'I', 592: 'II', 595: 'III', 525: 'IV', 490: 'V', 455: 'VI', 420: 'VII', 385: 'VIII'}
CONN = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'connsets')
for fn in sorted(glob.glob(os.path.join(CONN, 'srg1296_*_class*.txt')), key=lambda f: (-int(f.split('_')[1]), int(f.split('class')[1][:-4]))):
    (v, k, lam, mu), D = read_connset(fn); name = f'{roman[k]}.{int(fn.split("class")[1][:-4])}'
    inside = free = 0
    for U2 in sub2:
        for U3 in sub3:
            U = [(a, b) for a in U2 for b in U3 if (a, b) != (Z4, Z4)]
            n = sum(x in D for x in U)
            inside += n == 35; free += n == 0
    print(f'{name:6s} k={k}: subgroups of order 36 inside D: {inside}, free: {free}', flush=True)
