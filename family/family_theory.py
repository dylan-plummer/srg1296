#!/usr/bin/env python3
"""Finite checks behind Proposition 5.4 of the paper (connection sets invariant under rho = (alpha^12, beta^16)).
Plain Python 3, no packages.
 1. Trace conditions (Lemma 6.2) in the rho-model: |T_0| in 10Z, j in 5Z, |U_0| in 5Z, j' in 10Z. List every
    LS_m(36) and NL_m(36) parameter set that passes.
 2. Gauss periods of index 8 of F_81 (sums of psi over a class beta^i <beta^8>), and for every union T_0 of g classes:
    the congruence psi(T_0) in {s + 9t} for all psi != 1 (Lemma 6.1 at a = 0) versus strong regularity of Cay(F_81,T_0).
 3. Ring sizes |T_a| = (9/4)(m + W(a)), W the Walsh sum of S_1 (a union of j/5 ovoids mu_5-cosets in the dual).
usage: python3 family_theory.py > family_theory.out"""
import itertools
from collections import Counter
from fields import *
from duality import TRB, TRA, ELA

print('1. trace conditions in the rho-model')
passing = {}
for typ in ('LS', 'NL'):
    for m in range(1, 36):
        k, r, s = (35 * m, 36 - m, -m) if typ == 'LS' else (37 * m, m, m - 36)
        if not 0 < k <= 648: continue
        sol = [(d, j, dp, jp) for j in range(0, 16, 5) for jp in range(0, 81, 10)
               for d in [(k + 15 * s + 36 * j) / 16] for dp in [(k + 80 * s + 36 * jp) / 81]
               if d == int(d) and dp == int(dp) and int(d) % 10 == 0 and 0 <= d <= 80 and int(dp) % 5 == 0 and 0 <= dp <= 15]
        if sol:
            passing[(typ, m)] = [tuple(map(int, x)) for x in sol]; print(f'  {typ}_{m}(36), k = {k}: (|T0|, j, |U0|, j\') = {passing[(typ, m)]}')
print('  (all other LS_m / NL_m with k <= 648 fail; complements follow)')

print('2. Gauss periods of index 8 and unions of classes')
per = {i: {d: sum({0: 2, 1: -1, 2: -1}[TRB[(d, b)]] for b in B81 if logB[b] % 8 == i) // 2 for d in B81} for i in range(8)}
for i in range(8): print(f'  class {i}: values', dict(sorted(Counter(per[i].values()).items())))
sub = lambda x, y: tuple((p - q) % 3 for p, q in zip(x, y))
def is_srg(S):
    c = Counter(sub(x, y) for x in S for y in S if x != y)
    return len({c[x] for x in S}) == 1 and len({c[t] for t in B81 if t not in S}) == 1
for (typ, m), sol in passing.items():
    d = sol[0][0]; g = d // 10; s = -m if typ == 'LS' else m - 36
    allowed = {s + 9 * t for t in range(5)}
    ok = [J for J in itertools.combinations(range(8), g) if all(sum(per[i][x] for i in J) in allowed for x in B81)]
    srg = [J for J in ok if is_srg({b for b in B81 if logB[b] % 8 in J})]
    print(f'  {typ}_{m}: g = {g} classes; unions satisfying the congruence: {len(ok)}, all strongly regular: {len(srg) == len(ok)}'
          + ('' if ok else '  -> NO rho-invariant PDS'))

print('3. ring sizes from the Walsh sums of S_1')
mu5 = [A16[3 * t] for t in range(5)]
W = lambda S, a: sum((-1) ** TRA[(c, a)] for c in S)
for (typ, m), sol in passing.items():
    if typ != 'LS': continue
    j = sol[0][1]; d = sol[0][0]
    for cosets in itertools.combinations(range(3), j // 5):
        S1 = [mulA(u, x) for u in cosets for x in mu5]        # union of j/5 cosets of mu_5 (in the dual, via the trace form)
        sizes = Counter(9 * (m + W(S1, a)) / 4 for a in A16)
        print(f'  LS_{m}: S_1 = {j // 5} coset(s) of mu_5: |T_a| (a != 0) takes values {dict(sorted(sizes.items()))}')
        break
