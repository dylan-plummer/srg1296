#!/usr/bin/env python3
"""Shared local structure of the Latin square type graphs I, III-VIII (and II for contrast).  Plain Python 3.
For every graph and every ring a != 0 the ring set T_a is classified by its *local type*: |T_a|, whether 0 in T_a, and the
multiset of the 40 values psi_d(T_a)/9 (d over projective points; integers by Lemma 2.3).  Also reported:
the character values of T_0; the F_2-fibres U_b (b != 0) by weight and by whether they are affine subspaces;
the sizes |S_psi|; and 'matching': the share of characters psi with |S_psi| = 8 whose row (psi(T_a))_{a != 0} has a
single nonzero entry (one ring carries the whole positive Walsh pattern).
usage: python3 shared_structure.py > shared_structure.out"""
import os, glob, collections, math
from fields import *
from duality import TRA, TRB
HERE = os.path.dirname(os.path.abspath(__file__)); CONN = os.path.join(HERE, '..', 'connsets')
P = [b for b in B81 if logB[b] < 40]
def psi(d, T):  # character sum of a symmetric set T (contains b and -b) at d
    return sum({0: 2, 1: -1, 2: -1}[TRB[(d, x)]] for x in T if x in Pset) // 2 + (1 if Z4 in T else 0) * 0
Pset = set(B81)
def chi_sum(d, T):
    s2 = sum({0: 2, 1: -1, 2: -1}[TRB[(d, x)]] for x in T if x != Z4)
    return s2 // 2 + (1 if Z4 in T else 0)
def is_affine(U):
    U = set(U)
    if not U: return False
    t = next(iter(U)); V = {tuple((p + q) % 2 for p, q in zip(u, t)) for u in U}
    return all(tuple((p + q) % 2 for p, q in zip(x, y)) in V for x in V for y in V)
roman = {560: 'I', 592: 'II', 595: 'III', 525: 'IV', 490: 'V', 455: 'VI', 420: 'VII', 385: 'VIII'}
types = collections.defaultdict(set); lines = []
files = sorted(glob.glob(os.path.join(CONN, 'srg1296_*_class*.txt')), key=lambda f: (-int(f.split('_')[1]), int(f.split('class')[1][:-4])))
for fn in files:
    (v, k, lam, mu), D = read_connset(fn); name = f'{roman[k]}.{int(fn.split("class")[1][:-4])}'
    g = math.isqrt((lam - mu) ** 2 + 4 * (k - mu)); r, s = (lam - mu + g) // 2, (lam - mu - g) // 2
    T = collections.defaultdict(set); U = collections.defaultdict(set)
    for a, b in D: T[a].add(b); U[b].add(a)
    loc = collections.Counter()
    for a in A16:
        vals = collections.Counter(chi_sum(d, T[a]) // 9 for d in P)
        assert all(chi_sum(d, T[a]) % 9 == 0 for d in P)
        key = (len(T[a]), Z4 in T[a], tuple(sorted(vals.items())))
        loc[key] += 1; types[key].add(roman[k])
    t0 = collections.Counter(chi_sum(d, T[Z4]) for d in P)
    ub = collections.Counter((len(U[b]), is_affine(U[b])) for b in P)
    # S_psi sizes and matching
    Sz = collections.Counter(); match = [0, 0]
    for d in P:
        row = {a: chi_sum(d, T[a]) for a in [Z4] + A16}
        S = 0
        for c in [Z4] + A16:
            val = sum((-1) ** TRA[(c, a)] * row[a] for a in row)
            S += (val == r)
        Sz[S] += 1
        if S == 8:
            match[1] += 1; match[0] += (sum(1 for a in A16 if row[a] != 0) == 1)
    lines.append((name, k, dict(loc), dict(sorted(t0.items())), dict(sorted(ub.items())), dict(sorted(Sz.items())), match))
names = {key: f'L{i+1}' for i, key in enumerate(sorted(types, key=lambda x: (x[0], x[1], x[2])))}
print('Local types of ring sets T_a (a != 0): size, 0 in T, multiset of psi(T_a)/9 over the 40 projective psi')
for key, nm in names.items():
    print(f'  {nm:4s} |T|={key[0]:2d} 0in={str(key[1]):5s} values {dict(key[2])}   used by graphs: {sorted(types[key])}')
print()
for name, k, loc, t0, ub, Sz, match in lines:
    print(f'{name:6s} k={k}: rings {dict((names[x], n) for x, n in loc.items())}; T0 character sums {t0}; '
          f'U_b (weight, affine) {ub}; |S_psi| {Sz}; matching {match[0]}/{match[1]}')

# geometric identification of some local types (hyperplanes W = ker of a nonzero trace form b -> Tr(d b))
hyper = {d: frozenset([Z4] + [b for b in B81 if TRB[(d, b)] == 0]) for d in P}          # 27-element subspaces
def ident(Tset):
    T = frozenset(Tset); F = frozenset([Z4] + B81); out = []
    for d1, d2 in __import__('itertools').combinations(P, 2):
        W1, W2 = hyper[d1], hyper[d2]
        if T == F - (W1 | W2): out.append('F81 minus two hyperplanes')
        if T == (W1 | W2): out.append('union of two hyperplanes')
    for d in P:   # affine plane: (W minus a 2-dim subspace V) with V = W & W'
        for d2 in P:
            if d2 == d: continue
            if T == hyper[d] - (hyper[d] & hyper[d2]): out.append('hyperplane minus a 2-dim subspace (an affine plane of PG(3,3), +-)')
    return sorted(set(out))
print()
print('Geometric identification of local types (all rings of all graphs):')
seen = collections.defaultdict(collections.Counter)
for fn in files:
    (v, k, lam, mu), D = read_connset(fn)
    T = collections.defaultdict(set)
    for a, b in D: T[a].add(b)
    for a in A16:
        vals = collections.Counter(chi_sum(d, T[a]) // 9 for d in P)
        key = (len(T[a]), Z4 in T[a], tuple(sorted(vals.items())))
        seen[names[key]][tuple(ident(T[a]))] += 1
for nm in sorted(seen, key=lambda x: int(x[1:])):
    print(f'  {nm}: {dict(seen[nm])}')
