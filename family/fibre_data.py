#!/usr/bin/env python3
"""Theorem 2.7 of the paper: the binary fiber data of an LS_m(36) PDS in Z_2^4 x Z_3^4 (1 <= m <= 18).
For S subset of A* = F_2^4 \\ {0} let |T_a|(S) = (9/4)(m + W_S(a)) (Lemma 2.3) and Pi_m(S) = {a != 0 : |T_a|(S) odd}.
Necessary (Lemmas 2.3, 2.5): all |T_a|(S) are integers in [0, 81], U_0 = Pi_m(S_1) and S_1 = Pi_m(U_0).
This script lists all S surviving these conditions for each m, with |T_0| = (5m + 9|S|)/4 and the ring sizes, and
identifies every survivor geometrically: the sets {x != 0 : Q(x) = e} for the quadratic forms Q on F_2^4, and the
punctured planes and their complements.  Plain Python 3.
usage: python3 fibre_data.py > fibre_data.out"""
import itertools, collections
H = [[(-1) ** (bin(a & b).count('1') % 2) for b in range(16)] for a in range(16)]
def sizes(m, S):
    out = []
    for a in range(1, 16):
        n = 9 * (m + sum(H[a][b] for b in S))
        if n % 4 or not 0 <= n // 4 <= 81: return None
        out.append(n // 4)
    return out
def Pi(m, S):
    q = sizes(m, S); return None if q is None else frozenset(a for a in range(1, 16) if q[a - 1] % 2)
# geometric catalog
bit = lambda x, i: (x >> i) & 1
shapes = collections.defaultdict(set)
for coef in itertools.product(range(2), repeat=10):           # Q(x) = sum_{i<=j} q_ij x_i x_j
    idx = [(i, j) for i in range(4) for j in range(i, 4)]
    Q = lambda x: sum(c * bit(x, i) * bit(x, j) for c, (i, j) in zip(coef, idx)) % 2
    zeros = frozenset(x for x in range(1, 16) if Q(x) == 0); ones = frozenset(x for x in range(1, 16) if Q(x) == 1)
    nz = len(zeros) + 1                                          # number of zeros including 0
    # nondegenerate forms on F_2^4: elliptic (6 zeros) / hyperbolic (10 zeros); check nondegeneracy via polar form
    B = [[(Q(1 << i | 1 << j) - Q(1 << i) - Q(1 << j)) % 2 if i != j else 0 for j in range(4)] for i in range(4)]
    rows = [sum(B[i][j] << j for j in range(4)) for i in range(4)]; rank = 0
    for col in range(4):
        piv = next((r for r in range(rank, 4) if rows[r] >> col & 1), None)
        if piv is None: continue
        rows[rank], rows[piv] = rows[piv], rows[rank]
        for r in range(4):
            if r != rank and rows[r] >> col & 1: rows[r] ^= rows[rank]
        rank += 1
    det = int(rank == 4)
    if det == 1 and nz == 6: shapes['elliptic: Q=0 (ovoid, 5)'].add(zeros); shapes['elliptic: Q=1 (10)'].add(ones)
    if det == 1 and nz == 10: shapes['hyperbolic: Q=0 (3x3 grid, 9)'].add(zeros); shapes['hyperbolic: Q=1 (6)'].add(ones)
planes = {frozenset({a, b, a ^ b}) for a in range(1, 16) for b in range(1, 16) if a != b}
shapes['punctured plane (3)'] = planes
shapes['complement of a plane (12)'] = {frozenset(range(1, 16)) - p for p in planes}
shapes['empty'] = {frozenset()}; shapes['all of A* (15)'] = {frozenset(range(1, 16))}
def shape(S): return [n for n, F in shapes.items() if S in F]
for m in range(1, 19):
    surv = collections.defaultdict(list)
    for j in range(16):
        for S in itertools.combinations(range(1, 16), j):
            S = frozenset(S); U = Pi(m, S)
            if U is None or Pi(m, U) != S: continue
            surv[(len(S), (5 * m + 9 * len(S)) // 4, tuple(sorted(collections.Counter(sizes(m, S)).items())))].append((S, U))
    for (j, d, sz), L in surv.items():
        sh_S = collections.Counter(tuple(shape(S)) for S, U in L); sh_U = collections.Counter(tuple(shape(U)) for S, U in L)
        allshape = all(len(shape(S)) == 1 for S, U in L) and len(sh_S) == 1 and len(L) == len(shapes[shape(L[0][0])[0]])
        print(f'LS_{m:<2d}: |S_1| = |U_0| = {j:2d}, |T_0| = {d:2d}, ring sizes {dict(sz)}; {len(L)} choices of S_1, '
              f'all of shape {list(sh_S)[0]} (all sets of that shape: {allshape}); U_0 shape {list(sh_U)}')
