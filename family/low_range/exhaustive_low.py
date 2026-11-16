#!/usr/bin/env python3
"""Solver-free exhaustive search for rho-invariant LS_m(36) PDSs with m in {5,6,7} (Proposition 5.4, low range).
By Proposition 5.4 (and multiplication by F16* x F81*, which commutes with rho) we may fix
  T_0 = union of the index-8 classes J0 (an NL_g(9) set, g = 9 - m),  U_0 = union of the mu_5-cosets I0 (e = g - 1 cosets),
  T_(alpha^u) = {0} u (+-line of PG(3,3))  for u in I0   (|T| = 9: a 9-divisible 4-point set is a line),
  T_(alpha^u) = +-P, |P| = 9, all plane sections of size 0 mod 3   for u not in I0  (|T| = 18),
and T_(alpha^(3t+u)) = beta^(-16t) T_(alpha^u).  For each class representative d = beta^i (i < 8) and each c in F16 the
character sum  sum_a (-1)^Tr(ca) psi_d(T_a)  must lie in {r, s}.  Where psi_d(T_0) gives |S_psi| = 0, all 16 sums equal s;
these rows are matched exactly by a hash join over the last seed; all other rows and the full PDS property are then checked.
usage: exhaustive_low.py m J0 I0     e.g.  exhaustive_low.py 7 0,4 0"""
import sys, itertools, collections, time
import numpy as np
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '..'))
from fields import *
from duality import TRA, TRB
m = int(sys.argv[1]); J0 = [int(x) for x in sys.argv[2].split(',')]; I0 = [int(x) for x in sys.argv[3].split(',')]
r, s = 36 - m, -m; t0 = time.time()
P = [b for b in B81 if logB[b] < 40]                       # 40 projective points: beta^j, j < 40 (beta^40 = -1)
D80 = B81                                                  # nonzero d
C = np.array([[2 if TRB[(d, x)] == 0 else -1 for x in P] for d in D80], dtype=np.int64)   # psi_d(+-x), 80 x 40
didx = {d: i for i, d in enumerate(D80)}
# catalogues
mom = lambda S: tuple(sum(x[i] * x[j] for x in S) % 3 for i in range(4) for j in range(i, 4))
lines = [S for S in itertools.combinations(range(40), 4) if all((1 + C[:, list(S)].sum(1)) % 9 == 0)]
def nine_sets():
    L, R = range(20), range(20, 40); left = collections.defaultdict(list)
    for w in range(10):
        for S in itertools.combinations(L, w): left[(w, mom([P[i] for i in S]))].append(S)
    out = []
    for w in range(10):
        for S in itertools.combinations(R, w):
            mm = mom([P[i] for i in S]); need = tuple((-x) % 3 for x in mm)
            for T in left.get((9 - w, need), []): out.append(T + S)
    return [S for S in out if all(C[:, list(S)].sum(1) % 9 == 0)]
nines = nine_sets()
print(f'catalogue: {len(lines)} lines (9-sets with 0), {len(nines)} nine-point sets; {time.time() - t0:.0f}s', flush=True)
def charvec(S, zero):  # psi_d(T) for all 80 d
    return C[:, list(S)].sum(1) + int(zero)
T0 = [b for b in B81 if logB[b] % 8 in J0]
chiT0 = np.array([sum(2 if TRB[(d, x)] == 0 else -1 for x in T0 if logB[x] < 40) for d in D80]) // 1   # T0 = +-, count each pair once
# T0 contains both b and -b; the loop above counts only representatives with log < 40, giving the full psi_d(T0)
reps = [B81[i] for i in range(8)]
Avals = {u: [A16[(3 * t + u) % 15] for t in range(5)] for u in range(3)}
sign = lambda c, a: (-1) ** TRA[(c, a)]
CS = [Z4] + A16
def contrib(cv, u):
    """K_u[d_rep][c] = sum_t (-1)^Tr(c alpha^(3t+u)) psi_(d beta^(-16t))(T_u)"""
    K = np.zeros((8, 16), dtype=np.int64)
    for di, d in enumerate(reps):
        vals = [cv[didx[mulB(-16 * t, d)]] for t in range(5)]
        for ci, c in enumerate(CS): K[di, ci] = sum(sign(c, Avals[u][t]) * vals[t] for t in range(5))
    return K
KT0 = np.array([[chiT0[didx[d]]] * 16 for d in reps])
Sz = [(4 * (chiT0[didx[d]] - s)) // 9 for d in reps]; assert all((4 * (chiT0[didx[d]] - s)) % 9 == 0 for d in reps)
exact = [di for di in range(8) if Sz[di] == 0]
print('|S_psi| per class rep:', Sz, ' exact rows:', exact, flush=True)
seeds = [(u, lines if u in I0 else nines, u in I0) for u in range(3)]
KK = [[contrib(charvec(S, z), u) for S in cat] for u, cat, z in seeds]
# hash the last seed on the exact rows
last = collections.defaultdict(list)
for i, K in enumerate(KK[2]): last[K[exact].tobytes()].append(i)
found = []; tested = 0
for i0, K0 in enumerate(KK[0]):
    for i1, K1 in enumerate(KK[1]):
        tested += 1
        need = (s - KT0 - K0 - K1)[exact]
        for i2 in last.get(need.astype(np.int64).tobytes(), []):
            tot = KT0 + K0 + K1 + KK[2][i2]
            if np.all((tot == r) | (tot == s)): found.append((i0, i1, i2))
print(f'm={m}: pairs tested {tested}, candidates passing all class-representative rows: {len(found)}; {time.time() - t0:.0f}s', flush=True)
# full check of the first few candidates
def build(i0, i1, i2):
    D = {(Z4, b) for b in T0} | {(A16[i], Z4) for i in range(15) if i % 3 in I0}
    for u, idx in zip(range(3), (i0, i1, i2)):
        S = seeds[u][1][idx]; Tu = {P[p] for p in S} | {mulB(40, P[p]) for p in S} | ({Z4} if seeds[u][2] else set())
        for t in range(5): D |= {(A16[3 * t + u], mulB(-16 * t, b)) for b in Tu}
    return D
from duality import char_values
ok = 0
for f in found[:3]:
    D = build(*f); vals = collections.Counter(v for x, v in char_values(D).items() if x != (Z4, Z4))
    good = len(D) == 35 * m and set(vals) <= {r, s}; ok += good
    print('  candidate', f, '|D| =', len(D), 'character values', dict(vals), 'PDS:', good, flush=True)
    if good and m in (5, 6, 7):
        write_connset(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), f'srg1296_{35*m}_exhaustive.txt'),
                      (1296, 35 * m, 34 + (m - 1) * (m - 2), m * (m - 1)), D, f'found by exhaustive_low.py {" ".join(sys.argv[1:])}')
print('RESULT m=%d: %s' % (m, 'solutions exist' if ok else ('NO solution (exhaustive)' if not found else 'candidates failed full check')))
