#!/usr/bin/env python3
"""Structure of all connection sets in ../connsets (Section 6 of the paper).  Needs numpy and python-igraph.
For each Cay(G, D), G = F_16 x F_81, with restricted eigenvalues r > s, Delta = r - s, it reports:
  * ring graph   Cay(F_81, T_0),  T_0 = {b : (0,b) in D}, and whether T_0 is a union of cyclotomic classes of index 8;
  * column graph Cay(F_16, U_0),  U_0 = {a : (a,0) in D}, and whether U_0 is a union of cosets of mu_5 = <alpha^3>;
  * the ring sizes |T_a| (a != 0);
  * Lemma 6.1 (divisibility): for a != 0 every nontrivial character sum of T_a is divisible by Delta/gcd(Delta,16),
    and for b != 0 every Walsh sum of U_b = {a : (a,b) in D} is divisible by Delta/gcd(Delta,81);
  * the weights |S_beta| of the positive sets S_beta = {alpha : chi_(alpha,beta)(D) = r} (beta != 0);
  * |Aut Cay(G,D)| (bliss), and whether Cay(G, D*) is isomorphic to Cay(G, D), D* the dual set (duality.py);
  * the isomorphism classes among all connection sets with the same parameters (should be all distinct).
usage: python3 family_structure.py > family_structure.out"""
import os, glob, math, collections, sys
import numpy as np, igraph as ig
from fields import *
from duality import char_values, dual_set, TRA, TRB, ELA, ELB

HERE = os.path.dirname(os.path.abspath(__file__)); CONN = os.path.join(HERE, '..', 'connsets')
F2 = [tuple((t >> s) & 1 for s in range(4)) for t in range(16)]; F3 = [tuple((t // 3 ** s) % 3 for s in range(4)) for t in range(81)]
i2 = {a: i for i, a in enumerate(F2)}; i3 = {b: i for i, b in enumerate(F3)}
d2 = np.array([[i2[tuple((p - q) % 2 for p, q in zip(x, y))] for y in F2] for x in F2])
d3 = np.array([[i3[tuple((p - q) % 3 for p, q in zip(x, y))] for y in F3] for x in F3])


def srg_of(n, S, sub):
    S = set(S)
    if not S: return 'empty'
    if len(S) == n - 1: return f'K{n}'
    c = collections.Counter(sub(x, y) for x in S for y in S if x != y)
    lam = {c[s] for s in S}; mu = {c[t] for t in range(n) if t != 0 and t not in S}
    return f'srg({n},{len(S)},{lam.pop()},{mu.pop()})' if len(lam) == 1 and len(mu) == 1 else f'{len(S)}-set, not srg'


def cayley(D):
    M = np.zeros((16, 81), bool)
    for a, b in D: M[i2[a], i3[b]] = True
    ea, eb = np.divmod(np.arange(1296), 81)
    return ig.Graph.Adjacency(M[d2[ea][:, ea], d3[eb][:, eb]].tolist(), mode='undirected')


files = sorted(glob.glob(os.path.join(CONN, 'srg1296_*_class*.txt')), key=lambda f: (-int(f.split('_')[1]), int(f.split('class')[1][:-4])))
roman = {560: 'I', 592: 'II', 595: 'III', 525: 'IV', 490: 'V', 455: 'VI', 420: 'VII', 385: 'VIII'}
graphs = collections.defaultdict(list)
for fn in files:
    (v, k, lam, mu), D = read_connset(fn); name = f'{roman[k]}.{int(fn.split("class")[1][:-4])}'
    g = math.isqrt((lam - mu) ** 2 + 4 * (k - mu)); r, s = (lam - mu + g) // 2, (lam - mu - g) // 2; Delta = r - s
    typ = f'LS_{-s}(36)' if k == -s * 35 else f'NL_{r}(36)' if k == 37 * r else '?'
    T = collections.defaultdict(set); U = collections.defaultdict(set)
    for a, b in D: T[a].add(b); U[b].add(a)
    T0, U0 = T[Z4], U[Z4]
    ring = srg_of(81, [i3[b] for b in T0], lambda x, y: d3[x, y])
    col = srg_of(16, [i2[a] for a in U0], lambda x, y: d2[x, y])
    cyc8 = all(mulB(8, b) in T0 for b in T0); mu5 = all(mulA(3, a) in U0 for a in U0)
    g8 = len({logB[b] % 8 for b in T0}) if T0 and cyc8 else None
    sizes = collections.Counter(len(T[a]) for a in A16)
    # divisibility (Lemma 6.1): characters of each T_a, a != 0, and Walsh sums of each U_b, b != 0
    gA, gB = Delta // math.gcd(Delta, 16), Delta // math.gcd(Delta, 81)
    w3 = lambda d, S: sum({0: 2, 1: -1, 2: -1}[TRB[(d, b)]] for b in S)       # 2 * real character sum
    okA = all(w3(d, T[a]) % (2 * gA) == 0 for a in A16 for d in B81)
    okB = all(sum((-1) ** TRA[(c, a)] for a in U[b]) % gB == 0 for b in B81 for c in A16)
    chi = char_values(D)
    weights = collections.Counter(sum(1 for c in ELA if chi[(c, d)] == r) for d in B81)
    G = cayley(D); aut = G.count_automorphisms(sh='fl')
    Ds = dual_set(D, k); selfdual = G.isomorphic(cayley(Ds))
    graphs[k].append((name, G))
    print(f'{name:6s} {typ:10s} |Aut|={aut:7d}  ring {ring}' + (f' = {g8} classes of index 8' if g8 else '') +
          f'; column {col}' + (f' = {len(U0) // 5} cosets of mu5' if U0 and mu5 else '') +
          f'; ring sizes {dict(sorted(sizes.items()))}; divisibility {gA}|chi(T_a): {okA}, {gB}|W(U_b): {okB};'
          f' |S_beta| {dict(sorted(weights.items()))}; dual isomorphic: {selfdual}', flush=True)
print('pairwise non-isomorphic within each parameter set:')
for k, L in graphs.items():
    ok = all(not L[i][1].isomorphic(L[j][1]) for i in range(len(L)) for j in range(i + 1, len(L)))
    print(f'  k={k}: {len(L)} graphs, pairwise non-isomorphic: {ok}', flush=True)
