#!/usr/bin/env python3
r"""Do our srg(1296,595,274,272) Cayley graphs have Kovaľ's 'two deleted spreads' structure?
For every subgroup U = U2 x U3 of order 36 of F16 x F81 (U2 of order 4 in F_2^4, U3 of order 9 in F_3^4):
  - is U disjoint from D?            (then the cosets of U are a spread of 36-cocliques)
  - is D + (U \ 0) a PDS?             (LS_18(36) = (1296,630,306,306); 'adding one spread')
and for orthogonal pairs U, U' (U & U' = 0) of such subgroups:
  - is D + (U\0) + (U'\0) a PDS?      (LS_19(36) = (1296,665,340,342); 'adding two spreads')
Also run on the (1296,560,244,240) and (1296,592,268,272) graphs for comparison (with their own targets)."""
import os, glob, itertools
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
CONN = os.path.join(HERE, "..", "..", "connsets")

def powers(p, poly, n):
    cur = [1] + [0] * (n - 1); out = []
    for _ in range(p ** n - 1):
        out.append(tuple(cur)); top = cur[-1]; cur = [0] + cur[:-1]
        cur = [(c - top * poly[i]) % p for i, c in enumerate(cur)]
    return out
A16 = powers(2, [1, 0, 0, 1, 1], 4); B81 = powers(3, [2, 0, 0, 1, 1], 4)
SHAPE = (2, 2, 2, 2, 3, 3, 3, 3)

def load(fn):
    L = [l.split() for l in open(fn) if l.strip() and not l.startswith("#")]
    prm = tuple(map(int, L[0])); M = np.zeros(SHAPE, bool)
    for i, j in L[1:]:
        a = (0,) * 4 if i == "-" else A16[int(i)]; b = (0,) * 4 if j == "-" else B81[int(j)]
        M[a + b] = True
    return prm, M

def subspaces(p, n, k):
    """all k-dim subspaces of F_p^n as boolean masks over F_p^n (shape (p,)*n)."""
    vecs = list(itertools.product(range(p), repeat=n)); seen = set(); out = []
    for basis in itertools.combinations(vecs[1:], k):
        span = {tuple(sum(c * b[t] for c, b in zip(co, basis)) % p for t in range(n)) for co in itertools.product(range(p), repeat=k)}
        if len(span) != p ** k: continue
        key = frozenset(span)
        if key in seen: continue
        seen.add(key); m = np.zeros((p,) * n, bool)
        for s in span: m[s] = True
        out.append(m)
    return out

def is_pds(M, r, s):
    F = np.fft.fftn(M.astype(float)).ravel()[1:]
    return bool(np.all((np.abs(F - r) < 1e-6) | (np.abs(F - s) < 1e-6)))

U2s, U3s = subspaces(2, 4, 2), subspaces(3, 4, 2)
print(f"{len(U2s)} subgroups of order 4 in Z_2^4, {len(U3s)} of order 9 in Z_3^4")
Us = [(np.einsum("abcd,efgh->abcdefgh", u2, u3), a, b) for a, u2 in enumerate(U2s) for b, u3 in enumerate(U3s)]
if __name__ == "__main__":
    targets = {595: [((306, 306), "LS18"), ((340, 342), "LS19")], 560: [], 592: []}
    EV = {(630, 306, 306): (18, -18), (665, 340, 342): (17, -19)}
    for fn in sorted(glob.glob(os.path.join(CONN, "srg1296_*_class*.txt")), key=lambda f: (f.split("_class")[0], int(f.split("_class")[1][:-4]))):
        prm, D = load(fn)
        free = [(U, a, b) for U, a, b in Us if not (U & D).any()]
        line = f"{os.path.basename(fn)[8:-4]:22s} subgroups of order 36 missing D: {len(free):4d}"
        if prm[1] == 595:
            zero = np.zeros(SHAPE, bool); zero[(0,) * 8] = True
            one = [x for x in free if is_pds(D | (x[0] & ~zero), 18, -18)]
            two = 0; pairs = 0
            for (U, a, b), (V, c, d) in itertools.combinations(one, 2):
                if (U & V).sum() == 1:                      # orthogonal: U & V = {0}
                    pairs += 1
                    if is_pds(D | ((U | V) & ~zero), 17, -19): two += 1
            line += f"  | D+U* is LS18: {len(one):3d}  | orthogonal pairs: {pairs}, of which D+U*+V* is LS19: {two}"
        print(line, flush=True)
