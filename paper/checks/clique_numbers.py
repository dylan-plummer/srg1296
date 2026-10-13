#!/usr/bin/env python3
"""Exact clique numbers of the graphs in ../../connsets (Cayley graphs on F16 x F81).

Gamma is vertex-transitive, so a maximum clique may be assumed to contain the vertex 0 = (0,0) and an
edge {0, x}; hence omega(Gamma) = 2 + max_x omega(N(0) & N(x)). The neighbours x of 0 are reduced to orbit
representatives under the known automorphisms fixing 0: the maps (a,b) -> (g1 s1(a), g2 s2(b)) with
(g1, g2) in F16* x F81* and s1, s2 Frobenius powers, that preserve D. Each subproblem is solved exactly
with igraph's clique_number (Cliquer). A witness clique of size omega is produced and checked.

usage: python3 clique_numbers.py [NPROC]   (needs numpy, igraph)"""
import os, sys, glob, time
from multiprocessing import Pool
import numpy as np
import igraph as ig

HERE = os.path.dirname(os.path.abspath(__file__))
CONN = os.path.join(HERE, "..", "..", "connsets")

def powers(p, poly, n):
    cur = [1] + [0] * (n - 1); out = []
    for _ in range(p ** n - 1):
        out.append(tuple(cur)); top = cur[-1]; cur = [0] + cur[:-1]
        cur = [(c - top * poly[i]) % p for i, c in enumerate(cur)]
    return out
A16 = powers(2, [1, 0, 0, 1, 1], 4); B81 = powers(3, [2, 0, 0, 1, 1], 4)
encA = lambda t: sum(c << s for s, c in enumerate(t)); encB = lambda t: sum(c * 3 ** s for s, c in enumerate(t))
# vertices indexed by log coordinates: vid(i, j) with i in {-1..14}, j in {-1..79} (-1 = zero)
def vid(i, j): return (i + 1) * 81 + (j + 1)
ENC = np.array([[(0 if i < 0 else encA(A16[i])), (0 if j < 0 else encB(B81[j]))] for i in range(-1, 15) for j in range(-1, 80)])
SUBB = np.array([[sum(((x // 3 ** s) % 3 - (y // 3 ** s) % 3) % 3 * 3 ** s for s in range(4)) for y in range(81)] for x in range(81)])

def load(fn):
    L = [l.split() for l in open(fn) if l.strip() and not l.startswith("#")]
    prm = tuple(map(int, L[0]))
    D = {(-1 if i == "-" else int(i), -1 if j == "-" else int(j)) for i, j in L[1:]}
    inD = np.zeros((16, 81), bool)
    for i, j in D: inD[0 if i < 0 else encA(A16[i]), 0 if j < 0 else encB(B81[j])] = True
    a, b = ENC[:, 0], ENC[:, 1]
    A = inD[a[:, None] ^ a[None, :], SUBB[b[:, None], b[None, :]]]
    np.fill_diagonal(A, False)
    return prm, D, A

def autos_fixing_zero(D):
    """all maps (i,j) -> (i*2^u + x, j*3^w + y) (logs) preserving D; returned as vertex permutations."""
    perms = []
    for u in range(4):
        for w in range(4):
            Df = {((i * 2 ** u) % 15 if i >= 0 else -1, (j * 3 ** w) % 80 if j >= 0 else -1) for i, j in D}
            for x in range(15):
                for y in range(80):
                    if {((i + x) % 15 if i >= 0 else -1, (j + y) % 80 if j >= 0 else -1) for i, j in Df} == D:
                        perm = np.empty(1296, int)
                        for i in range(-1, 15):
                            for j in range(-1, 80):
                                ii = (i * 2 ** u + x) % 15 if i >= 0 else -1
                                jj = (j * 3 ** w + y) % 80 if j >= 0 else -1
                                perm[vid(i, j)] = vid(ii, jj)
                        perms.append(perm)
    return perms

def work(fn):
    t0 = time.time()
    prm, D, A = load(fn)
    k, lam, mu = prm[1:]
    assert (A.sum(1) == k).all()
    perms = autos_fixing_zero(D)
    for p in perms:                                   # each must be an automorphism fixing 0
        assert p[0] == 0 and (A[np.ix_(p, p)] == A).all()
    z = vid(-1, -1); assert z == 0
    nb = np.where(A[0])[0]
    seen, reps = set(), []
    for x in nb:
        if x in seen: continue
        orb = {int(p[x]) for p in perms}; seen |= orb; reps.append(int(x))
    best, bestx = -1, None
    for x in reps:
        H = np.where(A[0] & A[x])[0]
        g = ig.Graph.Adjacency(A[np.ix_(H, H)].astype(int).tolist(), mode="undirected")
        w = g.clique_number()
        if w > best: best, bestx = w, x
    omega = best + 2
    # witness
    H = np.where(A[0] & A[bestx])[0]
    g = ig.Graph.Adjacency(A[np.ix_(H, H)].astype(int).tolist(), mode="undirected")
    cl = [0, bestx] + [int(H[v]) for v in g.largest_cliques()[0]] if best > 0 else [0, bestx]
    assert len(cl) == omega and all(A[a, b] for a in cl for b in cl if a != b)
    return os.path.basename(fn)[8:-4], omega, len(perms), len(reps), time.time() - t0

if __name__ == "__main__":
    nproc = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    files = sorted(glob.glob(os.path.join(CONN, "srg1296_560*_class*.txt")) + glob.glob(os.path.join(CONN, "srg1296_592*_class*.txt")) + glob.glob(os.path.join(CONN, "srg1296_595*_class*.txt")),
                   key=lambda f: (f.split("_class")[0], int(f.split("_class")[1][:-4])))
    with Pool(nproc) as pool:
        for name, omega, naut, nrep, dt in pool.imap(work, files):
            print(f"{name:22s} omega = {omega:3d}   (|known Stab(0)| = {naut}, {nrep} edge orbits, {dt:.0f}s)", flush=True)
