"""Conjugacy-class representatives of GL(n, p) (p prime) via rational canonical forms:
block-diagonal companion matrices of f^e for irreducible monic f != x."""
import itertools
import numpy as np


def polymul(a, b, p):
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            r[i + j] = (r[i + j] + x * y) % p
    return r


def irreducibles(d, p):
    """monic irreducible polys of degree d over F_p, coefficient lists low->high, excluding x."""
    out = []
    # brute force: monic polys of degree d with nonzero constant term, test no factor of degree <= d/2
    smaller = {}
    def all_monic(deg):
        for c in itertools.product(range(p), repeat=deg):
            yield list(c) + [1]
    def divides(f, g):  # does f divide g
        g = g[:]
        while len(g) >= len(f):
            c = g[-1] % p
            if c:
                sh = len(g) - len(f)
                for i, x in enumerate(f):
                    g[sh + i] = (g[sh + i] - c * x) % p
            g.pop()
        return all(x % p == 0 for x in g)
    irr_small = []
    for e in range(1, d // 2 + 1):
        irr_small += [f for f in irreducibles(e, p)] if e < d else []
    for f in all_monic(d):
        if f[0] == 0:
            continue
        if d == 1 or not any(divides(g, f) for g in irr_small + ([[0, 1]] if d > 1 else [])):
            out.append(f)
    return out


def companion(f, p):
    d = len(f) - 1
    C = np.zeros((d, d), dtype=np.int64)
    for i in range(1, d):
        C[i, i - 1] = 1
    for i in range(d):
        C[i, d - 1] = (-f[i]) % p
    return C


def partitions(n, mx=None):
    if mx is None: mx = n
    if n == 0:
        yield []; return
    for i in range(min(n, mx), 0, -1):
        for q in partitions(n - i, i):
            yield [i] + q


def class_reps(n, p):
    irr = {d: irreducibles(d, p) for d in range(1, n + 1)}
    atoms = []  # (degree, poly)
    for d in range(1, n + 1):
        for f in irr[d]:
            atoms.append((d, tuple(f)))
    # choose for each atom a partition; total degree n
    reps = []
    def rec(i, remaining, chosen):
        if remaining == 0:
            reps.append(list(chosen)); return
        if i == len(atoms):
            return
        d, f = atoms[i]
        rec(i + 1, remaining, chosen)
        for tot in range(1, remaining // d + 1):
            for lam in partitions(tot):
                chosen.append((f, lam))
                rec(i + 1, remaining - d * tot, chosen)
                chosen.pop()
    rec(0, n, [])
    mats = []
    for rep in reps:
        blocks = []
        for f, lam in rep:
            for part in lam:
                g = [1]
                for _ in range(part):
                    g = polymul(g, list(f), p)
                blocks.append(companion(g, p))
        M = np.zeros((n, n), dtype=np.int64)
        o = 0
        for B in blocks:
            s = B.shape[0]
            M[o:o + s, o:o + s] = B
            o += s
        mats.append(M)
    return mats


def order_mod(M, p, cap=10 ** 7):
    n = M.shape[0]
    I = np.eye(n, dtype=np.int64)
    A = M % p
    k = 1
    while not (A == I).all():
        A = (A @ M) % p
        k += 1
        if k > cap:
            return None
    return k


if __name__ == "__main__":
    import sys
    n, p = int(sys.argv[1]), int(sys.argv[2])
    R = class_reps(n, p)
    print(len(R))
