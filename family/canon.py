"""Closed form of a connection set invariant under (alpha^12, beta^16) and -1 (the order-5 model of the paper):
   D = {(0, beta^j) : j mod 8 in J0} u {(alpha^i, 0) : i mod 3 in I0} u {(alpha^(3t+u), beta^j) : (j + 16 t) mod 40 in S_u},
with t in 0..4, u in 0..2.  canonical(D) picks the lexicographically smallest image under the maps
(a, b) -> (alpha^x a^(2^u), beta^y b^(3^w)) with 2^u = 3^w (mod 5), which preserve the model."""
from fields import *


def logs(D): return {(None if a == Z4 else logA[a], None if b == Z4 else logB[b]) for a, b in D}


def closed_form(L):
    J0 = sorted({j % 8 for i, j in L if i is None})
    I0 = sorted({i % 3 for i, j in L if j is None})
    S = [sorted({j % 40 for i, j in L if i == u and j is not None}) for u in range(3)]
    return J0, I0, S


def rebuild(J0, I0, S):
    L = {(None, j) for j in range(80) if j % 8 in J0} | {(i, None) for i in range(15) if i % 3 in I0}
    L |= {(3 * t + u, j) for t in range(5) for u in range(3) for j in range(80) if (j + 16 * t) % 40 in S[u]}
    return L


def transforms():
    for u in range(4):
        for w in range(4):
            if (2 ** u - 3 ** w) % 5: continue
            for x in range(15):
                for y in range(80): yield x, y, u, w


def apply(L, x, y, u, w):
    return {(None if i is None else (x + i * 2 ** u) % 15, None if j is None else (y + j * 3 ** w) % 80) for i, j in L}


def canonical(D):
    L0 = logs(D); best = None
    for t in transforms():
        L = apply(L0, *t); cf = closed_form(L)
        if rebuild(*cf) != L: continue          # not in closed form (cannot happen for invariant D)
        key = (cf[0], cf[1], cf[2])
        if best is None or key < best[0]: best = (key, t, L)
    return best
