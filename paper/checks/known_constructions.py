#!/usr/bin/env python3
"""Checks against known constructions (plain Python 3 + numpy).
1. Multiplicative stabiliser in F16* x F81* (and with Galois maps) of every connection set in ../../connsets.
2. Fernandez-Alcober-Kwashira-Martinez standard cyclotomy of order e | gcd(15,80) = 5 on F16 x F81:
   every union of orbits of size 560 or 592 is tested for being a PDS with our parameters."""
import os, glob, itertools, cmath
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

def load(fn):
    L = [l.split() for l in open(fn) if l.strip() and not l.startswith("#")]
    prm = tuple(map(int, L[0]))
    D = {(-1 if i == "-" else int(i), -1 if j == "-" else int(j)) for i, j in L[1:]}   # (log a, log b), -1 = zero
    return prm, D

def mult(D, x, y):            # multiply by (alpha^x, beta^y)
    return {((i + x) % 15 if i >= 0 else -1, (j + y) % 80 if j >= 0 else -1) for i, j in D}
def frob(D, u, w):            # (a,b) -> (a^(2^u), b^(3^w))
    return {((i * 2 ** u) % 15 if i >= 0 else -1, (j * 3 ** w) % 80 if j >= 0 else -1) for i, j in D}

print("1. multiplicative / semilinear stabilisers")
for fn in sorted(glob.glob(os.path.join(CONN, "srg1296_5[69]*_class*.txt")), key=lambda f: (f.split("_class")[0], int(f.split("_class")[1][:-4]))):
    prm, D = load(fn)
    stab = [(x, y) for x in range(15) for y in range(80) if mult(D, x, y) == D]
    full = [(x, y, u, w) for x in range(15) for y in range(80) for u in range(4) for w in range(4)
            if mult(frob(D, u, w), x, y) == D]
    print(f"  {os.path.basename(fn)[8:-4]:22s} |Stab_mult| = {len(stab):3d}   |Stab_semilinear| = {len(full):3d}")

# --- FKM standard cyclotomy of order 5: H = {(t1^r1, t2^r2): r1 + r2 = 0 mod 5} for primitive roots t1, t2.
# With t1 = alpha^a, t2 = beta^b (a, b units), H = {(alpha^x, beta^y): a' x + b' y = 0 mod 5} with a', b' units mod 5.
print("2. FKM standard cyclotomy of order e = 5 (e = 1 gives orbit sizes 1, 15, 80, 1200 only)")
TRA = [None] * 15; TRB = [None] * 80
for i in range(15):
    s = [0] * 4
    for t in range(4): s = [(x + y) % 2 for x, y in zip(s, A16[(i * 2 ** t) % 15])]
    TRA[i] = s[0]
for j in range(80):
    s = [0] * 4
    for t in range(4): s = [(x + y) % 3 for x, y in zip(s, B81[(j * 3 ** t) % 80])]
    TRB[j] = s[0]
w3 = cmath.exp(2j * cmath.pi / 3)
def char_values(D):
    vals = set()
    for c in range(-1, 15):
        for e in range(-1, 80):
            if c < 0 and e < 0: continue
            z = sum((-1) ** (0 if (c < 0 or i < 0) else TRA[(c + i) % 15]) * w3 ** (0 if (e < 0 or j < 0) else TRB[(e + j) % 80]) for i, j in D)
            vals.add(round(z.real, 6))
            if len(vals) > 2: return vals
    return vals
targets = {560: (20, -16), 592: (16, -20)}
tested = 0; found = 0
for ap in (1,):                      # scaling (a', b') by a unit gives the same H; fix a' = 1
    for bp in range(1, 5):
        orbits = {r: {(x, y) for x in range(15) for y in range(80) if (ap * x + bp * y) % 5 == r} for r in range(5)}
        F1 = {(x, -1) for x in range(15)}; F2 = {(-1, y) for y in range(80)}
        pieces = [orbits[r] for r in range(5)] + [F1, F2]
        for mask in range(1, 2 ** 7):
            D = set().union(*[pieces[t] for t in range(7) if mask >> t & 1])
            if len(D) in targets:
                tested += 1
                vals = char_values(D)
                if vals == set(map(float, targets[len(D)])): found += 1; print("   PDS!", bp, mask)
print(f"   unions of orbits of size 560 or 592 tested: {tested}; PDS with our parameters among them: {found}")
