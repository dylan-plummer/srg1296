#!/usr/bin/env python3
"""Standalone builder + verifier for the Cayley graphs on F_16 x F_81 found in this project.
Plain Python 3, no packages; shares no code with the search.

usage: python3 build_from_logs.py CONNSET.txt [out.g6]

CONNSET.txt: first line "v k lambda mu"; then one element of D per line as "i j", meaning
(alpha^i, beta^j), with "-" standing for the field element 0.
Fields: F_16 = F_2[x]/(x^4 + x^3 + 1), alpha = x;  F_81 = F_3[y]/(y^4 + y^3 + 2), beta = y.
Vertices are pairs (a, b); (a, b) ~ (a', b') iff (a - a', b - b') is in D.
The script checks D = -D, 0 not in D, |D| = k, then degree and the lambda/mu counts for every pair.
"""
import sys


def powers(p, poly, n):
    cur = [1] + [0] * (n - 1); out = []
    for _ in range(p ** n - 1):
        out.append(tuple(cur))
        top = cur[-1]; cur = [0] + cur[:-1]
        cur = [(c - top * poly[i]) % p for i, c in enumerate(cur)]
    assert len(set(out)) == p ** n - 1, "not primitive"
    return out


A16 = powers(2, [1, 0, 0, 1, 1], 4)
B81 = powers(3, [2, 0, 0, 1, 1], 4)
Z16, Z81 = (0,) * 4, (0,) * 4

lines = [l.split() for l in open(sys.argv[1]) if l.strip() and not l.startswith("#")]
v, k, lam, mu = map(int, lines[0])
D = set()
for i, j in lines[1:]:
    D.add((Z16 if i == "-" else A16[int(i) % 15], Z81 if j == "-" else B81[int(j) % 80]))
assert v == 1296 and len(D) == k, "wrong size"
neg = lambda x: (x[0], tuple((-c) % 3 for c in x[1]))
assert all(neg(d) in D for d in D) and (Z16, Z81) not in D

F2 = [tuple((t >> s) & 1 for s in range(4)) for t in range(16)]
F3 = [tuple((t // 3 ** s) % 3 for s in range(4)) for t in range(81)]
V = [(a, b) for a in F2 for b in F3]
sub = lambda x, y: (tuple((p - q) % 2 for p, q in zip(x[0], y[0])), tuple((p - q) % 3 for p, q in zip(x[1], y[1])))
n = len(V)
nb = [0] * n
for u in range(n):
    m = 0
    for w in range(n):
        if u != w and sub(V[u], V[w]) in D: m |= 1 << w
    nb[u] = m
ok = all(bin(m).count("1") == k for m in nb) and all(
    bin(nb[u] & nb[w]).count("1") == (lam if nb[u] >> w & 1 else mu) for u in range(n) for w in range(u + 1, n))
print(f"srg({v},{k},{lam},{mu}) from {sys.argv[1]}:", "VERIFIED" if ok else "FAILED")
if ok and len(sys.argv) > 2:
    bits = [(nb[u] >> w) & 1 for w in range(1, n) for u in range(w)]
    bits += [0] * (-len(bits) % 6)
    data = [126, 63 + (n >> 12 & 63), 63 + (n >> 6 & 63), 63 + (n & 63)]
    data += [63 + int("".join(map(str, bits[t:t + 6])), 2) for t in range(0, len(bits), 6)]
    with open(sys.argv[2], "w", newline="\n") as f: f.write(bytes(data).decode() + "\n")
sys.exit(0 if ok else 1)
