#!/usr/bin/env python3
"""Check Constructions 1, 2 and 3 of the paper exactly as printed (plain Python 3, no packages).
Builds the three Cayley graphs on F_16 x F_81 from the printed data and verifies, for every vertex pair,
the strongly regular parameters. Runtime: about a minute."""

# F_16 = F_2[x]/(x^4+x^3+1), alpha = x ;  F_81 = F_3[y]/(y^4+y^3+2), beta = y
def powers(p, poly, n):
    cur = [1] + [0] * (n - 1); out = []
    for _ in range(p ** n - 1):
        out.append(tuple(cur)); top = cur[-1]; cur = [0] + cur[:-1]
        cur = [(c - top * poly[i]) % p for i, c in enumerate(cur)]
    assert len(set(out)) == p ** n - 1
    return out

A16 = powers(2, [1, 0, 0, 1, 1], 4)
B81 = powers(3, [2, 0, 0, 1, 1], 4)
Z = (0, 0, 0, 0)

# Construction 1: srg(1296, 560, 244, 240)
S = {(0, 0), (0, 2), (0, 6), (0, 11), (0, 13), (0, 14), (0, 17), (0, 18), (0, 19),
     (1, 3), (1, 6), (1, 7), (1, 8), (1, 10), (1, 11), (1, 13), (1, 15), (1, 19),
     (2, 1), (2, 2), (2, 5), (2, 9), (2, 10), (2, 11), (2, 13), (2, 16), (2, 17)}
D1 = {(Z, B81[j]) for j in range(0, 80, 4)}
D1 |= {(A16[i], B81[j]) for i in range(15) for j in range(80) if (i % 3, (j - 4 * (i // 3)) % 20) in S}

# Construction 2: srg(1296, 592, 268, 272); J[a] for a = 0 (key None) and a = alpha^i (key i)
J = {None: [0, 1, 3, 4, 7, 8, 9, 11, 12, 13, 16, 17, 19],
     0: [0, 2, 6, 11, 13, 14, 17, 18, 19], 1: [4, 6, 7, 10, 11, 12, 13, 15, 19],
     2: [1, 5, 10, 12, 13, 16, 17, 18, 19], 3: [1, 2, 7, 9, 10, 14, 15, 18, 19],
     4: [3, 8, 10, 11, 14, 15, 16, 17, 19], 5: [0, 1, 5, 6, 9, 13, 14, 15, 17],
     6: [1, 2, 3, 5, 6, 7, 10, 14, 17], 7: [0, 1, 5, 6, 7, 8, 10, 13, 14],
     8: [2, 4, 5, 8, 9, 10, 11, 13, 17], 9: [3, 5, 6, 7, 9, 10, 13, 14, 18],
     10: [0, 2, 3, 5, 7, 11, 15, 18, 19], 11: [0, 2, 7, 9, 10, 11, 15, 16, 18],
     12: [1, 2, 3, 6, 9, 10, 11, 15, 18], 13: [0, 3, 5, 6, 9, 10, 12, 14, 17],
     14: [0, 1, 2, 3, 4, 10, 15, 18, 19]}
D2 = {(Z if a is None else A16[a], B81[j]) for a, js in J.items() for j in range(80) if j % 20 in js}

# Construction 3: srg(1296, 595, 274, 272); J3[a] for a = 0 (key None) and a = alpha^i (key i), plus (a, 0) for a in F_4^*
J3 = {None: [0, 4, 5, 8, 12, 15, 16],
      0: [0, 2, 5, 6, 11, 13, 14, 15, 17, 18, 19], 1: [5, 9, 10, 12, 13, 14, 16, 18, 19],
      2: [2, 7, 8, 10, 14, 15, 16, 17, 19], 3: [3, 4, 5, 6, 10, 12, 13, 14, 19],
      4: [1, 2, 4, 5, 6, 8, 10, 11, 17], 5: [0, 2, 3, 4, 6, 7, 11, 14, 16, 18, 19],
      6: [2, 9, 10, 12, 15, 16, 17, 18, 19], 7: [0, 2, 4, 5, 9, 10, 11, 16, 17],
      8: [3, 4, 6, 10, 11, 12, 13, 15, 18], 9: [1, 2, 4, 8, 10, 11, 13, 15, 18],
      10: [0, 1, 2, 6, 8, 9, 12, 13, 14, 17, 18], 11: [0, 3, 8, 10, 12, 14, 15, 17, 19],
      12: [5, 6, 7, 8, 10, 11, 14, 16, 17], 13: [0, 1, 4, 5, 10, 13, 16, 18, 19],
      14: [0, 6, 7, 8, 10, 11, 12, 13, 15]}
D3 = {(Z if a is None else A16[a], B81[j]) for a, js in J3.items() for j in range(80) if j % 20 in js}
D3 |= {(A16[i], Z) for i in (0, 5, 10)}          # F_4^* = {1, alpha^5, alpha^10}

F2 = [tuple((t >> s) & 1 for s in range(4)) for t in range(16)]
F3 = [tuple((t // 3 ** s) % 3 for s in range(4)) for t in range(81)]
V = [(a, b) for a in F2 for b in F3]
sub = lambda x, y: (tuple((p - q) % 2 for p, q in zip(x[0], y[0])), tuple((p - q) % 3 for p, q in zip(x[1], y[1])))

def check(D, k, lam, mu):
    assert len(D) == k and (Z, Z) not in D
    assert all((d[0], tuple((-c) % 3 for c in d[1])) in D for d in D)      # D = -D
    n = len(V); nb = []
    for u in range(n):
        m = 0
        for w in range(n):
            if u != w and sub(V[u], V[w]) in D: m |= 1 << w
        nb.append(m)
    ok = all(bin(m).count("1") == k for m in nb) and all(
        bin(nb[u] & nb[w]).count("1") == (lam if nb[u] >> w & 1 else mu) for u in range(n) for w in range(u + 1, n))
    print(f"srg(1296,{k},{lam},{mu}):", "VERIFIED" if ok else "FAILED")
    return ok

ok1 = check(D1, 560, 244, 240)
ok2 = check(D2, 592, 268, 272)
ok3 = check(D3, 595, 274, 272)
raise SystemExit(0 if ok1 and ok2 and ok3 else 1)
