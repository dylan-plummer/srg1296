#!/usr/bin/env python3
"""Check Constructions 4.1-4.6 of the paper exactly as printed (plain Python 3, no packages).
Builds the eight Cayley graphs on F_16 x F_81 (Constructions 4.5, 4.6: on F_2^4 x F_3^4) from the printed data and verifies, for every vertex pair,
the strongly regular parameters. Runtime: a few minutes."""

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

# Construction 4.1: srg(1296, 560, 244, 240)
S = {(0, 0), (0, 2), (0, 6), (0, 11), (0, 13), (0, 14), (0, 17), (0, 18), (0, 19),
     (1, 3), (1, 6), (1, 7), (1, 8), (1, 10), (1, 11), (1, 13), (1, 15), (1, 19),
     (2, 1), (2, 2), (2, 5), (2, 9), (2, 10), (2, 11), (2, 13), (2, 16), (2, 17)}
D1 = {(Z, B81[j]) for j in range(0, 80, 4)}
D1 |= {(A16[i], B81[j]) for i in range(15) for j in range(80) if (i % 3, (j - 4 * (i // 3)) % 20) in S}

# Construction 4.2: srg(1296, 592, 268, 272); J[a] for a = 0 (key None) and a = alpha^i (key i)
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

# Construction 4.3: srg(1296, 595, 274, 272); J3[a] for a = 0 (key None) and a = alpha^i (key i), plus (a, 0) for a in F_4^*
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

# Construction 4.4 (Table 3): D(J0, I0, S) =
#   {(0, beta^j): j mod 8 in J0} u {(alpha^i, 0): i mod 3 in I0} u {(alpha^(3t+u), beta^j): (j + 16t) mod 40 in S_u}
T44 = {525: ([0, 1, 6], [0], [[0, 2, 3, 8, 10, 12, 14, 15, 20, 30, 32, 37, 39],
                              [0, 2, 3, 5, 6, 11, 12, 13, 14, 19, 20, 21, 26, 28, 30, 31, 33, 35],
                              [3, 5, 6, 11, 12, 13, 16, 17, 20, 21, 24, 25, 29, 30, 31, 32, 35, 36]]),
       490: ([0, 1, 4, 5], [0, 1], [[0, 1, 2, 8, 14, 18, 19, 23, 25, 33, 34, 37, 38],
                                    [4, 5, 6, 9, 10, 13, 21, 23, 25, 28, 29, 30, 35],
                                    [0, 2, 3, 4, 5, 6, 9, 10, 14, 15, 18, 19, 22, 25, 29, 30, 31, 34]]),
       455: ([0, 1, 2, 3, 5], [0, 1, 2], [[0, 1, 2, 4, 5, 11, 21, 28, 29, 31, 32, 37, 38],
                                          [3, 5, 6, 7, 13, 14, 19, 22, 25, 31, 33, 35, 39],
                                          [6, 10, 13, 14, 19, 22, 26, 27, 28, 30, 32, 35, 36]])}
def con44(J0, I0, S):
    D = {(Z, B81[j]) for j in range(80) if j % 8 in J0} | {(A16[i], Z) for i in range(15) if i % 3 in I0}
    D |= {(A16[3 * t + u], B81[j]) for t in range(5) for u in range(3) for j in range(80) if (j + 16 * t) % 40 in S[u]}
    return D

# Construction 4.5 (Table 4), in coordinates F_2^4 x F_3^4 (vectors written as digit strings).  sigma = (M2, M3):
#   M2(b0,b1,b2,b3) = (b1, b0+b1, b3, b2+b3),  M3(v0,v1,v2,v3) = (v0+v2, v1+v3, v2, v3);  D = {0} x +-T0  u  sigma-orbits of seeds
T45_0 = "0001 0010 0011 0012 0100 0101 0102 0112 0121 1000 1010 1012 1020 1021 1100 1111 1112 1121 1122 1212 1221"
T45 = [("0100", False, "0001 0010 0011 0101 0102 0110 0111 0120 0122"),
       ("0001", True, "0001 0010 0011 0121 1021 1101 1102 1110 1111 1112 1120 1122 1200"),
       ("0101", True, "0001 0010 0011 0112 1001 1002 1010 1011 1012 1020 1022 1121 1200"),
       ("1001", True, "0001 0012 0120 0122 1001 1010 1101 1122 1200 1210 1212 1221 1222"),
       ("1101", True, "0001 0100 0110 0111 1000 1002 1020 1100 1102 1111 1200 1211 1220")]
def con45(T45_0=T45_0, T45=T45):
    vec = lambda w: tuple(int(ch) for ch in w)
    pm = lambda P: {vec(w) for w in P.split()} | {tuple((-x) % 3 for x in vec(w)) for w in P.split()}
    M2 = lambda b: (b[1], (b[0] + b[1]) % 2, b[3], (b[2] + b[3]) % 2)
    M3 = lambda v: ((v[0] + v[2]) % 3, (v[1] + v[3]) % 3, v[2], v[3])
    D = {((0, 0, 0, 0), t) for t in pm(T45_0)}
    for lab, z, P in T45:
        b, T = vec(lab), pm(P) | ({(0, 0, 0, 0)} if z else set())
        for _ in range(3):
            D |= {(b, t) for t in T}; b = M2(b); T = {M3(t) for t in T}
    # identify F_2^4 x F_3^4 with V below: b -> coefficient tuple (b0..b3) of F_16, v -> coefficient tuple (v0..v3) of F_81
    return D

# Construction 4.6 (Table 5), same coordinates and the same sigma as Construction 4.5
T46_0 = "0001 0010 0012 0100 0101 0102 0112 0121 1000 1010 1012 1020 1021 1112 1121 1212 1221"
T46 = [("0110", False, "0001 0010 0011 0101 0102 0110 0111 0120 0122"),
       ("0111", False, "0102 0110 0121 1001 1012 1020 1200 1211 1222"),
       ("0001", True, "0001 0010 0011 0012 1001 1002 1010 1011 1020 1022 1200 1212 1221"),
       ("0100", True, "0010 0100 0110 0120 1002 1020 1022 1111 1112 1121 1200 1201 1210"),
       ("0101", True, "0001 0101 0110 0112 1000 1011 1021 1101 1102 1122 1200 1202 1210")]

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
ok4 = [check(con44(*T44[k]), k, lam, mu) for k, lam, mu in ((525, 216, 210), (490, 190, 182), (455, 166, 156))]
ok5 = check(con45(), 420, 144, 132)
ok6 = check(con45(T46_0, T46), 385, 124, 110)
raise SystemExit(0 if ok1 and ok2 and ok3 and all(ok4) and ok5 and ok6 else 1)
