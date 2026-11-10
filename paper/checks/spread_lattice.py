#!/usr/bin/env python3
"""Spread structure of the 18 connection sets, checked from scratch (own field arithmetic, character criterion).

For each graph Cay(G, D), G = F_16 x F_81 = Z_2^4 x Z_3^4 (as vector spaces F_2^4 x F_3^4):
  * all 35 * 130 subgroups U = U_2 x U_3 of order 36 (2-dim subspaces), those with U \\ {0} disjoint from D;
  * the orthogonality graph on them (U, V adjacent iff U n V = {0});
  * for every set of pairwise orthogonal such subgroups U_1..U_j: is D u U_1* u ... u U_j* an LS_{m+j}(36) PDS?
    (character criterion: all non-trivial character values in {36-m-j, -(m+j)}; D = -D, 0 not in D);
  * Bush type (Kovaľ, Thm 4.3 / Section 4.10): for E = D u U* u V* u W* of type LS_19(36), the matrix K = 2(A_E + I) - J
    indexed by (U-coset, V-coset) is a symmetric Hadamard matrix of Bush type with unit-diagonal blocks, and the
    W-cosets form a third spread orthogonal to block rows and block columns.
usage: python spread_lattice.py [OUT.json]
"""
import os, sys, json, itertools
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CONN = os.path.join(HERE, "..", "..", "connsets")


def powers(p, poly, n):
    cur = [1] + [0] * (n - 1); out = []
    for _ in range(p ** n - 1):
        out.append(tuple(cur)); top = cur[-1]; cur = [0] + cur[:-1]
        cur = [(c - top * poly[i]) % p for i, c in enumerate(cur)]
    assert len(set(out)) == p ** n - 1
    return out


A16 = powers(2, [1, 0, 0, 1, 1], 4)          # alpha^i in F_2[x]/(x^4+x^3+1), coefficient vectors
B81 = powers(3, [2, 0, 0, 1, 1], 4)          # beta^j  in F_3[y]/(y^4+y^3+2)
Z2, Z3 = (0,) * 4, (0,) * 4
SHAPE = (2,) * 4 + (3,) * 4


def load(name):
    L = [l.split() for l in open(os.path.join(CONN, name + ".txt")) if l.strip() and not l.startswith("#")]
    prm = tuple(map(int, L[0]))
    D = set()
    for i, j in L[1:]:
        a = Z2 if i == "-" else A16[int(i)]
        b = Z3 if j == "-" else B81[int(j)]
        D.add(a + b)
    assert len(D) == prm[1]
    return prm, D


def neg(x): return x[:4] + tuple((-c) % 3 for c in x[4:])


def char_values(S):
    """all character sums of the set S (multiset of the 1296 values, index 0 = trivial character)"""
    f = np.zeros(SHAPE)
    for x in S: f[x] = 1
    return np.fft.fftn(f).ravel()


def is_ls(S, m):
    if Z2 + Z3 in S or any(neg(x) not in S for x in S): return False
    if len(S) != m * 35: return False
    v = char_values(S)[1:]
    r, s = 36 - m, -m
    return bool(np.all(np.minimum(np.abs(v - r), np.abs(v - s)) < 1e-6))


def subspaces(p):
    """all 2-dim subspaces of F_p^4, as frozensets of vectors"""
    vecs = [v for v in itertools.product(range(p), repeat=4) if any(v)]
    seen = set(); out = []
    for u, w in itertools.combinations(vecs, 2):
        S = frozenset(tuple((a * x + b * y) % p for x, y in zip(u, w)) for a in range(p) for b in range(p))
        if len(S) == p * p and S not in seen: seen.add(S); out.append(S)
    return out


SUB2, SUB3 = subspaces(2), subspaces(3)
assert (len(SUB2), len(SUB3)) == (35, 130)


def group(i, j): return frozenset(a + b for a in SUB2[i] for b in SUB3[j])


def orthogonal(U, V):  # U = (i2, i3)
    return len(SUB2[U[0]] & SUB2[V[0]]) == 1 and len(SUB3[U[1]] & SUB3[V[1]]) == 1


def bush_check(E, U, V, W):
    """K = 2(A_E+I)-J indexed by (U-coset, V-coset); Bush type, unit-diagonal blocks, W a third orthogonal spread."""
    Ug, Vg = group(*U), group(*V)
    Ul, Vl = sorted(Ug), sorted(Vg)
    def add(x, y): return tuple((a + b) % 2 for a, b in zip(x[:4], y[:4])) + tuple((a + b) % 3 for a, b in zip(x[4:], y[4:]))
    # G = U (+) V: vertex (a, b) = u_a + v_b
    verts = [add(u, v) for u in Ul for v in Vl]
    assert len(set(verts)) == 1296
    idx = {x: t for t, x in enumerate(verts)}
    def sub(x, y): return add(x, neg(y))
    K = np.empty((1296, 1296), dtype=np.int64)
    for s, x in enumerate(verts):
        for t, y in enumerate(verts):
            K[s, t] = 1 if (s == t or sub(x, y) in E) else -1
    ok_had = np.array_equal(K @ K, 1296 * np.eye(1296, dtype=np.int64)) and np.array_equal(K, K.T)
    B = K.reshape(36, 36, 36, 36)                      # B[a, b, c, d] = K[(a,b),(c,d)]
    ok_diag = all((B[a, :, a, :] == 1).all() for a in range(36))
    ok_off = all((B[a, :, c, :].sum(0) == 0).all() and (B[a, :, c, :].sum(1) == 0).all()
                 for a in range(36) for c in range(36) if a != c)
    ok_unit = all((B[a, b, c, b] == 1) for a in range(36) for b in range(36) for c in range(36))
    res = dict(hadamard=ok_had, bush_diag=ok_diag, bush_offdiag_zero_sums=ok_off, unit_diagonal_blocks=bool(ok_unit))
    if W is not None:
        Wg = group(*W)
        # W-cosets are cliques of Cay(G,E) (W* in E) meeting every U-coset and every V-coset once
        res["third_spread_W"] = bool(all(w == Z2 + Z3 or w in E for w in Wg) and orthogonal(W, U) and orthogonal(W, V))
    return res


def main():
    out = {}
    names = sorted(f[:-4] for f in os.listdir(CONN) if f.endswith(".txt"))
    for name in names:
        prm, D = load(name)
        m = {560: 16, 592: None, 595: 17}[prm[1]]
        if m is None:
            assert not any(True for _ in [])
        free = [(i, j) for i in range(35) for j in range(130)
                if not any(x in D for x in group(i, j) if x != Z2 + Z3)]
        res = {"k": prm[1], "free_subgroups": len(free)}
        if m is not None and free:
            orth_edges = [(a, b) for a, b in itertools.combinations(range(len(free)), 2) if orthogonal(free[a], free[b])]
            res["orth_edges"] = orth_edges
            levels = {}
            for j in range(1, len(free) + 1):
                for T in itertools.combinations(range(len(free)), j):
                    if any(not orthogonal(free[a], free[b]) for a, b in itertools.combinations(T, 2)): continue
                    S = set(D)
                    for t in T: S |= {x for x in group(*free[t]) if x != Z2 + Z3}
                    good = is_ls(S, m + j)
                    levels.setdefault(j, [0, 0]); levels[j][0] += 1; levels[j][1] += int(good)
            res["orthogonal_sets_by_size"] = {j: v[0] for j, v in levels.items()}
            res["of_which_LS_PDS"] = {j: v[1] for j, v in levels.items()}
            # Bush check on one LS_19 configuration
            need = 19 - m
            for T in itertools.combinations(range(len(free)), need):
                if all(orthogonal(free[a], free[b]) for a, b in itertools.combinations(T, 2)):
                    E = set(D)
                    for t in T: E |= {x for x in group(*free[t]) if x != Z2 + Z3}
                    if m == 16: U, V, W = (free[t] for t in T)
                    else:   # m = 17: only two added; the third spread must already be in D -> search for it
                        U, V = (free[t] for t in T); W = None
                        for i in range(35):
                            for jj in range(130):
                                Wc = (i, jj)
                                if orthogonal(Wc, U) and orthogonal(Wc, V) and all(x == Z2 + Z3 or x in E for x in group(*Wc)):
                                    W = Wc; break
                            if W: break
                    res["bush_LS19"] = bush_check(E, U, V, W)
                    if W is None: res["bush_LS19"]["third_spread_W"] = "none among subgroups (D has no subgroup of order 36)"
                    break
        out[name] = res
        print(name, json.dumps({k: v for k, v in res.items() if k != "orth_edges"}), flush=True)
    if len(sys.argv) > 1: json.dump(out, open(sys.argv[1], "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
