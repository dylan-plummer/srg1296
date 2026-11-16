"""LS_12(36) = srg(1296,420,144,132) as a Cayley set in F4^2 x F3^4, invariant under H = diag(F4*,F4*) acting on A=F4^2
through phi: H -> GL(4,3), h1=diag(w,1) -> g1, h2=diag(1,w) -> g2.  Template (d,j,d',j') = (42,12,12,42):
T_0 = 42-set, 18-fibres on the axis {(x,0)}, 27-fibres elsewhere.  Fibre reps: T0 (<g1,g2>-inv), Te1 = T_(1,0) (g2-inv),
Te2 = T_(0,1) (g1-inv), T11 = T_(1,1) (free, from the complete 9-divisible catalogue).
usage: ls420.py PAIRS_FILE [start end]   (each line: 'g1 ; g2' as 16 ints each, row-major)"""
import sys, os, itertools, json, time, numpy as np
m_LS = int(os.environ.get('M_LS', '12')); s_ev = -m_LS
T0PTS = int(os.environ.get('T0PTS', '21')); E1PTS, E1Z = map(int, os.environ.get('E1', '9,0').split(',')); E2PTS, E2Z = map(int, os.environ.get('E2', '13,1').split(','))
import os; MAXTRIP = int(os.environ.get('MAXTRIP', '200000'))
V3 = list(itertools.product(range(3), repeat=4)); P = [v for v in V3 if any(v) and next(x for x in v if x) == 1]
pidx = {p: i for i, p in enumerate(P)}
def norm(v):
    c = next(x for x in v if x); return tuple((x * c) % 3 for x in v)
DOT0 = np.array([[sum(a*b for a, b in zip(p, q)) % 3 == 0 for q in P] for p in P])      # [beta, point]
def chars(mask, zero):  # 40-vector of character sums of T = {0?} u +-P(mask)
    pts = np.array([(mask >> j) & 1 for j in range(40)], dtype=np.int64)
    n = pts.sum(); return int(zero) + 3 * (DOT0 @ pts) - n
def matmul(M, N): return tuple(tuple(sum(M[i][k]*N[k][j] for k in range(4)) % 3 for j in range(4)) for i in range(4))
I4 = tuple(tuple(int(i == j) for j in range(4)) for i in range(4))
def perm_pts(M):        # action of M on projective points (columns vectors)
    return [pidx[norm(tuple(sum(M[i][k]*p[k] for k in range(4)) % 3 for i in range(4)))] for p in P]
def perm_dual(M):       # beta -> M^T beta
    return [pidx[norm(tuple(sum(M[k][i]*p[k] for k in range(4)) % 3 for i in range(4)))] for p in P]
def orbits(perms):
    seen = [-1]*40; out = []
    for i in range(40):
        if seen[i] >= 0: continue
        O = {i}; fr = [i]
        while fr:
            x = fr.pop()
            for pm in perms:
                if pm[x] not in O: O.add(pm[x]); fr.append(pm[x])
        for x in O: seen[x] = len(out)
        out.append(sorted(O))
    return out
def invariant_sets(perms, npts, zero, check):
    O = orbits(perms); res = []
    for r in range(len(O) + 1):
        for comb in itertools.combinations(range(len(O)), r):
            if sum(len(O[c]) for c in comb) != npts: continue
            mask = sum(1 << x for c in comb for x in O[c]); ch = chars(mask, zero)
            if check(ch): res.append((mask, ch))
    return res
# A = F4^2 as 4 bits; F4 element c0 + c1 w -> (c0,c1); multiplication by w: (c0,c1) -> (c1, c0+c1)
def wmul(c, i):
    for _ in range(i): c = (c[1], (c[0] + c[1]) % 2)
    return c
F4 = [(0, 0), (1, 0), (0, 1), (1, 1)]
abits = lambda x, y: x + y
Avec = [abits(x, y) for x in F4 for y in F4]; aid = {a: i for i, a in enumerate(Avec)}
H16 = np.array([[(-1) ** (sum(p*q for p, q in zip(al, a)) % 2) for a in Avec] for al in Avec], dtype=np.int64)
one, zer = (1, 0), (0, 0)
a0 = aid[zer + zer]
ax1 = [aid[wmul(one, i) + zer] for i in range(3)]            # T_(w^i,0) = g1^i Te1
ax2 = [aid[zer + wmul(one, j)] for j in range(3)]            # T_(0,w^j) = g2^j Te2
free = [(aid[wmul(one, i) + wmul(one, j)], i, j) for i in range(3) for j in range(3)]
cat = None
def load_catalogue():
    global cat
    masks = np.array([int(l.split()[1], 16) for l in open('cat_0.tsv') if l.startswith('27')], dtype=np.uint64)
    pts = ((masks[:, None] >> np.arange(40, dtype=np.uint64)) & np.uint64(1)).astype(np.int8)
    C = (1 + 3 * (pts.astype(np.int16) @ DOT0.T.astype(np.int16)) - 13)
    assert np.all(C % 9 == 0)
    cat = (masks, (C // 9).astype(np.int8))
def run(g1, g2, log):
    g1 = tuple(tuple(r) for r in g1); g2 = tuple(tuple(r) for r in g2)
    pw = lambda M, i: I4 if i == 0 else matmul(M, pw(M, i - 1))
    G = {(i, j): matmul(pw(g1, i), pw(g2, j)) for i in range(3) for j in range(3)}
    assert matmul(g1, g2) == matmul(g2, g1) and pw(g1, 3) == I4 and pw(g2, 3) == I4
    allowed0 = {-m_LS + 9*t for t in range(5)}
    D0 = invariant_sets([perm_pts(g1), perm_pts(g2)], T0PTS, False, lambda ch: set(ch.tolist()) <= allowed0)
    D1 = invariant_sets([perm_pts(g2)], E1PTS, bool(E1Z), lambda ch: np.all(ch % 9 == 0))
    D2 = invariant_sets([perm_pts(g1)], E2PTS, bool(E2Z), lambda ch: np.all(ch % 9 == 0))
    rec = dict(g1=g1, g2=g2, D0=len(D0), D1=len(D1), D2=len(D2))
    if not (D0 and D1 and D2): rec['triples'] = 0; log.write(json.dumps(rec) + '\n'); log.flush(); return []
    if len(D0) * len(D1) * len(D2) > MAXTRIP: rec['skipped_triples'] = len(D0) * len(D1) * len(D2); log.write(json.dumps(rec) + '\n'); log.flush(); return []
    dual = {ij: perm_dual(G[ij]) for ij in G}          # beta -> (phi(h)^T beta)
    masks, C = cat; found = []; ntrip = 0; t0 = time.time()
    # unknown-part functionals: for each beta, the 9 catalogue coordinates dual[(i,j)][beta] with sign (-1)^{alpha.a_ij}
    Hfree = H16[:, [a for a, _, _ in free]]                                 # 16 x 9
    cols = np.array([[dual[(i, j)][b] for _, i, j in free] for b in range(40)])   # 40 x 9
    # order betas: most constraining first (fewest distinct catalogue coordinates)
    border = sorted(range(40), key=lambda b: len(set(cols[b])))
    Hf16 = Hfree.T.astype(np.int16)
    U = [C[:, cols[b]].astype(np.int16) @ Hf16 for b in range(40)]          # 40 x (N x 16)
    UQ = [np.unique(u, axis=0, return_inverse=True) for u in U]
    for m0, c0 in D0:
        for m1, c1 in D1:
            for m2, c2 in D2:
                ntrip += 1
                Yk = np.zeros((40, 16), dtype=np.int64)
                Yk[:, a0] = c0 + m_LS
                for i in range(3): Yk[:, ax1[i]] = c1[dual[(i, 0)]]
                for j in range(3): Yk[:, ax2[j]] = c2[dual[(0, j)]]
                K = Yk @ H16.T
                if np.any(K % 9): continue
                Kq = K // 9
                mask = None
                for b in border[:4]:
                    uq, inv = UQ[b]; v = uq + Kq[b]
                    okr = np.all((v == 0) | (v == 4), axis=1)
                    if not okr.any(): mask = None; break
                    mb = okr[inv.ravel()]
                    mask = mb if mask is None else (mask & mb)
                    if not mask.any(): break
                if mask is None or not mask.any(): continue
                cand = np.flatnonzero(mask)
                for b in border[4:]:
                    v = U[b][cand] + Kq[b]
                    cand = cand[np.all((v == 0) | (v == 4), axis=1)]
                    if not len(cand): break
                if len(cand):
                    for ci in cand[:5]:
                        found.append(dict(T0=hex(m0), Te1=hex(m1), Te2=hex(m2), T11=hex(int(masks[ci]))))
    rec.update(triples=ntrip, found=len(found), seconds=round(time.time() - t0, 1))
    if found: rec['solutions'] = found
    log.write(json.dumps(rec) + '\n'); log.flush()
    return found
if __name__ == '__main__':
    pairs = [l.split(';') for l in open(sys.argv[1]) if l.strip()]
    lo, hi = (int(sys.argv[2]), int(sys.argv[3])) if len(sys.argv) > 3 else (0, len(pairs))
    load_catalogue()
    with open(os.environ.get('TAG', '') + f'log_{lo}_{hi}.jsonl', 'a') as log:
        for line in pairs[lo:hi]:
            g1, g2 = [np.array(list(map(int, x.split()))).reshape(4, 4).tolist() for x in line]
            if run(g1, g2, log): print('FOUND', g1, g2, flush=True)
    print('DONE', lo, hi, flush=True)
