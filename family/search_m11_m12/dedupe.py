import itertools, numpy as np, random
exec(open('pairs.py').read().split("out = []; seen = set()")[0])
random.seed(0)
pairs = [tuple(np.array(list(map(int, x.split()))).reshape(4, 4) for x in l.split(';')) for l in open('pairs_all.txt')]
V = np.array(list(itertools.product(range(3), repeat=4))).T
P = [v for v in itertools.product(range(3), repeat=4) if any(v) and next(x for x in v if x) == 1]
pid = {p: i for i, p in enumerate(P)}
def norm(v): c = next(x for x in v if x); return tuple((x*c) % 3 for x in v)
def norbits(gens):
    seen = set(); n = 0
    for p in P:
        if p in seen: continue
        n += 1; fr = [p]; seen.add(p)
        while fr:
            x = fr.pop()
            for g in gens:
                y = norm(tuple(int(t) for t in (g @ np.array(x)) % 3))
                if y not in seen: seen.add(y); fr.append(y)
    return n
def inv3(M):
    # brute inverse via powers: M^(order-1); order divides |GL| ; use adjugate mod 3
    d = int(round(np.linalg.det(M))) % 3
    adj = np.round(np.linalg.inv(M) * np.linalg.det(M)).astype(np.int64) % 3
    return (adj * d) % 3   # d^{-1} = d in F3
def centraliser_gens(g, n=25):
    E = [np.zeros((4,4), dtype=np.int64) for _ in range(16)]
    for i in range(16): E[i].flat[i] = 1
    basis = []
    L = np.array([((g @ e - e @ g) % 3).ravel() for e in E]).T % 3
    A = L.copy(); piv = []; r = 0
    for c in range(16):
        p = next((i for i in range(r, 16) if A[i, c]), None)
        if p is None: continue
        A[[r, p]] = A[[p, r]]; A[r] = (A[r] * A[r, c]) % 3
        for i in range(16):
            if i != r and A[i, c]: A[i] = (A[i] - A[i, c] * A[r]) % 3
        piv.append(c); r += 1
    for f in [c for c in range(16) if c not in piv]:
        x = np.zeros(16, dtype=np.int64); x[f] = 1
        for i, c in enumerate(piv): x[c] = (-A[i, f]) % 3
        basis.append(x.reshape(4, 4))
    gens = []
    while len(gens) < n:
        X = sum(random.randrange(3) * b for b in basis) % 3
        if int(round(np.linalg.det(X))) % 3: gens.append(X)
    return gens
key = lambda a, b: (tuple(a.ravel()), tuple(b.ravel()))
keep = []; done = set()
cg = {}
for a, b in pairs:
    k = key(a, b)
    if k in done: continue
    # conjugation by C(a) ∩ C(b) is not needed; we use C(a) if a canonical (non-identity), else C(b)
    base = a if not np.array_equal(a, np.eye(4, dtype=np.int64)) and k[0] in [tuple(v.ravel()) for v in canon.values()] else b
    gk = tuple(base.ravel())
    if gk not in cg: cg[gk] = [(h, inv3(h)) for h in centraliser_gens(base)]
    orb = {k}; fr = [(a, b)]
    while fr:
        x, y = fr.pop()
        for h, hi in cg[gk]:
            nx, ny = (h @ x @ hi) % 3, (h @ y @ hi) % 3; kk = key(nx, ny)
            if kk not in orb: orb.add(kk); fr.append((nx, ny))
    done |= orb
    n12, n1, n2 = norbits([a, b]), norbits([a]), norbits([b])
    keep.append((a, b, n12, n1, n2, len(orb)))
print('classes', len(keep))
ok = [t for t in keep if max(t[2], t[3], t[4]) <= 22]
print('with all orbit counts <= 22:', len(ok))
with open('pairs_cls.txt', 'w') as f:
    for a, b, *_ in sorted(ok, key=lambda t: t[2] + t[3] + t[4]):
        f.write(' '.join(map(str, a.ravel())) + ' ; ' + ' '.join(map(str, b.ravel())) + '\n')
import collections; print(collections.Counter((t[2], t[3], t[4]) for t in ok))
