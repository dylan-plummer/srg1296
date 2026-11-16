"""Independent check of srg(1296,420,144,132): D rebuilt from construction.json (explicit points + the two matrices);
connection set and matrix read as data only; none of the construction's code is used."""
import sys, itertools, json, numpy as np
R = sys.argv[1]; c = json.load(open(R + '/construction.json')); v_, k, lam, mu = c['parameters']
MB, MT = np.array(c['binary_matrix']), np.array(c['ternary_matrix'])
apB = lambda b: tuple(int(x) for x in MB @ np.array(b) % 2); apT = lambda t: tuple(int(x) for x in MT @ np.array(t) % 3)
assert all(np.array_equal(np.linalg.matrix_power(M, 3) % p, np.eye(4, dtype=int)) for M, p in ((MB, 2), (MT, 3)))
neg = lambda t: tuple((-x) % 3 for x in t)
T0 = {tuple(p) for p in c['fibre_projective_points']}; T0 |= {neg(t) for t in T0}
D = {((0,) * 4, t) for t in T0}
for s in c['seeds']:
    b = tuple((s['binary_label'] >> (3 - j)) & 1 for j in range(4))
    T = {tuple(p) for p in s['projective_points']}; T |= {neg(t) for t in T}
    if s['include_zero']: T.add((0,) * 4)
    for _ in range(3):
        D |= {(b, t) for t in T}; b = apB(b); T = {apT(t) for t in T}
assert len(D) == k, len(D)
crt = lambda b, t: tuple((3 * x + 4 * y) % 6 for x, y in zip(b, t))
Z = {crt(b, t) for b, t in D}; assert len(Z) == k and (0,) * 4 not in Z and all(tuple((-x) % 6 for x in z) in Z for z in Z)
cnt = {}
for d in Z:
    for e in Z:
        g = tuple((x - y) % 6 for x, y in zip(d, e)); cnt[g] = cnt.get(g, 0) + 1
print('exact PDS identity (all %d^2 differences):' % k, all(cnt.get(g, 0) == (k if g == (0,) * 4 else lam if g in Z else mu) for g in itertools.product(range(6), repeat=4)))
I = np.zeros((6,) * 4)
for z in Z: I[z] = 1
ev = np.round(np.fft.fftn(I).real, 6); print('spectrum:', dict(zip(*[x.tolist() for x in np.unique(ev, return_counts=True)])))
Zf = {tuple(int(x) for x in l.split()) for l in open(R + '/connections_Z6_4.txt') if l.strip()}
print('connection file = rebuilt D:', Zf == Z)
A = np.loadtxt(R + '/srg1296_420_144_132.txt', dtype=np.int64)
idx = lambda b, t: 81 * int(''.join(map(str, b)), 2) + int(''.join(map(str, t)), 3)
G = np.zeros((1296, 1296), dtype=np.int64)
for b in itertools.product(range(2), repeat=4):
    for t in itertools.product(range(3), repeat=4):
        i = idx(b, t)
        for db, dt in D: G[i, idx(tuple((x + y) % 2 for x, y in zip(b, db)), tuple((x + y) % 3 for x, y in zip(t, dt)))] = 1
print('matrix = Cay(G, D):', bool((G == A).all()))
J = np.ones_like(A); Id = np.eye(1296, dtype=np.int64)
print('A^2 = %dI + %dA + %dJ:' % (k - mu, lam - mu, mu), bool((A @ A == (k - mu) * Id + (lam - mu) * A + mu * J).all()))
sizes = {}
for b, t in D: sizes[b] = sizes.get(b, 0) + 1
print('fibre sizes over F2^4:', dict(sorted(sizes.items())))
print('U0 = {b : (b,0) in D}:', sorted(b for b, t in D if t == (0,) * 4))
