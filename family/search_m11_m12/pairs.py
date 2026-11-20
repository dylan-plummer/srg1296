import itertools, numpy as np
def mm(A, B): return (A @ B) % 3
I = np.eye(4, dtype=np.int64)
def J(blocks):
    M = np.zeros((4, 4), dtype=np.int64); i = 0
    for b in blocks:
        for t in range(b):
            M[i+t, i+t] = 1
            if t + 1 < b: M[i+t, i+t+1] = 1
        i += b
    return M
canon = {'J2': J([2,1,1]), 'J22': J([2,2]), 'J31': J([3,1])}
def det3(M): return int(round(np.linalg.det(M))) % 3
def key(M): return tuple(M.ravel().tolist())
out = []; seen = set()
for name, g1 in canon.items():
    # centralizer in M4(F3): solve g1 X - X g1 = 0 by brute force on a basis (linear algebra mod 3)
    basis = []
    E = [np.zeros((4,4), dtype=np.int64) for _ in range(16)]
    for i in range(16): E[i].flat[i] = 1
    # build linear map X -> g1 X - X g1 as 16x16 matrix mod 3, nullspace by Gaussian elimination
    L = np.array([((g1 @ e - e @ g1) % 3).ravel() for e in E]).T % 3
    # nullspace mod 3
    A = L.copy(); piv = []; r = 0
    for c in range(16):
        p = next((i for i in range(r, 16) if A[i, c]), None)
        if p is None: continue
        A[[r, p]] = A[[p, r]]; A[r] = (A[r] * A[r, c]) % 3    # inverse of 1 is 1, of 2 is 2
        for i in range(16):
            if i != r and A[i, c]: A[i] = (A[i] - A[i, c] * A[r]) % 3
        piv.append(c); r += 1
    free = [c for c in range(16) if c not in piv]
    for f in free:
        x = np.zeros(16, dtype=np.int64); x[f] = 1
        for i, c in enumerate(piv): x[c] = (-A[i, f]) % 3
        basis.append(x.reshape(4, 4))
    cnt = 0
    for coeffs in itertools.product(range(3), repeat=len(basis)):
        X = sum(c * b for c, b in zip(coeffs, basis)) % 3 if basis else np.zeros((4,4), dtype=np.int64)
        if det3(X) == 0: continue
        if not np.array_equal(mm(mm(X, X), X), I): continue
        for pair in [(g1, X), (X, g1)]:
            k = (key(pair[0]), key(pair[1]))
            if k not in seen: seen.add(k); out.append(pair); cnt += 1
    print(name, 'centralizer dim', len(basis), 'pairs', cnt)
for g in canon.values():
    for pair in [(I, g), (g, I)]:
        k = (key(pair[0]), key(pair[1]))
        if k not in seen: seen.add(k); out.append(pair)
out.insert(0, (I, I))
with open('pairs_all.txt', 'w') as f:
    for a, b in out: f.write(' '.join(map(str, a.ravel())) + ' ; ' + ' '.join(map(str, b.ravel())) + '\n')
print('total', len(out))
