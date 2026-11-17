"""Shared field arithmetic for the family scripts: F_16 = F_2[x]/(x^4+x^3+1), alpha = x; F_81 = F_3[y]/(y^4+y^3+2), beta = y.
Elements are coefficient tuples (c_0, c_1, c_2, c_3), lowest degree first, as in ../build_from_logs.py."""


def powers(p, poly, n):
    cur = [1] + [0] * (n - 1); out = []
    for _ in range(p ** n - 1):
        out.append(tuple(cur))
        top = cur[-1]; cur = [0] + cur[:-1]
        cur = [(c - top * poly[i]) % p for i, c in enumerate(cur)]
    assert len(set(out)) == p ** n - 1
    return out


A16 = powers(2, [1, 0, 0, 1, 1], 4)          # A16[i] = alpha^i
B81 = powers(3, [2, 0, 0, 1, 1], 4)          # B81[j] = beta^j
Z4 = (0, 0, 0, 0)
logA = {a: i for i, a in enumerate(A16)}
logB = {b: j for j, b in enumerate(B81)}


def mulA(i, a):   # alpha^i * a
    return a if a == Z4 else A16[(logA[a] + i) % 15]


def mulB(j, b):
    return b if b == Z4 else B81[(logB[b] + j) % 80]


def frobA(u, a):  # a^(2^u)
    return a if a == Z4 else A16[(logA[a] * 2 ** u) % 15]


def frobB(w, b):  # b^(3^w)
    return b if b == Z4 else B81[(logB[b] * 3 ** w) % 80]


def read_connset(path):
    """connsets/*.txt -> (params, set of (a, b) coefficient-tuple pairs)."""
    lines = [l.split() for l in open(path) if l.strip() and not l.startswith('#')]
    params = tuple(map(int, lines[0])); D = set()
    for i, j in lines[1:]:
        D.add((Z4 if i == '-' else A16[int(i) % 15], Z4 if j == '-' else B81[int(j) % 80]))
    assert len(D) == params[1]
    return params, D


def write_connset(path, params, D, comment):
    key = lambda e: (e[0] != '-', int(e[0]) if e[0] != '-' else -1, int(e[1]) if e[1] != '-' else -1)
    rows = sorted([('-' if a == Z4 else str(logA[a]), '-' if b == Z4 else str(logB[b])) for a, b in D], key=key)
    with open(path, 'w', newline='\n') as f:
        f.write(' '.join(map(str, params)) + '\n')
        for c in comment.splitlines(): f.write('# ' + c + '\n')
        for i, j in rows: f.write(f'{i} {j}\n')
