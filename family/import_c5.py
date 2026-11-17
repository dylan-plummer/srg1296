#!/usr/bin/env python3
"""Import a C5 'fibre-join' connection set given in (Z/6)^4 coordinates (z_i = 3 b_i + 4 v_i mod 6, b in F_2^4,
v in F_3^4) into the F_16 x F_81 coordinates of this repository.

The construction's order-5 map is R = (B, F) with B(b) = (b1+b2+b3+b4, b1, b2, b3) on F_2^4 and
F(v) = (-v1-v2-v3-v4, v1, v2, v3) on F_3^4. Both have characteristic polynomial x^4+x^3+x^2+x+1, irreducible over
F_2 and F_3, so for any nonzero e the maps  B^k e -> zeta^k,  F^k e' -> eta^k  (k = 0..3) extend to linear
isomorphisms phi_A : F_2^4 -> F_16, phi_B : F_3^4 -> F_81 with phi_A B = zeta phi_A and phi_B F = eta phi_B, where
zeta = alpha^12 and eta = beta^16 are the elements of order 5 of the model <(alpha^12, beta^16)> of the paper's search.
usage: import_c5.py CONNECTIONS_Z6_4.txt k lam mu OUT.txt"""
import sys, itertools
from fields import A16, B81, Z4, mulA, mulB, write_connset

src, k, lam, mu, out = sys.argv[1], *map(int, sys.argv[2:5]), sys.argv[5]
Bm = lambda b: ((b[0] + b[1] + b[2] + b[3]) % 2, b[0], b[1], b[2])
Fm = lambda v: ((-v[0] - v[1] - v[2] - v[3]) % 3, v[0], v[1], v[2])


def iso(p, step, e, mul_zeta):
    """linear iso phi with phi(step^k e) = zeta^k (as field elements), checked by enumeration."""
    basis = [e]
    for _ in range(3): basis.append(step(basis[-1]))
    img = [mul_zeta(k, (1, 0, 0, 0)) for k in range(4)]
    phi = {}
    for c in itertools.product(range(p), repeat=4):
        x = tuple(sum(ci * bk[t] for ci, bk in zip(c, basis)) % p for t in range(4))
        y = tuple(sum(ci * ik[t] for ci, ik in zip(c, img)) % p for t in range(4))
        phi[x] = y
    assert len(phi) == p ** 4 and len(set(phi.values())) == p ** 4
    assert all(phi[step(x)] == mul_zeta(1, phi[x]) for x in phi)
    return phi


phiA = iso(2, Bm, (0, 0, 0, 1), lambda n, a: mulA(12 * n, a))
phiB = iso(3, Fm, (0, 0, 0, 1), lambda n, b: mulB(16 * n, b))
D = set()
for l in open(src):
    if not l.strip() or l.startswith('#'): continue
    z = tuple(int(x) for x in l.split())
    D.add((phiA[tuple(x % 2 for x in z)], phiB[tuple(x % 3 for x in z)]))
assert len(D) == k and (Z4, Z4) not in D
assert all((a, mulB(40, b)) in D for a, b in D)                       # D = -D
assert all((mulA(12, a), mulB(16, b)) in D for a, b in D)            # invariant under (alpha^12, beta^16)
write_connset(out, (1296, k, lam, mu), D,
              f'imported by family/import_c5.py from a (Z/6)^4 connection set; invariant under (alpha^12, beta^16) and -1')
print('wrote', out, len(D))
