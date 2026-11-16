#!/usr/bin/env python3
"""Import a degree-420 or degree-385 connection set (graphs VII, VIII) (Z/6)^4 -> F_16 x F_81.  Its order-3 symmetry is R = (B, U) with
B(b) = (b1, b0+b1, b3, b2+b3) on F_2^4 (fixed-point-free) and U(v) = (v0+v2, v1+v3, v2, v3) on F_3^4 (unipotent, Jordan 2+2).
phi_A: B^k e_i -> alpha^(5k) c_i (e_1 = (1,0,0,0), c_1 = 1; e_2 = (0,0,1,0), c_2 = alpha), so phi_A B = alpha^5 phi_A;
phi_B: coefficient vectors (v0..v3) -> v0 + v1 y + v2 y^2 + v3 y^3 in F_81.
usage: import_c3.py CONNECTIONS_Z6_4.txt k lam mu NAME OUT.txt"""
import sys, itertools
from fields import A16, B81, Z4, mulA, write_connset
Bm = lambda b: (b[1], (b[0] + b[1]) % 2, b[3], (b[2] + b[3]) % 2)
basis = [(1, 0, 0, 0), Bm((1, 0, 0, 0)), (0, 0, 1, 0), Bm((0, 0, 1, 0))]
img = [A16[0], A16[5], A16[1], A16[6]]
phiA = {}
for c in itertools.product(range(2), repeat=4):
    x = tuple(sum(ci * bk[t] for ci, bk in zip(c, basis)) % 2 for t in range(4))
    phiA[x] = tuple(sum(ci * ik[t] for ci, ik in zip(c, img)) % 2 for t in range(4))
assert len(set(phiA.values())) == 16 and all(phiA[Bm(x)] == mulA(5, phiA[x]) for x in phiA)
D = set()
for l in open(sys.argv[1]):
    if not l.strip(): continue
    z = tuple(int(x) for x in l.split())
    D.add((phiA[tuple(x % 2 for x in z)], tuple(x % 3 for x in z)))
k, lam, mu, name, out = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5], sys.argv[6]
assert len(D) == k and (Z4, Z4) not in D
write_connset(out, (1296, k, lam, mu), D,
              f'graph {name}: fibre-join construction with a C3 symmetry (16 Nov 2026), imported by family/import_c3.py\n'
              'invariant under (a, b) -> (alpha^5 a, U b), U(v0,v1,v2,v3) = (v0+v2, v1+v3, v2, v3) on coefficient vectors')
print('wrote', out)
