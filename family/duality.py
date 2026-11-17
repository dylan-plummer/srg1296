"""Characters of G = F_16 x F_81 via the trace forms, and the dual set of a PDS.
chi_(c,d)(a,b) = (-1)^Tr_{16/2}(c a) * w^Tr_{81/3}(d b),  w = exp(2 pi i/3)."""
from collections import Counter
from fields import *

ELA = [Z4] + A16; ELB = [Z4] + B81


def mul(x, y, p, poly):
    r = [0] * 7
    for i in range(4):
        for j in range(4): r[i + j] = (r[i + j] + x[i] * y[j]) % p
    for d in range(6, 3, -1):
        c = r[d]
        if c:
            r[d] = 0
            for t in range(4): r[d - 4 + t] = (r[d - 4 + t] - c * poly[t]) % p
    return tuple(r[:4])


def trace(x, p, poly):
    s, t = x, x
    for _ in range(3):
        y = t
        for _ in range(p - 1): y = mul(y, t, p, poly)
        t = y; s = tuple((u + v) % p for u, v in zip(s, t))
    assert all(c == 0 for c in s[1:])
    return s[0]


P16, P81 = [1, 0, 0, 1], [2, 0, 0, 1]
TRA = {(c, a): trace(mul(c, a, 2, P16), 2, P16) for c in ELA for a in ELA}
TRB = {(d, b): trace(mul(d, b, 3, P81), 3, P81) for d in ELB for b in ELB}


def char_values(D):
    """exact character values: chi(D) = sum (-1)^tA * w^tB, returned as integers (2*real part / 2, checked)."""
    out = {}
    for c in ELA:
        sa = {a: (-1) ** TRA[(c, a)] for a in ELA}
        for d in ELB:
            n = Counter()
            for a, b in D: n[(sa[a], TRB[(d, b)])] += 1
            # sum of s * w^t : real part = sum_s s*(n[s,0] - (n[s,1]+n[s,2])/2), imag part ~ (n[s,1]-n[s,2])
            im = sum(s * (n[(s, 1)] - n[(s, 2)]) for s in (1, -1))
            re2 = sum(s * (2 * n[(s, 0)] - n[(s, 1)] - n[(s, 2)]) for s in (1, -1))
            assert im == 0 and re2 % 2 == 0
            out[(c, d)] = re2 // 2
    return out


def dual_set(D, k):
    chi = char_values(D); ev = Counter(v for x, v in chi.items() if x != (Z4, Z4))
    theta = next(v for v, n in ev.items() if n == k)
    return {x for x, v in chi.items() if x != (Z4, Z4) and v == theta}
