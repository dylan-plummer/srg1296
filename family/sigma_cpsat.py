#!/usr/bin/env python3
"""CP-SAT in the model <sigma, -1>, sigma = (a -> alpha^5 a, b -> u(b)), u(v0,v1,v2,v3) = (v0+v2, v1+v3, v2, v3) on coefficient
vectors of F_81 (the symmetry of graphs VII and VIII), with the binary data forced by Theorem 2.7 fixed: U_0 = union of the
F_4-lines alpha^i F_4^* with i mod 5 in I0, ring sizes |T_a| = odd value on U_0, even value elsewhere.
usage: sigma_cpsat.py k lam mu seconds I0 odd even      e.g. sigma_cpsat.py 315 90 72 1200 0 27 18"""
import sys, os, math, time
from ortools.sat.python import cp_model
from fields import *
from duality import TRA, TRB, ELA, ELB
k, lam, mu, secs = map(int, sys.argv[1:5]); I0 = [int(x) for x in sys.argv[5].split(',')] if sys.argv[5] != '-' else []
odd, even = int(sys.argv[6]), int(sys.argv[7])
g = math.isqrt((lam - mu) ** 2 + 4 * (k - mu)); r, s = (lam - mu + g) // 2, (lam - mu - g) // 2
u = lambda b: ((b[0] + b[2]) % 3, (b[1] + b[3]) % 3, b[2], b[3])
act = lambda x: (mulA(5, x[0]), u(x[1])); neg = lambda x: (x[0], tuple((-c) % 3 for c in x[1]))
inU0 = lambda a: a != Z4 and logA[a] % 5 in I0
orbit_of = {}; orbits = []
for a in ELA:
    for b in ELB:
        x = (a, b)
        if x == (Z4, Z4) or x in orbit_of: continue
        O = set(); fr = [x]
        while fr:
            y = fr.pop()
            if y in O: continue
            O.add(y); fr += [act(y), neg(y)]
        for y in O: orbit_of[y] = len(orbits)
        orbits.append(sorted(O))
m = cp_model.CpModel(); X = [m.NewBoolVar('') for _ in orbits]
m.Add(sum(len(O) * X[i] for i, O in enumerate(orbits)) == k)
for i, O in enumerate(orbits):
    a, b = O[0]
    if b == Z4: m.Add(X[i] == int(inU0(a)))
for a in A16:
    m.Add(sum(X[orbit_of[(a, b)]] for b in B81) + int(inU0(a)) == (odd if inU0(a) else even))
RE = {0: 2, 1: -1, 2: -1}; IM = {0: 0, 1: 1, 2: -1}; seen = set()
for c in ELA:
    for d in ELB:
        if (c, d) == (Z4, Z4): continue
        re = [0] * len(orbits); im = [0] * len(orbits)
        for i, O in enumerate(orbits):
            for a, b in O:
                sg = (-1) ** TRA[(c, a)]; t = TRB[(d, b)]; re[i] += sg * RE[t]; im[i] += sg * IM[t]
        key = (tuple(re), tuple(im))
        if key in seen: continue
        seen.add(key); y = m.NewBoolVar('')
        m.Add(sum(cf * X[i] for i, cf in enumerate(re) if cf) == 2 * s + 2 * g * y)
        if any(im): m.Add(sum(cf * X[i] for i, cf in enumerate(im) if cf) == 0)
sv = cp_model.CpSolver(); sv.parameters.max_time_in_seconds = secs; sv.parameters.num_workers = int(os.environ.get('WORKERS', '4'))
t0 = time.time(); st = sv.Solve(m)
print('RESULT', (1296, k, lam, mu), len(orbits), 'orbits', sv.StatusName(st), f'{time.time() - t0:.1f}s', flush=True)
if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
    D = {x for i, O in enumerate(orbits) if sv.Value(X[i]) for x in O}
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'low_range', f'srg1296_{k}_{lam}_{mu}_sigma.txt')
    write_connset(out, (1296, k, lam, mu), D, f'sigma-invariant PDS found by family/sigma_cpsat.py ({" ".join(sys.argv[1:])})'); print('wrote', out)
