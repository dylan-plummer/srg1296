#!/usr/bin/env python3
"""CP-SAT search for rho-invariant PDSs in G = F_16 x F_81, rho = (alpha^12, beta^16), with D = -D, optionally with the
structure forced by Proposition 5.4 fixed in advance: T_0 = union of the index-8 classes J0, U_0 = union of the mu_5-cosets
I0, and |T_a| = big - 9 [a in U_0] for a != 0.  Exact integer character constraints (2*Re and Im/ (sqrt3/2) parts).
usage: rho_cpsat.py k lam mu seconds [J0 I0 big]   e.g.  rho_cpsat.py 245 64 42 600 0,4 0 18
Writes ../connsets-format file low_range/srg1296_k_lam_mu_rho.txt on success (checked by ../build_from_logs.py)."""
import sys, os, math, time
from ortools.sat.python import cp_model
from fields import *
from duality import TRA, TRB, ELA, ELB
k, lam, mu, secs = map(int, sys.argv[1:5])
fix = len(sys.argv) > 5
if fix:
    J0 = [int(x) for x in sys.argv[5].split(',')]; I0 = [int(x) for x in sys.argv[6].split(',')] if sys.argv[6] != '-' else []; big = int(sys.argv[7])
g = math.isqrt((lam - mu) ** 2 + 4 * (k - mu)); r, s = (lam - mu + g) // 2, (lam - mu - g) // 2
act = lambda x: (mulA(12, x[0]), mulB(16, x[1])); neg = lambda x: (x[0], mulB(40, x[1]))
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
if fix:
    for i, O in enumerate(orbits):
        a, b = O[0]
        if a == Z4: m.Add(X[i] == int(logB[b] % 8 in J0))
        elif b == Z4: m.Add(X[i] == int(logA[a] % 3 in I0))
    for a in A16:
        m.Add(sum(X[orbit_of[(a, b)]] for b in B81) + int(logA[a] % 3 in I0) == big - 9 * int(logA[a] % 3 in I0))
RE = {0: 2, 1: -1, 2: -1}; IM = {0: 0, 1: 1, 2: -1}
seen = set()
for c in ELA:
    for d in ELB:
        if (c, d) == (Z4, Z4): continue
        re = [0] * len(orbits); im = [0] * len(orbits)
        for i, O in enumerate(orbits):
            for a, b in O:
                sg = (-1) ** TRA[(c, a)]; t = TRB[(d, b)]
                re[i] += sg * RE[t]; im[i] += sg * IM[t]
        key = (tuple(re), tuple(im))
        if key in seen: continue
        seen.add(key); y = m.NewBoolVar('')
        m.Add(sum(cf * X[i] for i, cf in enumerate(re) if cf) == 2 * s + 2 * g * y)
        if any(im): m.Add(sum(cf * X[i] for i, cf in enumerate(im) if cf) == 0)
sv = cp_model.CpSolver(); sv.parameters.max_time_in_seconds = secs; sv.parameters.num_workers = int(os.environ.get("WORKERS", "4")); sv.parameters.random_seed = int(os.environ.get("SEED", "0"))
t0 = time.time(); st = sv.Solve(m)
print('RESULT', (1296, k, lam, mu), 'fixed' if fix else 'free', len(orbits), 'orbits', sv.StatusName(st), f'{time.time() - t0:.1f}s', flush=True)
if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
    D = {x for i, O in enumerate(orbits) if sv.Value(X[i]) for x in O}
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'low_range', f'srg1296_{k}_{lam}_{mu}_rho.txt')
    write_connset(out, (1296, k, lam, mu), D, f'rho-invariant PDS found by family/rho_cpsat.py ({" ".join(sys.argv[1:])})')
    print('wrote', out)
