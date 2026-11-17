"""CP-SAT for C5-invariant Cayley SRGs on F2^4 x F3^4, with/without the divisibility cuts forced by r-s.
usage: divcut.py k lam mu cuts(0/1) seconds [workers] [model]  ; model = c5 (default) | c5only_b (C5 on F3^4 only)"""
import sys, itertools, math, time, numpy as np
from ortools.sat.python import cp_model
k, lam, mu, cuts, secs = map(int, sys.argv[1:6]); workers = int(sys.argv[6]) if len(sys.argv) > 6 else 8
model_name = sys.argv[7] if len(sys.argv) > 7 else 'c5'
gap = math.isqrt((lam-mu)**2 + 4*(k-mu)); r, s = (lam-mu+gap)//2, (lam-mu-gap)//2; Dl = r - s
A = list(itertools.product(range(2), repeat=4)); B = list(itertools.product(range(3), repeat=4)); O2, O3 = A[0], B[0]
Bm = lambda b: ((b[0]+b[1]+b[2]+b[3]) % 2, b[0], b[1], b[2]); Fm = lambda v: ((-v[0]-v[1]-v[2]-v[3]) % 3, v[0], v[1], v[2])
neg = lambda v: tuple((-x) % 3 for x in v)
G = [(a, v) for a in A for v in B if (a, v) != (O2, O3)]
def act(g):
    a, v = g
    return (Bm(a), Fm(v)) if model_name == 'c5' else (a, Fm(v))
orb = {}; reps = []
for g in G:
    if g in orb: continue
    O = set(); fr = [g]
    while fr:
        x = fr.pop()
        if x in O: continue
        O.add(x); fr += [act(x), (x[0], neg(x[1]))]
    for x in O: orb[x] = len(reps)
    reps.append(sorted(O))
n = len(reps); print('orbits', n, 'r,s', r, s, flush=True)
m = cp_model.CpModel(); X = [m.NewBoolVar(f'x{i}') for i in range(n)]
x = lambda a, v: 0 if (a, v) == (O2, O3) else X[orb[(a, v)]]
c3 = lambda b, v: 2 if sum(p*q for p, q in zip(b, v)) % 3 == 0 else -1     # 2*cos
c2 = lambda al, a: -1 if sum(p*q for p, q in zip(al, a)) % 2 else 1
m.Add(sum(len(O) * X[i] for i, O in enumerate(reps)) == k)
# character constraints, one per orbit of characters (characters orbits = same structure under dual action; use all, dedup by coefficient vector)
seen = set(); ncon = 0
for al in A:
    for be in B:
        if (al, be) == (O2, O3): continue
        coef = [0]*n
        for i, O in enumerate(reps):
            coef[i] = sum(c2(al, a) * c3(be, v) for a, v in O)
        key = tuple(coef)
        if key in seen: continue
        seen.add(key); y = m.NewBoolVar('')
        m.Add(sum(c * X[i] for i, c in enumerate(coef) if c) == 2*s + 2*Dl*y); ncon += 1
print('character constraints', ncon, flush=True)
if cuts:
    g1 = Dl // math.gcd(16, Dl); g2 = Dl // math.gcd(81, Dl); nc = 0; seen = set()
    proj = [v for v in B if v != O3 and next(t for t in v if t) == 1]
    for a in A:                                  # f_beta(a) = sum_v x(a,v) cos  in g1*Z (a != 0), = s mod g1 (a = 0)
        for be in proj:
            terms = [(c3(be, v), x(a, v)) for v in B]
            key = ('A', a, tuple(orb.get((a, v), -1) for v in B), be)
            z = m.NewIntVar(-200, 200, '')
            off = 0 if a != O2 else 2*s
            m.Add(sum(c * t for c, t in terms if not isinstance(t, int)) == off + 2*g1*z); nc += 1
    for v in [O3] + proj:                        # Walsh(U_v)(alpha) in g2*Z (v != 0), = s mod g2 (v = 0, alpha != 0)
        for al in A:
            if v == O3 and al == O2: continue
            z = m.NewIntVar(-40, 40, ''); off = 0 if v != O3 else s
            m.Add(sum(c2(al, a) * x(a, v) for a in A if not isinstance(x(a, v), int)) == off + g2*z); nc += 1
    if Dl % 36 == 0 and (k - s) % Dl == 0:      # LS-like: beta = 0 rows and alpha = 0 columns too
        for a in A[1:]:
            z = m.NewIntVar(0, 20, ''); m.Add(sum(x(a, v) for v in B) == g1*z); nc += 1
    print('cuts', nc, 'g1', g1, 'g2', g2, flush=True)
sv = cp_model.CpSolver(); sv.parameters.max_time_in_seconds = secs; sv.parameters.num_workers = workers
t = time.time(); st = sv.Solve(m)
print('RESULT', (k, lam, mu), 'cuts', cuts, 'model', model_name, sv.StatusName(st), f'{time.time()-t:.1f}s', flush=True)
if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
    D = [g for i, O in enumerate(reps) if sv.Value(X[i]) for g in O]
    crt = lambda a, v: tuple((3*p + 4*q) % 6 for p, q in zip(a, v))
    Z = {crt(*g) for g in D}; I = np.zeros((6,)*4)
    for z in Z: I[z] = 1
    ev = np.round(np.fft.fftn(I).real, 6); print('check spectrum', dict(zip(*np.unique(ev, return_counts=True))))
    open(f'/mnt/jinstore/JinLab02/dmp131/paut_runs/divcut/sol_{k}_{lam}_{mu}_{model_name}.txt', 'w').write(''.join(' '.join(map(str, z)) + '\n' for z in sorted(Z)))
