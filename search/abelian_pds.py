"""Exact search for (v,k,lam,mu) partial difference sets in an abelian group R = Z_{n1} x ... x Z_{nt}
that are unions of orbits of a group K of automorphisms containing all multipliers x -> t x,
gcd(t, v) = 1 (necessary for integral-eigenvalue PDS by Ma's multiplier theorem).

Condition (character criterion): chi(D) in {r, s} for every non-principal character chi, |D| = k,
0 not in D (D = -D is automatic as -1 is a multiplier).

usage as a library:  search(shape, k, lam, mu, extra_auts=[], time_limit=..)
extra_auts: list of integer matrices M (t x t) acting x -> M x (mod shape); must be automorphisms.
"""
import itertools, math, sys, time
import numpy as np
from ortools.sat.python import cp_model


def abelian_groups(v):
    """All abelian groups of order v as lists of prime-power cyclic factors."""
    fac = {}
    n = v; d = 2
    while d * d <= n:
        while n % d == 0: fac[d] = fac.get(d, 0) + 1; n //= d
        d += 1
    if n > 1: fac[n] = fac.get(n, 0) + 1
    def partitions(a, mx=None):
        if mx is None: mx = a
        if a == 0: yield []; return
        for i in range(min(a, mx), 0, -1):
            for p in partitions(a - i, i): yield [i] + p
    per = [[[p ** e for e in part] for part in partitions(a)] for p, a in fac.items()]
    for combo in itertools.product(*per):
        yield sorted(sum(combo, []))


class AbelianPDS:
    def __init__(self, shape, extra_auts=()):
        self.shape = tuple(shape)
        self.v = v = int(np.prod(shape))
        self.e = e = int(np.lcm.reduce(shape))
        grid = np.indices(shape).reshape(len(shape), -1).T  # v x t, row-major order
        self.grid = grid
        mult = np.array(shape)
        def index(X):
            return np.ravel_multi_index(tuple((X % mult).T), shape)
        # generators of K: multipliers (units mod e) + extra automorphisms
        perms = []
        units = [t for t in range(1, e) if math.gcd(t, e) == 1] if e > 1 else [1]
        # generate unit group by a few generators is enough, but using all is cheap
        for t in units:
            perms.append(index(grid * t))
        for M in extra_auts:
            M = np.array(M)
            perms.append(index(grid @ M.T))
        # union-find orbits on R
        parent = np.arange(v)
        def find(a):
            while parent[a] != a:
                parent[a] = parent[parent[a]]; a = parent[a]
            return a
        for p in perms:
            for x in range(v):
                a, b = find(x), find(int(p[x]))
                if a != b: parent[max(a, b)] = min(a, b)
        lab = np.array([find(x) for x in range(v)])
        reps, orb = np.unique(lab, return_inverse=True)
        self.orb = orb
        self.norb = len(reps)
        self.osize = np.bincount(orb)
        self.zero_orb = orb[0]
        # dual action: characters a -> chi_a(x) = exp(2 pi i sum a_j x_j / n_j). Automorphism M on R
        # induces a permutation of characters; orbits of the dual action = orbits of K^T on R^.
        dperms = [index(grid * t) for t in units]
        for M in extra_auts:
            M = np.array(M)
            # chi_a(M x) = chi_{a'}(x) with a' = (diag(n)^-1-scaled) transpose; for each factor size
            # use the general formula a'_j = sum_i a_i M_ij * (n_j / n_i)  (valid for automorphisms)
            n = np.array(shape)
            W = (M * (n[None, :] / n[:, None]))  # W_ij = M_ij * n_j / n_i
            Wi = np.rint(W).astype(np.int64)
            assert np.allclose(W, Wi), "automorphism not compatible"
            dperms.append(index(grid @ Wi))
        parent = np.arange(v)
        for p in dperms:
            for x in range(v):
                a, b = find(x), find(int(p[x]))
                if a != b: parent[max(a, b)] = min(a, b)
        dlab = np.array([find(x) for x in range(v)])
        dreps, dorb = np.unique(dlab, return_inverse=True)
        self.dorb = dorb
        self.ndorb = len(dreps)
        self.dsize = np.bincount(dorb)
        self.drep = dreps  # representative character index for each dual orbit (min element)
        # character sums: C[j, i] = sum_{x in O_i} chi_{a_j}(x)
        C = np.zeros((self.ndorb, self.norb))
        for i in range(self.norb):
            f = (orb == i).astype(float).reshape(shape)
            F = np.fft.fftn(f)          # F[a] = sum_x f(x) exp(-2 pi i <a,x>)
            Fa = F.reshape(-1)[dreps]
            C[:, i] = Fa.real
        Ci = np.rint(C).astype(np.int64)
        if np.abs(C - Ci).max() > 1e-6:
            raise ValueError("non-integral character sums (K too small)")
        self.C = Ci

    def search(self, k, lam, mu, time_limit=60.0, workers=4, max_solutions=1, forbid_full=True, solver_params=None):
        D = (lam - mu) ** 2 + 4 * (k - mu)
        sq = math.isqrt(D)
        if sq * sq != D:
            return [], "IRRATIONAL"
        r = ((lam - mu) + sq) // 2; s = ((lam - mu) - sq) // 2
        fmult = (-k - s * (self.v - 1)) // (r - s)
        m = cp_model.CpModel()
        x = [m.NewBoolVar(f"x{i}") for i in range(self.norb)]
        m.Add(x[self.zero_orb] == 0)
        m.Add(sum(int(self.osize[i]) * x[i] for i in range(self.norb)) == k)
        trivial_dual = self.dorb[0]
        bs = []
        for j in range(self.ndorb):
            if j == trivial_dual: continue
            b = m.NewBoolVar("")
            m.Add(sum(int(self.C[j, i]) * x[i] for i in range(self.norb) if self.C[j, i]) == s + (r - s) * b)
            bs.append((int(self.dsize[j]), b))
        m.Add(sum(c * b for c, b in bs) == fmult)
        solver = cp_model.CpSolver()
        solver.parameters.max_time_in_seconds = time_limit
        solver.parameters.num_workers = workers
        for key, val in (solver_params or {}).items():
            setattr(solver.parameters, key, val)
        sols = []; st = None
        while len(sols) < max_solutions:
            st = solver.Solve(m)
            if st not in (cp_model.OPTIMAL, cp_model.FEASIBLE): break
            sol = [i for i in range(self.norb) if solver.Value(x[i])]
            sols.append(sol)
            m.AddBoolOr([x[i].Not() if i in sol else x[i] for i in range(self.norb)])
        return sols, solver.StatusName(st)

    def elements(self, sol):
        return np.where(np.isin(self.orb, sol))[0]

    def cayley(self, sol):
        D = self.elements(sol)
        v = self.v; shape = np.array(self.shape)
        G = self.grid
        A = np.zeros((v, v), dtype=np.uint8)
        Dg = G[D]
        for x in range(v):
            nb = np.ravel_multi_index(tuple(((G[x] + Dg) % shape).T), self.shape)
            A[x, nb] = 1
        return A
