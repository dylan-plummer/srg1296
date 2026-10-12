"""Search for strongly regular graphs that are unions of symmetrised orbitals of a
permutation group acting (possibly intransitively) on {0..v-1}.

Core API
  cc = Config(gens, v)               # gens: list of permutation arrays (images of 0..v-1)
  sols = cc.solve(k, lam, mu, ...)   # list of solutions (lists of relation indices)
  A = cc.adjacency(sol)              # 0/1 numpy matrix
  is_srg(A, k, lam, mu)
"""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from ortools.sat.python import cp_model


class Config:
    def __init__(self, gens, v):
        self.v = v
        gens = [np.asarray(g, dtype=np.int64) for g in gens]
        n = v * v
        idx = np.arange(n, dtype=np.int64)
        x, y = idx // v, idx % v
        rows, cols = [], []
        for g in gens:
            rows.append(idx)
            cols.append(g[x] * v + g[y])
        rows = np.concatenate(rows); cols = np.concatenate(cols)
        G = coo_matrix((np.ones(len(rows), dtype=np.int8), (rows, cols)), shape=(n, n))
        nc, lab = connected_components(G, directed=True, connection="weak")
        O = lab.reshape(v, v)
        # fibers = orbits on points (diagonal orbitals)
        dlab = np.unique(O[np.arange(v), np.arange(v)])
        fiber_of_orb = {int(d): i for i, d in enumerate(dlab)}
        self.fiber = np.array([fiber_of_orb[int(O[i, i])] for i in range(v)])
        self.nfib = len(dlab)
        # symmetrise: relation id = min(orbital, transpose orbital) over off-diagonal orbitals
        OT = O.T
        S = np.minimum(O, OT)
        diagset = set(int(d) for d in dlab)
        offd = sorted(set(np.unique(S).tolist()) - diagset)
        remap = -np.ones(nc, dtype=np.int64)
        for i, s in enumerate(offd):
            remap[s] = i
        R = remap[S]                      # -1 on diagonal
        self.R = R
        self.r = r = len(offd)
        # representative pair of each relation, fiber pair, valencies
        flat = R.ravel()
        first = np.full(r, -1, dtype=np.int64)
        order = np.argsort(flat, kind="stable")
        fl = flat[order]
        starts = np.searchsorted(fl, np.arange(r))
        first = order[starts]
        self.rep = [(int(f // v), int(f % v)) for f in first]
        self.fibrep = [int(np.where(self.fiber == f)[0][0]) for f in range(self.nfib)]
        self.fibsize = [int((self.fiber == f).sum()) for f in range(self.nfib)]
        # deg[f][i] = number of neighbours via relation i of a vertex in fiber f
        self.deg = np.zeros((self.nfib, r), dtype=np.int64)
        for f, a in enumerate(self.fibrep):
            row = R[a]
            row = row[row >= 0]
            self.deg[f] += np.bincount(row, minlength=r)
        # structure constants at representatives: c[l] = dict (i,j)->count
        self.coef = []
        for l, (a, b) in enumerate(self.rep):
            ra = R[a]; cb = R[:, b]
            m = (ra >= 0) & (cb >= 0)
            key = ra[m] * r + cb[m]
            u, cnt = np.unique(key, return_counts=True)
            d = {}
            for kk, c in zip(u.tolist(), cnt.tolist()):
                i, j = divmod(kk, r)
                if i > j: i, j = j, i
                d[(i, j)] = d.get((i, j), 0) + c
            self.coef.append(d)
        self.size = np.array([int((R == i).sum()) for i in range(r)])  # number of ordered pairs
        # ordered orbitals relabelled 0..N-1 (for the spectral decomposition)
        u, O2 = np.unique(O, return_inverse=True)
        self.O = O2.reshape(v, v)
        self.N = len(u)
        self.blocks = None

    # ------------------------------------------------------------------ spectral data
    def spectral(self, seed=1):
        """Decompose the centraliser algebra W (span of orbital matrices) into its simple
        blocks. For every block j returns (chi_j(1), m_j, trivial?, T_j) where T_j[l] =
        tr(A_l e_j) for the symmetrised relation l, e_j the central idempotent."""
        if self.blocks is not None:
            return self.blocks
        v, N, O = self.v, self.N, self.O
        cnt = np.bincount(O.ravel(), minlength=N)
        # structure constants p[i,j,m] = (A_i A_j)[a_m, b_m]
        first = np.zeros(N, dtype=np.int64)
        fl = O.ravel(); order = np.argsort(fl, kind="stable")
        first = order[np.searchsorted(fl[order], np.arange(N))]
        P = np.zeros((N, N, N))
        for mm in range(N):
            a, b = divmod(int(first[mm]), v)
            key = O[a] * N + O[:, b]
            P[:, :, mm] = np.bincount(key, minlength=N * N).reshape(N, N)
        # centre: sum_i z_i (p[i,l,m] - p[l,i,m]) = 0  for all l, m
        C = (P - P.transpose(1, 0, 2)).transpose(1, 2, 0).reshape(N * N, N)
        _, sv, vt = np.linalg.svd(C, full_matrices=False)
        tol = 1e-8 * max(1.0, sv[0])
        Z = vt[sv < tol] if (sv < tol).any() else np.zeros((0, N))
        # nullity may be under-reported by full_matrices=False when N*N >= N (fine)
        rng = np.random.default_rng(seed)
        z = rng.standard_normal(len(Z)) @ Z
        Mz = z[O]
        Mz = Mz + Mz.T
        w, U = np.linalg.eigh(Mz)
        # cluster eigenvalues
        groups = []
        start = 0
        scale = max(1.0, np.abs(w).max())
        for i in range(1, v + 1):
            if i == v or w[i] - w[i - 1] > 1e-7 * scale:
                groups.append((start, i)); start = i
        ones_f = np.zeros((v, self.nfib))
        ones_f[np.arange(v), self.fiber] = 1.0
        # transpose pairing of ordered orbitals and map ordered -> symmetrised relation
        tr_of = self.O.T
        symrel = np.full(N, -1)
        symrel_flat = self.R.ravel()
        symrel[self.O.ravel()] = symrel_flat
        blocks = []
        for (s0, s1) in groups:
            Uj = U[:, s0:s1]
            d = s1 - s0
            E = Uj @ Uj.T
            t = np.bincount(O.ravel(), weights=E.ravel(), minlength=N)   # tr(A_o e_j) (E symmetric)
            # m_j^2 = dim W e_j = rank of the matrix of e_j * A_i in orbital coordinates
            c = t / cnt
            L = np.einsum("o,oim->im", c, P)       # e_j A_i in coords (rows i)
            m2 = np.linalg.matrix_rank(L, tol=1e-6 * max(1.0, np.abs(L).max()))
            mj = int(round(np.sqrt(m2)))
            chi = d / mj
            triv = bool(np.abs(Uj.T @ ones_f).max() > 1e-6)
            T = np.zeros(self.r)
            for o in range(N):
                l = symrel[o]
                if l >= 0:
                    T[l] += t[o]
            Qo = np.einsum("abm,m->ab", P, t)          # tr(A_a A_b e_j), ordered orbitals
            Q = np.zeros((self.r, self.r))
            for a in range(N):
                la = symrel[a]
                if la < 0: continue
                for b_ in range(N):
                    lb = symrel[b_]
                    if lb < 0: continue
                    Q[la, lb] += Qo[a, b_]
            blocks.append((chi, mj, triv, T, d, Q))
        self.blocks = blocks
        return blocks

    def spectral_constraints(self):
        """Integral linear constraints: list of (coef vector (int), chi, m, triv).
        Irrational blocks are merged per (chi, m) (Galois-stable unions)."""
        out = []
        pend = {}
        for chi, mj, triv, T, d, Q in self.spectral():
            Tc = T / chi
            Qc = Q / chi
            if (np.abs(Tc - np.round(Tc)).max() < 1e-5 and np.abs(Qc - np.round(Qc)).max() < 1e-5
                    and abs(chi - round(chi)) < 1e-6):
                out.append((np.round(Tc).astype(np.int64), int(round(chi)), mj, triv, 1,
                            np.round(Qc).astype(np.int64)))
            else:
                key = (round(chi, 6), mj, triv)
                if key in pend:
                    p = pend[key]; p[0] += T; p[1] += 1; p[2] += Q
                else:
                    pend[key] = [T.copy(), 1, Q.copy()]
        for (chi, mj, triv), (T, nb, Q) in pend.items():
            Tc = T / chi; Qc = Q / chi
            if np.abs(Tc - np.round(Tc)).max() < 1e-5:
                Qi = np.round(Qc).astype(np.int64) if np.abs(Qc - np.round(Qc)).max() < 1e-5 else None
                out.append((np.round(Tc).astype(np.int64), int(round(chi)), mj, triv, nb, Qi))
        return out

    def solve(self, k, lam, mu, max_solutions=1, time_limit=60.0, workers=8, forbid=None,
              spectral=True):
        r = self.r
        m = cp_model.CpModel()
        x = [m.NewBoolVar(f"x{i}") for i in range(r)]
        for f in range(self.nfib):
            m.Add(sum(int(self.deg[f, i]) * x[i] for i in range(r) if self.deg[f, i]) == k)
        if spectral:
            import math
            D = (lam - mu) ** 2 + 4 * (k - mu)
            sq = math.isqrt(D)
            if sq * sq == D:
                rr = ((lam - mu) + sq) // 2; ss = ((lam - mu) - sq) // 2
                f_mult = []
                spec2 = []
                for coef, chi, mj, triv, nb, Qi in self.spectral_constraints():
                    M = mj * nb - (1 if triv else 0)       # number of eigenvalue slots r/s
                    b = m.NewIntVar(0, M, "")
                    rhs = ss * M + (k if triv else 0)
                    # sum coef*x = rhs + (rr-ss) b
                    m.Add(sum(int(c) * x[i] for i, c in enumerate(coef) if c) == rhs + (rr - ss) * b)
                    f_mult.append((chi, b))
                    if Qi is not None:
                        spec2.append((Qi, b, ss * ss * M + (k * k if triv else 0), rr * rr - ss * ss))
                fmul = (-k - ss * (self.v - 1)) // (rr - ss)
                m.Add(sum(chi * b for chi, b in f_mult) == fmul)
        prods = {}
        def prod(i, j):
            if i == j: return x[i]
            if (i, j) not in prods:
                y = m.NewBoolVar(f"y{i}_{j}")
                m.AddBoolAnd([x[i], x[j]]).OnlyEnforceIf(y)
                m.AddBoolOr([x[i].Not(), x[j].Not()]).OnlyEnforceIf(y.Not())
                prods[(i, j)] = y
            return prods[(i, j)]
        if spectral and 'spec2' in dir():
            for Qi, b, rhs2, dlt in spec2:
                terms = []
                for i in range(r):
                    for j in range(i, r):
                        c = int(Qi[i, j] + (Qi[j, i] if j != i else 0))
                        if c: terms.append(c * prod(i, j))
                m.Add(sum(terms) == rhs2 + dlt * b)
        for l in range(r):
            m.Add(sum(c * prod(i, j) for (i, j), c in self.coef[l].items()) == mu + (lam - mu) * x[l])
        if forbid:
            for sol in forbid:
                s = set(sol)
                m.AddBoolOr([x[i].Not() if i in s else x[i] for i in range(r)])
        sols = []
        solver = cp_model.CpSolver()
        solver.parameters.max_time_in_seconds = time_limit
        solver.parameters.num_workers = workers
        status = None
        while len(sols) < max_solutions:
            status = solver.Solve(m)
            if status not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
                break
            sol = [i for i in range(r) if solver.Value(x[i])]
            sols.append(sol)
            m.AddBoolOr([x[i].Not() if i in sol else x[i] for i in range(r)])
        st = solver.StatusName(status) if status is not None else "NA"
        return sols, st

    def adjacency(self, sol):
        mask = np.zeros(self.r + 1, dtype=np.uint8)
        mask[np.array(sol)] = 1
        return mask[self.R]  # R==-1 picks mask[-1]=0


def is_srg(A, k, lam, mu):
    A = A.astype(np.int64)
    v = A.shape[0]
    if (A != A.T).any() or A.diagonal().any(): return False
    if (A.sum(1) != k).any(): return False
    A2 = A @ A
    J = np.ones((v, v), dtype=np.int64); I = np.eye(v, dtype=np.int64)
    return bool((A2 == k * I + lam * A + mu * (J - I - A)).all())


def graph6(A):
    v = A.shape[0]
    bits = []
    for j in range(1, v):
        for i in range(j):
            bits.append(int(A[i, j]))
    bits += [0] * (-len(bits) % 6)
    if v < 63: head = [63 + v]
    elif v < 258048: head = [126, 63 + (v >> 12 & 63), 63 + (v >> 6 & 63), 63 + (v & 63)]
    out = head + [63 + int("".join(map(str, bits[t:t + 6])), 2) for t in range(0, len(bits), 6)]
    return bytes(out).decode()
