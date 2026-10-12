"""Re-run the UNKNOWN (job, parameter set) cases of the field-product searches with more
CP-SAT workers and a longer time limit.

python rerun_fields.py LOG NPROC WORKERS TLIM
Cases are taken from the logs listed in QUEUE, in that order (within a log: fewest orbits first).
Resumable: cases already in LOG are skipped. Hits are verified and written to hits_abelian/.
"""
import sys, os, ast, time, math
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"): os.environ[_v] = "1"
from multiprocessing import Pool
import numpy as np
from abelian_pds import AbelianPDS
from run_fields import field_mats, blockdiag, matpow_mod
from fusion import is_srg, graph6

QUEUE = [("logs/fields_2x2_3x4.log", "2^2,3^4"),
         ("logs/fields_2x2_17x2.log", "2^2,17^2"),
         ("logs/fields_5x2_7x2.log", "5^2,7^2"),
         ("logs/fields_2x4_5x2.log", "2^4,5^2"),
         ("logs/fields_2x4_7x2.log", "2^4,7^2"),
         ("logs/fields_2x6_3x2.log", "2^6,3^2"),
         ("logs/fields_2x4_3x4.log", "2^4,3^4"),
         ("logs/fieldsdeep_2x4_3x4.log", "2^4,3^4"),
         ("logs/fields_2x2_3x2_5x2.log", "2^2,3^2,5^2")]

_MATS = {}
def mats(spec):
    if spec not in _MATS:
        fields = [tuple(map(int, s.split("^"))) for s in spec.split(",")]
        Cs, Frs, shape = [], [], []
        for p, e in fields:
            C, F = field_mats(p, e); Cs.append(C); Frs.append(F); shape += [p] * e
        _MATS[spec] = (fields, Cs, Frs, tuple(shape))
    return _MATS[spec]

def work(case):
    spec, jobrepr, prm, W, T = case
    fields, Cs, Frs, shape = mats(spec)
    gs, g = ast.literal_eval(jobrepr)
    auts = [blockdiag([matpow_mod(C, a, p) for C, a, (p, e) in zip(Cs, gen, fields)]) for gen in gs]
    if g is not None:
        auts.append(blockdiag([matpow_mod(F, a, p) for F, a, (p, e) in zip(Frs, g, fields)]))
    P = AbelianPDS(shape, extra_auts=auts)
    k, la, mu = ast.literal_eval(prm)
    t = time.time()
    sols, st = P.search(k, la, mu, time_limit=T, workers=W, solver_params={"linearization_level": 2})
    v = int(np.prod(shape))
    if sols:
        A = P.cayley(sols[0]); ok = is_srg(A, k, la, mu)
        st = "FOUND" if ok else "FOUND-BAD"
        if ok:
            os.makedirs("hits_abelian", exist_ok=True)
            with open(f"hits_abelian/srg{v}_{k}_{la}_{mu}.g6", "a") as f: f.write(graph6(A) + "\n")
            with open(f"hits_abelian/srg{v}_{k}_{la}_{mu}.pds", "a") as f:
                f.write(f"{fields} job={jobrepr}\t{P.elements(sols[0]).tolist()}\n")
    return f"{spec}\t{jobrepr}\t{prm}\t{st}\torb={P.norb}\t{time.time()-t:.1f}s"

if __name__ == "__main__":
    LOG, NPROC, W, T = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), float(sys.argv[4])
    done = set()
    if os.path.exists(LOG):
        for l in open(LOG):
            p = l.rstrip("\n").split("\t")
            if len(p) >= 3: done.add((p[0], p[1], p[2]))
    cases = []
    for logf, spec in QUEUE:
        rows = []
        for l in open(logf):
            p = l.rstrip("\n").split("\t")
            if len(p) >= 4 and p[2] == "UNKNOWN":
                orb = int(p[3].split("=")[1])
                rows.append((orb, p[0], p[1]))
        rows.sort()
        for orb, jobrepr, prm in rows:
            if (spec, jobrepr, prm) not in done:
                cases.append((spec, jobrepr, prm, W, T))
    print(f"{len(cases)} cases to run", flush=True)
    with Pool(NPROC) as pool, open(LOG, "a") as f:
        for line in pool.imap(work, cases):
            f.write(line + "\n"); f.flush()
