"""Tier 3: re-run the cases still UNKNOWN in $SRC (default logs/rerun_fields.log; $SKIP = logs whose cases to skip), one at a time with many
CP-SAT workers and a long time limit.   python rerun_tier3.py LOG WORKERS TLIM
Resumable (cases already in LOG are skipped). Hits are verified and written to hits_abelian/."""
import sys, os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"): os.environ[_v] = "1"
from rerun_fields import work

V_OF = {"2^4,3^4": 1296, "2^6,3^2": 576, "2^4,5^2": 400, "2^4,7^2": 784, "2^2,3^4": 324, "2^2,3^2,5^2": 900}
ORDER = ["2^4,3^4", "2^6,3^2", "2^4,5^2", "2^4,7^2", "2^2,3^4", "2^2,3^2,5^2"]

if __name__ == "__main__":
    LOG, W, T = sys.argv[1], int(sys.argv[2]), float(sys.argv[3])
    done = set()
    for lf in [LOG] + [x for x in os.environ.get("SKIP", "").split(",") if x]:
        if os.path.exists(lf):
            done |= {tuple(l.split("\t")[:3]) for l in open(lf) if l.count("\t") >= 3}
    open_ = None
    if os.environ.get("OPEN_ONLY"):   # skip parameter sets no longer open in open_sets.tsv
        open_ = {tuple(map(int, l.split())) for l in open("open_sets.tsv") if l.strip()}
    rows = []
    for l in open(os.environ.get("SRC", "logs/rerun_fields.log")):
        p = l.rstrip("\n").split("\t")
        if len(p) >= 5 and p[3] == "UNKNOWN" and tuple(p[:3]) not in done:
            if open_ is not None and (V_OF[p[0]], *map(int, p[2].strip("()").split(","))) not in open_:
                continue
            rows.append((ORDER.index(p[0]) if p[0] in ORDER else 99, int(p[4].split("=")[1]), p[0], p[1], p[2]))
    rows.sort()
    print(f"{len(rows)} cases to run", flush=True)
    with open(LOG, "a") as f:
        for _, _, spec, job, prm in rows:
            f.write(work((spec, job, prm, W, T)) + "\n"); f.flush()
