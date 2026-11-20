"""Independent exact-cover count and cross-checks for the completed extension census."""
import json
from hashlib import sha256
from collections import Counter
from ortools.sat.python import cp_model
from build import HERE, V, elements, kernel, partition, sources, subgroups
from spread_census import enumerate_spreads

def load_lines(filename):
    return [json.loads(s) for s in (HERE/filename).read_text().splitlines()]

def main():
    p=next(p for p in sources() if p.name=='srg1296_560_244_240_class10.txt')
    groups=subgroups(kernel(p),0)
    cells=[block for u in groups for block in partition(elements(u),V)]
    model=cp_model.CpModel()
    take=[model.new_bool_var(f'c{i}') for i in range(len(cells))]
    meet=[[] for _ in V]
    for i,block in enumerate(cells):
        for v in block:meet[v].append(i)
    for ids in meet:model.add_exactly_one([take[i] for i in ids])
    class Solutions(cp_model.CpSolverSolutionCallback):
        def __init__(self):
            super().__init__()
            self.covers=set()
        def on_solution_callback(self):
            s=sum(1<<i for i,v in enumerate(take) if self.value(v))
            assert s.bit_count()==36 and s not in self.covers
            self.covers.add(s)
    solver=cp_model.CpSolver();solver.parameters.enumerate_all_solutions=True
    solver.parameters.num_workers=1
    sol=Solutions()
    status=solver.solve(model,sol)
    assert status==cp_model.OPTIMAL and sol.covers==enumerate_spreads(cells)
    digest=sha256(b''.join(m.to_bytes(36,'little') for m in sorted(sol.covers))).hexdigest()
    census=json.loads((HERE/'census.json').read_text())
    refs=[r for r in load_lines('compare595.jsonl') if r['file']!='mixed595.g6']
    classes=load_lines('class_checks.jsonl')+load_lines('uniform6_checks.jsonl')
    fp=lambda r:tuple(map(tuple,r['K4_per_vertex']))
    assert len(refs)==11 and len(set(r['bliss_sha256'] for r in refs))==10
    assert len(classes)==25 and len(set(map(fp,classes)))==25
    assert not (set(map(fp,classes))&set(map(fp,refs)))
    lookup={r['file']:r for r in classes}
    for rec in census['records']:
        out=lookup[rec['file']]
        assert out['bliss_sha256']==rec['canonical_sha256']
        assert out['automorphisms']*rec['orbit']==census['host_aut']
    assert sum(r['orbit'] for r in census['records'])==len(sol.covers)==1124
    assert len({r['file'] for r in census['records']})==24
    assert len(load_lines('uniform6.jsonl'))==5
    print(json.dumps(dict(exact_covers=len(sol.covers),cover_digest=digest,cp_sat_status=solver.status_name(status),
        I10_graphs=len(census['records']),all_graphs=len(classes),shared_reference_classes=10,
        distinct_clique_profiles=len(set(map(fp,classes))),shared_profile_matches=0,
        automorphism_orders=dict(sorted(Counter(r['automorphisms'] for r in classes).items()))),indent=2))
if __name__=='__main__':main()
