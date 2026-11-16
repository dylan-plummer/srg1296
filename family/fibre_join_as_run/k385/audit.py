"""Exact C3 quotient, dual-parity, moment and invariant-fibre transfer gates.

Only the specified elementary-abelian Cayley family is screened. A local
survivor is not a constructed graph or a completed compatibility problem.
"""
import collections
import itertools
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent


def cycles(perm):
    unseen = set(range(len(perm)))
    out = []
    while unseen:
        i = min(unseen)
        cyc = [i]
        while perm[cyc[-1]] != i:
            cyc.append(perm[cyc[-1]])
        unseen -= set(cyc)
        out.append(cyc)
    return out


def moment_solutions(d, allowed):
    n = d//2
    target2 = (81*d-d*d)//2
    out = []
    def rec(i, left, first, second, counts):
        if i == len(allowed):
            if left == 0 and first == -n and second == target2:
                out.append(counts)
            return
        e = allowed[i]
        for a in range(left+1):
            if second+a*e*e > target2:
                continue
            rec(i+1,left-a,first+a*e,second+a*e*e,counts+[a])
    rec(0,40,0,0,[])
    return out


def main():
    V = list(itertools.product(range(2),repeat=4))
    ix = {v:i for i,v in enumerate(V)}
    p = [ix[(v[1],(v[0]+v[1])%2,v[3],(v[2]+v[3])%2)] for v in V]
    cc = [c for c in cycles(p) if c != [0]]
    H = np.array([[(-1)**((a&b).bit_count()%2) for b in range(16)] for a in range(16)],dtype=np.int16)
    P = np.array([v for v in itertools.product(range(3),repeat=4)
                  if any(v) and next(x for x in v if x)==1],dtype=np.int16)
    pd = {tuple(v):i for i,v in enumerate(P)}
    def norm(v):
        return tuple(v*next(int(x) for x in v if x)%3)
    planes = np.array([sum(int(a)<<j for j,a in enumerate(row)) for row in (P@P.T%3==0)],dtype=np.uint64)
    fine = {}
    for name,G in {
        '22': np.array([[1,0,1,0],[0,1,0,1],[0,0,1,0],[0,0,0,1]],dtype=np.int16),
        '31': np.array([[1,1,0,0],[0,1,1,0],[0,0,1,0],[0,0,0,1]],dtype=np.int16),
    }.items():
        perm = [pd[norm(G@v%3)] for v in P]
        dual = [pd[norm(G.T@v%3)] for v in P]
        po = np.array([sum(1<<j for j in c) for c in cycles(perm)],dtype=np.uint64)
        masks = np.zeros(1<<len(po),dtype=np.uint64)
        for i,m in enumerate(po):
            masks[(np.arange(len(masks))>>i&1).astype(bool)] |= m
        sizes = np.bitwise_count(masks).astype(np.int16)
        E = 3*np.bitwise_count(masks[:,None]&planes).astype(np.int16)-sizes[:,None]
        fixed = [i for i,j in enumerate(dual) if i == j]
        fine[name] = dict(masks=masks,sizes=sizes,E=E,fixed=fixed)
    rows = json.loads((ROOT/'work/fibre_revisit_20261116/ma_gate.json').read_text())['records']
    records = []
    for row in rows:
        v,k,lam,mu = row['parameters']
        if v != 1296 or row['status'] != 'SURVIVES_MA_GATES':
            continue
        disc = (lam-mu)**2+4*(k-mu)
        if disc != 36**2:
            continue
        r,s = (lam-mu+36)//2,(lam-mu-36)//2
        f = (-k-1295*s)//36
        rd,sd = (1296-k+s)//36,-(k-s)//36
        groups = {}
        trace_count = dual_count = 0
        for selector in range(32):
            S = {a for i,c in enumerate(cc) if selector>>i&1 for a in c}
            eigen = np.array([k]+[r if a in S else s for a in range(1,16)],dtype=np.int16)
            raw = H@eigen
            if np.any(raw%16):
                continue
            q = raw//16
            d = int(q[0])
            if d%2 or not 0<=d<=80 or np.any(q<0) or np.any(q>81):
                continue
            trace_count += 1
            U = {a for a in range(1,16) if q[a]%2}
            ed = np.array([f]+[rd if a in U else sd for a in range(1,16)],dtype=np.int16)
            qd = H@ed
            if np.any(qd%16):
                continue
            qd //= 16
            if qd[0]%2 or not 0<=qd[0]<=80 or np.any(qd<0) or np.any(qd>81):
                continue
            if {a for a in range(1,16) if qd[a]%2} != S:
                continue
            dual_count += 1
            key = (d,tuple(sorted(collections.Counter(map(int,q[1:])).items())))
            if key in groups:
                groups[key]['labelled_quotients'] += 1
                continue
            n = d//2
            allowed = sorted({3*h-n for h in range(min(n,13)+1)
                              if 0<=16*(3*h-n-s)<=16*36 and 16*(3*h-n-s)%36==0})
            moments = moment_solutions(d,allowed)
            tests = {}
            for name,loc in fine.items():
                keep = (loc['sizes']==n)&np.all(np.isin(loc['E'],allowed),axis=1)
                raw_count = int(keep.sum())
                ranks = (16*(loc['E'][:,loc['fixed']]-s))//36
                keep &= np.all(np.isin(ranks,[0,4,12,16]),axis=1)
                hist = collections.Counter(tuple(collections.Counter(map(int,e))[a] for a in allowed)
                                           for e in loc['E'][keep])
                tests[name] = dict(raw_invariant_fibres=raw_count,fixed_character_compatible=int(keep.sum()),
                                   projective_spectrum_histogram=[dict(counts=list(c),fibres=count) for c,count in sorted(hist.items())])
            groups[key] = dict(zero_degree=d,j=len(S),cross_degree_counts=dict(key[1]),
                               seed_sizes=[int(q[c[0]]) for c in cc],quotient_sizes=list(map(int,q)),
                               labelled_quotients=1,allowed_fibre_eigenvalues=allowed,moment_solutions=moments,
                               fine_actions=tests)
        records.append(dict(parameters=[v,k,lam,mu],r=r,s=s,trace_quotients=trace_count,
                            dual_parity_quotients=dual_count,cases=list(groups.values()),
                            locally_constructed=k in (385,420,455,490,525)))
    control = next(rec for rec in records if rec['parameters'][1]==420)
    assert len(control['cases'])==1 and control['cases'][0]['zero_degree']==42
    assert control['cases'][0]['fine_actions']['22']['fixed_character_compatible']==108
    out = dict(scope='Fixed-point-free coarse C3 on F2^4; fine Jordan 2+2 or 3+1 on F3^4',
               local_preflight_only=True,records=records)
    (HERE/'audit.json').write_text(json.dumps(out,indent=2)+'\n')
    for rec in records:
        print(rec['parameters'], 'dual quotients',rec['dual_parity_quotients'],
              [(c['zero_degree'],c['cross_degree_counts'],
                {name:loc['fixed_character_compatible'] for name,loc in c['fine_actions'].items()}) for c in rec['cases']])


if __name__=='__main__':
    main()
