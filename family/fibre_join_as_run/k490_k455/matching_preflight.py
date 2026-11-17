"""Exact degree-455 local preflight for the successful five-matching architecture."""
from audit import *


def main():
    started=time.monotonic();D=data();P=D['P'];planes=np.array([sum(int(a)<<j for j,a in enumerate(row)) for row in D['N']],dtype=np.uint64)
    masks=np.array([int(line.split()[1],16) for line in (HERE/'local_sets_0.tsv').open() if line.startswith('27\t')],dtype=np.uint64)
    assert len(masks)==353340
    hits=np.bitwise_count(masks[:,None]&planes).astype(np.int16)
    assert np.all((hits-4)%3==0)
    C=((hits-4)//3).astype(np.int8)
    zz=np.zeros(41,dtype=np.int16)
    for i,cycle in enumerate(D['fine_cycles']):
        if 124>>i&1:zz[np.array(cycle)+1]=1
    eigen=(D['character'][:41,:41]@zz)[1:]
    assert np.all(np.abs(C[:,eigen==-4])<=1) and np.all(np.abs(C[:,eigen==5])<=2)
    positive=eigen==5;hist=collections.Counter(map(int,(C[:,positive].astype(np.int16)**2).sum(1)))
    cert=json.loads((SOURCE/'result/construction.json').read_text());F=np.array(cert['ternary_matrix'],dtype=np.int16)
    lookup={tuple(p):i for i,p in enumerate(P)};perm=[]
    for p in P:
        v=F.T@p%3;v=v*next(int(x) for x in v if x)%3;perm.append(lookup[tuple(v)])
    freq=[cycle for cycle in cycles(perm) if eigen[cycle[0]]==5]
    assert len(freq)==3
    cases=[]
    for i,cycle in enumerate(freq):
        zero=np.flatnonzero(np.all(C[:,cycle]==0,axis=1))
        special=np.flatnonzero(np.all(C[:,cycle]==np.array([2,0,0,0,0]),axis=1))
        np.savez_compressed(HERE/f'matching455_{i}.npz',zero_masks=masks[zero],zero_coefficients=C[zero],
            cap_masks=masks[special],cap_coefficients=C[special],frequency_cycle=np.array(cycle))
        cases.append(dict(index=i,matching_frequency_cycle=cycle,ordinary_local_sets=len(zero),special_local_sets=len(special)))
    out=dict(parameters=[1296,455,166,156],fibre_mask=124,complete_local_size27_sets=len(masks),
        positive_energy_histogram=dict(sorted(hist.items())),required_three_seed_positive_energy=12,
        matching_cases=cases,seconds=time.monotonic()-started,
        scope='Five rank-eight matching channels imposed as an optional construction ansatz; not required in all Cayley graphs')
    (HERE/'matching_preflight.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))


if __name__=='__main__':main()
