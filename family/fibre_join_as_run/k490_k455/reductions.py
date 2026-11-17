"""Centralizer equivalences and exact local specifications for the two rows."""
from audit import *


def polynomial_matrices(matrix,p):
    powers=[np.eye(4,dtype=np.int16)]
    for _ in range(3):powers.append(powers[-1]@matrix%p)
    return [sum((a*M for a,M in zip(co,powers)),np.zeros((4,4),dtype=np.int16))%p
            for co in itertools.product(range(p),repeat=4) if any(co)]


def main():
    D=data();P=D['P'];lookup={tuple(p):i for i,p in enumerate(P)}
    certificate=json.loads((SOURCE/'result/construction.json').read_text())
    fibre_cycles=D['fine_cycles'];where={j:i for i,cycle in enumerate(fibre_cycles) for j in cycle}
    actions=[]
    for matrix in polynomial_matrices(np.array(certificate['ternary_matrix'],dtype=np.int16),3):
        perm=[]
        for p in P:
            v=matrix@p%3;v=v*next(int(x) for x in v if x)%3;perm.append(lookup[tuple(v)])
        assert len(set(perm))==40
        action=[where[perm[cycle[0]]] for cycle in fibre_cycles]
        assert len(set(action))==8
        actions.append(action)
    actions=sorted(set(map(tuple,actions)));assert len(actions)==8
    def image(mask,action):return sum(1<<action[i] for i in range(8) if mask>>i&1)
    audit=json.loads((HERE/'audit.json').read_text());targets=[]
    current_mask=sum(1<<i for i,cycle in enumerate(fibre_cycles) if D['z'][0,cycle[0]+1])
    binary_actions=[]
    for matrix in polynomial_matrices(np.array(certificate['binary_matrix'],dtype=np.int16),2):
        perm=[]
        for b in range(16):
            digits=matrix@np.array([(b>>(3-j))&1 for j in range(4)],dtype=np.int16)%2
            perm.append(sum(int(x)<<(3-j) for j,x in enumerate(digits)))
        assert len(set(perm))==16
        binary_actions.append(perm)
    for row in audit['open_row_screens']:
        if row['parameters'][1] not in [455,490]:continue
        masks={f['fibre_mask'] for f in row['fibres']};todo=set(masks);classes=[]
        while todo:
            seed=min(todo);orbit={image(seed,a) for a in actions}
            assert orbit<=todo;todo-=orbit;classes.append(sorted(orbit))
        representatives=[(current_mask^255 if current_mask^255 in orbit else min(orbit)) for orbit in classes]
        quotients=row['fibres'][0]['quotients']
        qkeys={tuple(q['cross_sizes']) for q in quotients}
        first=quotients[0]['cross_sizes']
        orbit={tuple(first[b] for b in perm) for perm in binary_actions}
        assert orbit==qkeys
        targets.append(dict(parameters=row['parameters'],fibre_centralizer_classes=classes,
            representative_fibre_masks=representatives,quotient_centralizer_classes=1,
            representative_cross_sizes=first,essential_fibre_quotient_cases=len(classes),
            cross_seed_sizes=[first[cycle[0]] for cycle in D['binary_cycles']]))
    out=dict(current_degree525_fibre_mask=current_mask,ternary_projective_centralizer_actions=actions,
        binary_centralizer_order=len(binary_actions),targets=targets,
        local_formulas=dict(size27='T={0} union (+/- P), |P|=13; c_q=(|P intersect q-perp|-4)/3',
            size36='T=(+/- P), |P|=18; c_q=(|P intersect q-perp|-6)/3',
            rank4='For each q with rank 4, the 15 coefficients c are Walsh sums of a 4-subset with XOR zero (140 patterns)',
            rank8='For each q with rank 8, the 15 coefficients c are Walsh sums of an 8-subset with XOR zero (870 patterns)'),
        symmetry_scope='Only the stated C5 and inversion invariant abelian Cayley ansatz; these equivalences do not classify all SRGs')
    (HERE/'reductions.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()
