"""Exact quotient/dual-parity obstruction; standard library only."""
import itertools
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
H = [[(-1)**((a & b).bit_count() % 2) for b in range(16)] for a in range(16)]


def cross(m, S):
    sizes = []
    for a in range(1, 16):
        n = 9*(m+sum(H[a][b] for b in S))
        if n % 4 or not 0 <= n//4 <= 81:
            return None
        sizes.append(n//4)
    return sizes


def gate(m, S):
    q = cross(m, S)
    if q is None:
        return False
    U = {a for a in range(1, 16) if q[a-1] % 2}
    dq = cross(m, U)
    return dq is not None and {a for a in range(1, 16) if dq[a-1] % 2} == S


def main():
    classes = {}
    for j in (4, 12):
        initial = final = 0
        for tup in itertools.combinations(range(1, 16), j):
            S = set(tup)
            if cross(12, S) is None:
                continue
            initial += 1
            final += gate(12, S)
        classes[str(j)] = dict(initial=initial, survives_dual_parity=final)
    assert classes == {'4': {'initial': 105, 'survives_dual_parity': 0},
                       '12': {'initial': 35, 'survives_dual_parity': 35}}
    controls = {}
    for degree, m, fn in [
            (525, 15, ROOT/'work/construct525_20261116/result/connections_Z6_4.txt'),
            (455, 13, ROOT/'work/transfer525_20261116/results/455/connections_Z6_4.txt'),
            (490, 14, ROOT/'work/transfer525_20261116/results/490/connections_Z6_4.txt')]:
        D = [tuple(map(int, line.split())) for line in fn.open()]
        size = [0]*16
        U = set()
        for v in D:
            a = sum((x % 2) << (3-i) for i, x in enumerate(v))
            size[a] += 1
            if all(x % 3 == 0 for x in v):
                U.add(a)
        ev = [sum(H[a][b]*size[b] for b in range(16)) for a in range(16)]
        S = {a for a in range(1, 16) if ev[a] == 36-m}
        assert gate(m, S)
        assert U == {a for a in range(1, 16) if cross(m, S)[a-1] % 2}
        controls[str(degree)] = True
    out = dict(parameters=[1296,420,144,132], quotient_classes=classes,
               known_graph_controls=controls, degree24_template_excluded=True,
               arbitrary_srg_parameters_excluded=False,
               scope='Inverse-closed Cayley graph on F2^4 x F3^4')
    (HERE/'dual_parity_gate.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
