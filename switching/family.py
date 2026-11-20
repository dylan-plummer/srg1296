"""Verify binary switches on arbitrary ternary fibres and make a reproducible sample."""
import argparse
from hashlib import shake_256
import json
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from construct import build_cut, graph6, invariant
from scan import ROOT, Y, data, eligible, signs
from verify_graph6 import verify


GAUGES = {560:700, 592:933}


def run(degree, count):
    w = GAUGES[degree]
    path = next((ROOT.parent/'connsets').glob(f'srg1296_{degree}_*_class1.txt'))
    param, rows = data(path)
    changed = []
    for y in Y[1:]:
        valid, different = eligible(rows[y])
        assert valid[w]
        if different[w]:changed.append(y)
    # Starting with one vertex, traverse all changed fibre differences.
    reachable = {(0,0,0,0)}
    while True:
        expanded = reachable | {tuple((a+b)%3 for a,b in zip(r,d)) for r in reachable for d in changed}
        if expanded == reachable:break
        reachable = expanded
    assert len(reachable) == 81
    header = dict(source=path.name, parameters=param, negative_frequencies=np.flatnonzero(signs[w]<0).tolist(),
                  changed_nonzero_differences=len(changed), components=1,
                  distinct_labelled_graphs=str(2**80))
    (ROOT/f'family_{degree}.json').write_text(json.dumps(header,indent=2)+'\n')
    with TemporaryDirectory() as tmp:
        for i in range(count):
            label=f'SRG binary fibre switch {degree} sample {i}'
            mask = int.from_bytes(shake_256(label.encode()).digest(10),'little') << 1
            side = [bool(mask >> j & 1) for j in range(81)]
            _,a=build_cut(path,side,w)
            binary=graph6(a)
            outfile=Path(tmp)/'candidate.g6'
            outfile.write_bytes(binary)
            cert=verify(outfile)
            assert cert['parameters']==list(param)
            local=cert['four_cliques_through_vertex']
            total=sum(k*v for k,v in local.items())//4
            cert.pop('file')
            cert.pop('representative_vertices')
            cert.update(index=i, fibre_mask=hex(mask), total_four_cliques=total)
            if i==0:
                cert.update(invariant(a))
                (ROOT/f'arbitrary_{degree}.g6').write_bytes(binary)
            print(json.dumps(cert),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('degree',type=int,choices=GAUGES)
    parser.add_argument('count',type=int)
    args=parser.parse_args()
    run(args.degree,args.count)
