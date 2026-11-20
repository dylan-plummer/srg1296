"""Rebuild every archived graph from two sparse Cayley seeds and the partition recipes."""
import argparse,json
from pathlib import Path
from zipfile import ZipFile

def ternary(x,y):
    return sum(((x//3**j+y//3**j)%3)*3**j for j in range(4))
T=[[ternary(x,y) for y in range(81)] for x in range(81)]
def add(x,y):
    return 81*((x//81)^(y//81))+T[x%81][y%81]
def cosets(group):
    out=[];left=set(range(1296))
    while left:
        v=min(left)
        block=sorted(add(v,x) for x in group)
        left.difference_update(block)
        out.append(block)
    return out
def plane_cosets(plane):
    return cosets([81*a+b for a in plane[0] for b in plane[1]])
def encode(rows):
    data=bytearray(b'~?SO')
    value=bits=0
    for j,row in enumerate(rows):
        for i in range(j):
            value=(value<<1)|((row>>i)&1)
            bits+=1
            if bits==6:
                data.append(63+value)
                value=bits=0
    assert not bits
    return bytes(data)+b'\n'
def build(z,out):
    census=json.loads(z.read('census.json'))
    blocks=[b for plane in census['free_groups'] for b in plane_cosets(plane)]
    cases={'classes/'+r['file']:[blocks[i] for i in r['blocks']] for r in census['records']}
    u=json.loads(z.read('uniform6.jsonl').splitlines()[0])
    cases[u['file']]=plane_cosets(u['subgroup'])
    for source in ('I.10','I.6'):
        seed=json.loads(z.read('seeds/'+source+'.json'))
        assert len(set(seed))==560 and 0 not in seed
        rows=[sum(1<<add(v,d) for d in seed) for v in range(1296)]
        for name,spread in cases.items():
            if (name.startswith('classes/'))!=(source=='I.10'):continue
            a=rows.copy()
            assert sorted(x for b in spread for x in b)==list(range(1296))
            for b in spread:
                mask=sum(1<<v for v in b)
                for v in b:
                    assert not mask&a[v]
                    a[v]|=mask^(1<<v)
            raw=encode(a)
            assert raw==z.read(name),name
            if out:
                p=out/name
                p.parent.mkdir(exist_ok=True,parents=True)
                if p.exists() and p.read_bytes()!=raw:raise FileExistsError(p)
                p.write_bytes(raw)
    print('Rebuilt all 25 graph6 files byte-for-byte from two degree-560 seeds.')
if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('archive')
    ap.add_argument('--out',type=Path)
    args=ap.parse_args()
    with ZipFile(args.archive) as z:build(z,args.out)
