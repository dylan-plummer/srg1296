"""Verify the switch from the archived Cayley graph6 files, using integer sums."""
from functools import cache
from itertools import product
from pathlib import Path
from tempfile import TemporaryDirectory
from zipfile import ZipFile
import json
from verify_graph6 import read_graph

HERE=Path(__file__).resolve().parent
Y=list(product(range(3),repeat=4))
DIFF=[[Y.index(tuple((z-v)%3 for v,z in zip(a,b))) for b in Y] for a in Y]
def chi(a,x):return (-1)**((a&x).bit_count()%2)

@cache
def kernels(source,freq):
    # The archived source files number x first, and the first ternary
    # coefficient varies fastest there. The new files order ternary tuples
    # lexicographically and then place all 16 binary points consecutively.
    path=HERE.parent/'graphs'/Path(source).with_suffix('.g6')
    rows,_=read_graph(path)
    out=[]
    for y in Y:
        old=sum(v*3**i for i,v in enumerate(y))
        b=[(rows[0]>>(81*x+old))&1 for x in range(16)]
        spec={a:sum(t*chi(a,x) for x,t in enumerate(b)) for a in freq}
        nums=[16*t-2*sum(chi(a,x)*v for a,v in spec.items()) for x,t in enumerate(b)]
        bp=[t//16 for t in nums] if all(t in (0,16) for t in nums) else None
        def translate(v):
            return [sum(v[x^z]<<z for z in range(16)) for x in range(16)]
        out.append((translate(b),translate(bp) if bp is not None else None))
    return out

def main():
    counts={}
    with TemporaryDirectory() as tmp, ZipFile(HERE/'non_cayley_68.zip') as z:
        scratch=Path(tmp)/'test.g6'
        for r in json.loads(z.read('manifest.json')):
            table=kernels(r['source'],tuple(r['negative_frequencies']))
            if 'fibre_mask' in r:
                side=[int(r['fibre_mask'],16)>>j&1 for j in range(81)]
                assert all(b is not None for a,b in table[1:])
                changed=sum(a!=b for a,b in table[1:])
                assert changed==52
                counts[r['parameters'][1]]=changed
            else:
                side=[sum(a*b for a,b in zip(r['normal'],y))%3==0 for y in Y]
            scratch.write_bytes(z.read(r['file']))
            rows,_=read_graph(scratch)
            for i in range(81):
                blocks=[table[d][int(side[i]!=side[j])] for j,d in enumerate(DIFF[i])]
                assert all(b is not None for b in blocks)
                for x in range(16):
                    row=sum(b[x]<<(16*j) for j,b in enumerate(blocks))
                    assert row==rows[16*i+x],r['file']
    print(json.dumps(dict(archive_graphs=68,changed_nonzero_fibre_differences=counts)))
if __name__=='__main__':main()
