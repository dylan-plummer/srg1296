"""Verify the square of the explicit order-four automorphism on class 24."""
from pathlib import Path
import sys
from verify_bundle import decode,check_translation

def digits(y):
    return tuple(y//3**i%3 for i in range(4))
def encode(vector):
    return sum(int(x)%3*3**i for i,x in enumerate(vector))
def linear(y):
    columns=tuple(map(digits,(47,34,76,14)))
    return encode(sum(t*c[i] for t,c in zip(digits(y),columns)) for i in range(4))
def check(raw):
    rows=decode(raw)
    perm=[81*((v//81)^1)+linear(v%81) for v in range(1296)]
    square=[perm[x] for x in perm]
    assert sorted(perm)==list(range(1296))
    check_translation(rows,perm)
    assert all(square[v]==81*(v//81)+encode(-t for t in digits(v%81)) for v in range(1296))
    assert sum(v==x for v,x in enumerate(square))==16
    print('Verified order-four automorphism; its square is ternary negation and fixes 16 vertices.')
if __name__=='__main__':
    path=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent/'classes/class_24.g6'
    check(path.read_bytes())
