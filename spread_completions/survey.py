"""Find compatible coset spreads in the local Cayley examples on F2^4 x F3^4."""
from collections import Counter
from functools import cache
from itertools import combinations
import json
from pathlib import Path
import sys
import numpy as np
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT/'family'))
sys.path.insert(0, str(HERE.parent/'fibre_switch_20261120'))
from fields import read_connset
from verify_graph6 import read_graph

def coefficients(x, p):
    return tuple((x // p**i) % p for i in range(4))

def number(a, p):
    return sum(x*p**i for i,x in enumerate(a))

def addition(p):
    vectors = [coefficients(i,p) for i in range(p**4)]
    return np.array([[number([(a+b)%p for a,b in zip(x,y)],p) for y in vectors] for x in vectors])
ADD2, ADD3 = addition(2), addition(3)
def planes(add, p):
    subs = set()
    for i,j in combinations(range(1,len(add)),2):
        a,b=0,0
        points=set()
        for _ in range(p):
            for _ in range(p):
                points.add(int(add[a,b]))
                b=add[b,j]
            a=add[a,i]
        if len(points)==p*p:
            subs.add(tuple(sorted(points)))
    return sorted(subs)
PLANES2, PLANES3 = planes(ADD2,2), planes(ADD3,3)
assert (len(PLANES2),len(PLANES3))==(35,130)

def kernel(path):
    if path.suffix=='.g6':
        rows,_=read_graph(path)
        return np.array([rows[0]>>i & 1 for i in range(1296)],dtype=np.uint8).reshape(16,81)
    _,d=read_connset(path)
    out=np.zeros((16,81),dtype=np.uint8)
    for a,b in d:
        out[number(a,2),number(b,3)]=1
    return out

def subgroups(d, size):
    out=[]
    for a in PLANES2:
        piece=d[list(a)].sum(axis=0)
        for b in np.flatnonzero(piece[np.array(PLANES3)].sum(axis=1)==size):
            out.append((a,PLANES3[b]))
    return out

def meet(u,v):
    return len(set(u[0])&set(v[0]))*len(set(u[1])&set(v[1]))

def sources():
    return sorted((ROOT/'connsets').glob('*.txt'))+sorted((ROOT/'graphs/koval').glob('*.g6'))

if __name__=='__main__':
    for p in sources():
        d=kernel(p)
        if d.sum()%35:continue
        for size in (0,35):
            us=subgroups(d,size)
            if us:
                print(json.dumps(dict(source=p.name,within=size,subgroups=len(us),intersections=dict(Counter(meet(u,v) for u,v in combinations(us,2))))),flush=True)
