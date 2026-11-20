"""Search spectral sign changes on one side of a hyperplane of ternary fibres."""
from pathlib import Path
from itertools import product
from functools import cache
import json,sys
import numpy as np
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'family'))
from fields import read_connset
X=list(range(16));Y=list(product(range(3),repeat=4))
H=np.array([[1-2*((i&j).bit_count()%2) for j in X] for i in X],dtype=np.int16)
# Mod out all 16 character gauges by fixing their signs on the four standard
# basis vectors. These gauges are simply permutations of the binary points.
free=[x for x in X if x not in (0,1,2,4,8)]
signs=np.ones((1<<len(free),16),dtype=np.int16)
for i,j in enumerate(free):signs[:,j]=1-2*((np.arange(len(signs))>>i)&1)
def bilinear(u,v):return sum(a*b for a,b in zip(u,v))%3
normal=[u for u in Y[1:] if next(t for t in u if t)==1]
@cache
def eligible(s):
    b=np.array([(s>>i)&1 for i in X],dtype=np.int16);spec=b@H
    transformed=(signs*spec)@H # sixteen times the new binary vector
    return np.all((transformed==0)|(transformed==16),axis=1),np.any(transformed!=16*b,axis=1)
def data(path):
    param,D=read_connset(path);rows={y:0 for y in Y}
    for a,y in D:rows[y]|=1<<sum(t<<i for i,t in enumerate(a))
    return param,rows
def scan(path):
    param,rows=data(path);found=[];counts=[]
    for u in normal:
        used={rows[y] for y in Y if bilinear(u,y)!=0}
        valid=np.ones(len(signs),dtype=bool);changed=np.zeros(len(signs),dtype=bool)
        for s in used:
            v,c=eligible(s);valid&=v;changed|=c
        changes=np.flatnonzero(valid&changed)
        counts.append(len(changes))
        if len(changes):found.append([u,changes[:16].tolist()])
    return dict(file=path.name,parameters=param,counts=counts,examples=found)
if __name__=='__main__':
    for p in sorted((ROOT.parent/'connsets').glob('*.txt')):
        r=scan(p);print(json.dumps(r),flush=True)
