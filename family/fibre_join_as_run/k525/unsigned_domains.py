"""Project complete local code domains onto shared-frame unsigned patterns."""
from pathlib import Path
import json,itertools,collections,time
import numpy as np
HERE=Path(__file__).resolve().parent


def cycles(p,subset):
    todo=set(subset);out=[]
    while todo:
        a=min(todo);x=a;orbit=[]
        while x not in orbit:orbit.append(x);x=p[x]
        assert set(orbit)<=todo;todo-=set(orbit);out.append(orbit)
    return out


def load():
    D=json.loads((HERE/'source.json').read_text());P=np.array(D['points'],dtype=np.int16)
    N=P@P.T%3==0;plane=np.array([sum(int(x)<<j for j,x in enumerate(row)) for row in N],dtype=np.uint64)
    masks={27:[],36:[]}
    for line in (HERE/'local_sets.tsv').open():
        size,which,mask=line.split();masks[int(size)].append(int(mask,16))
    arrays={};counts={}
    for size,values in masks.items():
        z=np.array(values,dtype=np.uint64);C=np.empty((len(z),40),dtype=np.int8)
        for first in range(0,len(z),10000):
            hits=np.bitwise_count(z[first:first+10000,None]&plane[None,:]).astype(np.int16)
            assert np.all((hits-(4 if size==27 else 6))%3==0)
            C[first:first+10000]=(hits-(4 if size==27 else 6))//3
        np.savez_compressed(HERE/f'local{size}.npz',masks=z,coefficients=C)
        hist=collections.Counter(map(int,(C[:,np.array(D['diagonal_fourier'])==3].astype(np.int16)**2).sum(1)))
        arrays[size]=(z,C);counts[str(size)]=dict(sets=len(z),positive_square_histogram=dict(sorted(hist.items())))
    return D,arrays,counts


def main():
    started=time.monotonic();D,arrays,counts=load()
    action=next(g for g in D['group_actions'] if g['binary']!=list(range(16)))
    binary=action['binary'];dual=action['dual'];inv=np.argsort(dual)
    diffcycles=cycles(binary,range(1,16));freqcycles=cycles(dual,range(40));cap=set(D['binary_cap'])
    assert len(diffcycles)==3 and len(freqcycles)==8 and all(len(o)==5 for o in diffcycles+freqcycles)
    lines={tuple(sorted((a,b,a^b))) for a,b in itertools.combinations(range(1,16),2)}
    rank4=[];rank8=[]
    for line in lines:
        pattern=[int(d in line) for d in range(1,16)];rank4.append(pattern)
        U={0,*line}
        for a in range(1,16):
            if a not in U:rank8.append([int(d in {a^b for b in U}) for d in range(1,16)])
    rank8=[list(p) for p in sorted(set(map(tuple,rank8)))]+[[2*int(d==a) for d in range(1,16)] for a in range(1,16)]
    assert len(rank4)==35 and len(rank8)==120
    templates=[]
    for fo in freqcycles:
        p=fo[0]
        if p in D['matching_projective_channels']:
            d=D['binary_cap'][D['matching_projective_channels'].index(p)]
            patterns=[[2*int(a==d) for a in range(1,16)]]
        else:patterns=rank4 if D['diagonal_fourier'][p]==-6 else rank8
        triples=[]
        for pattern in patterns:triples.append([sum(pattern[d-1]*3**k for k,d in enumerate(dc)) for dc in diffcycles])
        templates.append(sorted(set(map(tuple,triples))))
    outputs=[]
    for b,dc in enumerate(diffcycles):
        size=27 if set(dc)==cap else 36;z,C=arrays[size];keep=np.ones(len(z),dtype=bool)
        codes=np.empty((len(z),8),dtype=np.int16)
        for f,fo in enumerate(freqcycles):
            p=fo[0];ids=[]
            for _ in range(5):ids.append(int(p));p=inv[p]
            codes[:,f]=np.abs(C[:,ids]).astype(np.int16)@np.array([3**i for i in range(5)],dtype=np.int16)
            allowed={t[b] for t in templates[f]};keep &= np.isin(codes[:,f],list(allowed))
        selected=np.flatnonzero(keep);unique,inverse=np.unique(codes[keep],axis=0,return_inverse=True)
        np.savez_compressed(HERE/f'bucket{b}.npz',codes=unique,local_indices=selected,option_indices=inverse)
        outputs.append(dict(index=b,differences=dc,cross_size=size,local_sets_surviving=len(selected),unsigned_options=len(unique)))
    out=dict(scope='Complete local domains with fixed five matching channels and shared Sylvester frame; extra matching channels permitted',
             local_counts=counts,buckets=outputs,frequency_cycles=freqcycles,difference_cycles=diffcycles,
             frequency_templates=[[list(t) for t in ts] for ts in templates],
             action=action,seconds=time.monotonic()-started)
    (HERE/'unsigned_domains.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ['action','frequency_templates']},indent=2),flush=True)


if __name__=='__main__':main()
