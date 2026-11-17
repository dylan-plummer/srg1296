"""Parameter-changing unions of seventeen natural pieces of each witness."""
import csv
import itertools
import json
import subprocess
from pathlib import Path
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent.parent


def main():
    subprocess.run(['g++','-O3','-std=c++17',str(HERE/'fusions.cpp'),'-o','/tmp/srg420_fusions'],check=True)
    P=np.array([v for v in itertools.product(range(3),repeat=4) if any(v) and next(x for x in v if x)==1],dtype=np.int16)
    pd={tuple(p):j+1 for j,p in enumerate(P)}
    H=np.array([[(-1)**((a&b).bit_count()%2) for b in range(16)] for a in range(16)],dtype=np.int16)
    fine=np.ones((41,41),dtype=np.int16);fine[0,1:]=2;fine[1:,1:]=np.where(P@P.T%3==0,2,-1)
    C=np.kron(H,fine)
    binary=list(itertools.product(range(2),repeat=4));bi={v:i for i,v in enumerate(binary)}
    perm=[bi[(v[1],(v[0]+v[1])%2,v[3],(v[2]+v[3])%2)] for v in binary]
    unseen=set(range(1,16));coarse={}
    while unseen:
        a=min(unseen);cyc={a,perm[a],perm[perm[a]]}
        for x in cyc:coarse[x]=a
        unseen-=cyc
    table={tuple(int(r[x]) for x in ('v','k','lam','mu')):r['status'] for r in csv.DictReader((ROOT/'work/sweep_20261115/brouwer_live.csv').open())}
    for k in [385,420]:
        source=HERE/'results/385' if k==385 else ROOT/'work/construct420_20261116/result'
        z=np.zeros((16,41),dtype=np.int16)
        for line in (source/'connections_Z6_4.txt').open():
            v=tuple(map(int,line.split()));a=bi[tuple(x%2 for x in v)];p=np.array([x%3 for x in v],dtype=np.int16)
            if not np.any(p):z[a,0]=1
            else:z[a,pd[tuple(p*next(int(x) for x in p if x)%3)]]=1
        labels={};parts=[]
        for a in range(16):
            for j in range(41):
                if a==j==0:continue
                key=('within',int(z[a,j])) if a==0 else ('cross',coarse[a],j==0,int(z[a,j]))
                if key not in labels:labels[key]=len(parts);parts.append(np.zeros(656,dtype=np.int16))
                parts[labels[key]][41*a+j]=1
        spectra=np.array(parts)@C.T
        text=f'{len(parts)} 656\n'+'\n'.join(' '.join(map(str,row)) for row in spectra)+'\n'
        result=subprocess.run(['/tmp/srg420_fusions'],input=text,capture_output=True,text=True,check=True)
        hits=[]
        for line in result.stdout.splitlines():
            mask,degree,lam,mu=map(int,line.split());param=(1296,degree,lam,mu)
            hits.append(dict(mask=mask,parameters=param,table_status=table.get(param),locally_constructed=degree in [385,420,455,490,525]))
        out=dict(source_degree=k,pieces=len(parts),unions_examined=2**len(parts),hits=hits,parameter_rows=sorted({tuple(h['parameters']) for h in hits}))
        (HERE/f'fusion_{k}.json').write_text(json.dumps(out,indent=2)+'\n')
        print(k,'pieces',len(parts),'unions',2**len(parts),'parameter rows',out['parameter_rows'])


if __name__=='__main__':main()
