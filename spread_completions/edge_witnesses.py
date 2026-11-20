"""Certify an edge outside every 36-clique by coloring its common neighborhood."""
import json
import numpy as np
from compare import HERE,read_graph
def coloring(rows,vertices):
    adj={v:{u for u in vertices if rows[v]>>u&1} for v in vertices}
    left=set(vertices);colors={};used={v:set() for v in vertices}
    while left:
        v=max(left,key=lambda x:(len(used[x]),len(adj[x]&left),-x));left.remove(v)
        c=0
        while c in used[v]:c+=1
        colors[v]=c
        for u in adj[v]&left:used[u].add(c)
    return [[v for v in vertices if colors[v]==c] for c in range(max(colors.values())+1)]
def reduce_colors(rows,vertices,classes,target=33):
    n=len(vertices)
    index={v:i for i,v in enumerate(vertices)}
    neighbors=[np.array([index[u] for u in vertices if rows[v]>>u&1]) for v in vertices]
    labels=np.full(n,-1)
    counts=np.zeros((n,target),dtype=int)
    tabu=np.zeros_like(counts)
    rng=np.random.default_rng(595)
    for group in sorted(classes,key=len,reverse=True):
        for v in group:
            x=index[v]
            best=np.flatnonzero(counts[x]==counts[x].min())
            c=rng.choice(best)
            labels[x]=c
            counts[neighbors[x],c]+=1
    ix=np.arange(n)
    best=int(counts[ix,labels].sum())//2
    for t in range(12000):
        current=counts[ix,labels]
        bad=ix[current>0]
        score=int(current.sum())//2
        if not score:
            return [[vertices[i] for i in ix[labels==c]] for c in range(target)]
        delta=counts[bad]-current[bad,None]
        delta[(tabu[bad]>t)&(score+delta>=best)]=n
        delta[np.arange(len(bad)),labels[bad]]=n
        ys,cs=np.where(delta==delta.min())
        pick=rng.integers(len(ys))
        x,c=bad[ys[pick]],cs[pick]
        old=labels[x]
        counts[neighbors[x],old]-=1
        counts[neighbors[x],c]+=1
        tabu[x,old]=t+int(.6*len(bad))+int(rng.integers(1,11))
        labels[x]=c
        best=min(best,score)
    return None
def main():
    for p in sorted((HERE/'classes').glob('*.g6'))+sorted(HERE.glob('uniform6_*.g6')):
        rows,_=read_graph(p)
        for j in range(1,len(rows)):
            if not rows[0]>>j&1:continue
            common=rows[0]&rows[j];vertices=[i for i in range(len(rows)) if common>>i&1]
            classes=coloring(rows,vertices)
            if len(classes)>33:
                classes=reduce_colors(rows,vertices,classes)
            if classes is not None:
                assert len(vertices)==274 and sorted(sum(classes,[]))==vertices
                assert all(not (rows[v]>>w&1) for block in classes for v in block for w in block)
                print(json.dumps(dict(file=p.name,edge=[0,j],colors=len(classes),partition=classes)),flush=True)
                break
        else:raise RuntimeError(p)
if __name__=='__main__':main()
