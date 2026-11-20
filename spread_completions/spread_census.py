"""Enumerate the full set of partitions by the free subgroup cosets of class I.10."""
from collections import Counter
from pathlib import Path
from hashlib import sha256
import json
import numpy as np
import igraph as ig
from build import HERE, sources, kernel, graph, V, partition, elements, subgroups, graph6, check
from compare import canon, encoding

def bits(indices):
    return sum(1<<int(i) for i in indices)

def indices(mask):
    while mask:
        bit=mask&-mask
        yield bit.bit_length()-1
        mask^=bit

def enumerate_spreads(cells):
    blocks=[bits(b) for b in cells]
    legal=[bits(i for i,b in enumerate(blocks) if (b>>x)&1) for x in range(1296)]
    conflicts=[bits(i for i,b in enumerate(blocks) if b&m) for m in blocks]
    # Points lying in exactly the same candidate blocks impose the same equation.
    constraints=sorted(set(legal))
    hit=[bits(j for j,s in enumerate(constraints) if (s>>i)&1) for i in range(len(cells))]
    def walk(avail,missing,chosen):
        if not missing:
            assert chosen.bit_count()==36
            yield chosen
            return
        best=None
        for j in indices(missing):
            choices=constraints[j]&avail
            if not choices:return
            if best is None or choices.bit_count()<best.bit_count():
                best=choices
                if best.bit_count()==1:break
        for i in indices(best):
            yield from walk(avail&~conflicts[i],missing&~hit[i],chosen|(1<<i))
    return set(walk((1<<len(cells))-1,(1<<len(constraints))-1,0))

def host_permutations(a,cells):
    i,j=np.triu(a,1).nonzero()
    host=ig.Graph(n=len(a),edges=list(zip(i.tolist(),j.tolist())))
    lookup={bits(b):i for i,b in enumerate(cells)}
    permutations=[]
    for gen in host.automorphism_group(sh='fl'):
        permutations.append([lookup[bits(gen[x] for x in b)] for b in cells])
    return host.count_automorphisms(sh='fl'),permutations

def host_classes(all_spreads,perms):
    unseen=all_spreads.copy()
    while unseen:
        start=min(unseen)
        pending=[start]
        unseen.remove(start)
        size=1
        while pending:
            s=list(indices(pending.pop()))
            for perm in perms:
                image=bits(perm[i] for i in s)
                if image in unseen:
                    unseen.remove(image)
                    pending.append(image)
                    size+=1
        yield start,size

def main():
    p=next(p for p in sources() if p.name=='srg1296_560_244_240_class10.txt')
    d=kernel(p);a=graph(d);u=subgroups(d,0)
    cells=[s for z in u for s in partition(elements(z),V)]
    covers=enumerate_spreads(cells)
    print('covers',len(covers),flush=True)
    aut,perms=host_permutations(a,cells)
    baselines=[json.loads(l) for l in (HERE/'compare595.jsonl').read_text().splitlines() if 'mixed595.g6' not in l]
    baseline={r['bliss_sha256']:r['file'] for r in baselines}
    cert={};records=[]
    (HERE/'classes').mkdir(exist_ok=True)
    for rep,size in host_classes(covers,perms):
        chosen=list(indices(rep))
        b=a.copy()
        for block in chosen:
            ix=np.ix_(cells[block],cells[block])
            assert not b[ix].any()
            b[ix]=1
        np.fill_diagonal(b,0)
        check(b,(1296,595,274,272))
        i,j=np.triu(b,1).nonzero()
        g=ig.Graph(n=1296,edges=list(zip(i.tolist(),j.tolist())))
        key=sha256(encoding(canon(g))).hexdigest()
        if key not in cert:
            name='class_%02d.g6'%(len(cert)+1)
            cert[key]=name
            (HERE/'classes'/name).write_bytes(graph6(b))
        rec=dict(blocks=chosen,orbit=size,file=cert[key],canonical_sha256=key,baseline=baseline.get(key))
        records.append(rec)
    data=dict(source=p.name,free_groups=u,cosets=len(cells),spreads=len(covers),host_aut=aut,
              directions_histogram=sorted(Counter(len({i//36 for i in indices(rep)}) for rep in covers).items()),
              host_orbits=len(records),graph_classes=len(cert),records=records)
    (HERE/'census.json').write_text(json.dumps(data,indent=2)+'\n')
    print('host orbits',len(records),'graph classes',len(cert),'directions',data['directions_histogram'],flush=True)
if __name__=='__main__':main()
