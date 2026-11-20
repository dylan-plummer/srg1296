"""Run with the zip path to verify graph parameters, clique counts and non-OA certificates."""
from collections import Counter
from hashlib import sha256
import json,sys
from zipfile import ZipFile

def decode(raw):
    enc=raw.rstrip(b'\n')
    assert enc[:4]==bytes([126,63,83,79]) and len(enc)==4+839160//6
    values=bytes(x-63 for x in enc[4:])
    assert max(values)<64
    rows=[0]*1296
    bit=0
    for j in range(1296):
        for i in range(j):
            if values[bit//6]&(32>>(bit%6)):
                rows[i]|=1<<j
                rows[j]|=1<<i
            bit+=1
    return rows
def translation(t):
    def tern(y):
        return sum(3**j*((y//(3**j)+t%81//(3**j))%3) for j in range(4))
    table=[tern(y) for y in range(81)]
    return [81*((i//81)^(t//81))+table[i%81] for i in range(1296)]
def check_translation(rows,p):
    tables=[]
    for a in range(0,len(rows),8):
        tab=[0]*256
        for m in range(1,256):
            bit=m&-m
            tab[m]=tab[m^bit]|(1<<p[a+bit.bit_length()-1])
        tables.append(tab)
    for i,r in enumerate(rows):
        assert sum(tab[b] for tab,b in zip(tables,r.to_bytes(162,'little')))==rows[p[i]]
def clique_count(rows,v):
    remaining=rows[v];cnt=0
    while remaining:
        x=remaining&-remaining;remaining^=x
        linked=remaining&rows[x.bit_length()-1]
        while linked:
            y=linked&-linked;linked^=y
            cnt+=(linked&rows[y.bit_length()-1]).bit_count()
    return cnt
def verify(z,r,witness):
    raw=z.read(r['path'])
    assert sha256(raw).hexdigest()==r['sha256']
    rows=decode(raw)
    assert all(m.bit_count()==595 for m in rows)
    for j,mask in enumerate(rows):
        for i in range(j):
            assert (rows[i]&mask).bit_count()==272+2*((mask>>i)&1)
    i,j=witness['edge'];assert rows[i]>>j&1
    classes=witness['partition'];assert len(classes)<=33
    vertices=sum(classes,[])
    assert len(vertices)==len(set(vertices))==274 and sum(1<<x for x in vertices)==rows[i]&rows[j]
    for block in classes:
        m=sum(1<<x for x in block)
        assert all(not (rows[x]&m) for x in block)
    perms=[translation(t) for t in r['translations']]
    for p in perms:check_translation(rows,p)
    left=set(range(1296));counts=Counter()
    while left:
        v=min(left);orbit={v};todo=[v]
        while todo:
            x=todo.pop()
            for p in perms:
                y=p[x]
                if y not in orbit:orbit.add(y);todo.append(y)
        assert len(orbit)==r['translation_order'] and orbit<=left
        left-=orbit
        counts[clique_count(rows,v)]+=len(orbit)
    assert sorted(map(list,counts.items()))==r['K4_per_vertex']
    total=sum(k*n for k,n in counts.items())
    assert total==4*r['four_cliques']
    return total//4
if __name__=='__main__':
    with ZipFile(sys.argv[1]) as z:
        manifest=json.loads(z.read('manifest.json'))
        witnesses={a['file']:a for a in map(json.loads,z.read('edge_certificates.jsonl').splitlines())}
        totals=[verify(z,r,witnesses[r['file']]) for r in manifest]
        assert len(totals)==len(set(totals))==25
        print('Verified 25 graph6 files, all 20,979,000 unordered pairs, 25 distinct four-clique totals, and 25 non-OA certificates.')
