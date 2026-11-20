"""Check sample seeds, stored graphs and pairwise distinctions."""
from hashlib import shake_256,sha256
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parent

def key(census):
    return tuple(sorted((int(k),v) for k,v in census.items()))

def read(path):
    return [json.loads(x) for x in path.read_text().splitlines()]

def main():
    baseline={}
    for p in sorted(ROOT.glob('independent_*.jsonl')):
        for r in read(p):
            _,k,_,_=r['parameters']
            assert r['sha256']==sha256((ROOT/r['file']).read_bytes()).hexdigest()
            baseline[k]=key(r['four_cliques_through_vertex'])
    result={}
    for degree,hist in sorted(baseline.items()):
        seen={hist}
        assert len(hist)>1
        if degree in (560,592):
            sample=read(ROOT/f'sample_{degree}.jsonl')
            assert len(sample)==32
            for i,r in enumerate(sample):
                label=f'SRG binary fibre switch {degree} sample {i}'.encode()
                mask=int.from_bytes(shake_256(label).digest(10),'little')<<1
                assert r['index']==i and int(r['fibre_mask'],16)==mask
                h=key(r['four_cliques_through_vertex'])
                assert len(h)>1 and sum(n for _,n in h)==1296
                assert h not in seen and sum(k*n for k,n in h)==4*r['total_four_cliques']
                assert r['adjacent_pairs']+r['nonadjacent_pairs']==839160
                seen.add(h)
            assert sample[0]['sha256']==sha256((ROOT/f'arbitrary_{degree}.g6').read_bytes()).hexdigest()
            assert sample[0]['orbit_sizes']==[16]*81 and sample[0]['automorphism_order']==16
        result[degree]=len(seen)
    assert result=={385:1,420:1,560:33,592:33}
    print(json.dumps(dict(distinct_non_cayley_classes=result,total=sum(result.values()))))

if __name__=='__main__':main()
