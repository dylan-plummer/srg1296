"""Package graph6 representatives, difference sets and portable verification data."""
from hashlib import sha256
import json
from zipfile import ZIP_DEFLATED, ZipFile, ZipInfo
from build import HERE,V,coset,elements,kernel,partition,sources,subgroups
def records(filename):
    return [json.loads(s) for s in (HERE/filename).read_text().splitlines()]
def translation_generators(blocks):
    cover=set(map(tuple,blocks))
    group=[t for t in V if all(tuple(coset(b,t)) in cover for b in blocks)]
    reached={0};gens=[]
    for t in group:
        if t in reached:continue
        gens.append(int(t))
        new=reached.copy()
        for _ in range(6):
            new=set(map(int,coset(list(new),t)))
            reached.update(new)
    assert set(group)==reached
    return gens,len(group)
def write():
    records10=records('class_checks.jsonl')
    info=json.loads((HERE/'census.json').read_text())
    cells=[b for z in info['free_groups'] for b in partition(elements(z),V)]
    manifest=[]
    for r,g in zip(info['records'],records10,strict=True):
        assert r['file']==g['file']
        blocks=[cells[i] for i in r['blocks']]
        gens,order=translation_generators(blocks)
        manifest.append(dict(**g,path='classes/'+r['file'],source='I.10',directions=len({t//36 for t in r['blocks']}),
                             translations=gens,translation_order=order))
    g,=records('uniform6_checks.jsonl')
    gens,order=translation_generators(partition(elements(records('uniform6.jsonl')[0]['subgroup']),V))
    manifest.append(dict(**g,path=g['file'],source='I.6',directions=1,translations=gens,translation_order=order))
    for r in manifest:
        r['four_cliques']=sum(x*n for x,n in r['K4_per_vertex'])//4
    assert len({r['four_cliques'] for r in manifest})==25
    (HERE/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    files=[r['path'] for r in manifest]+['result.txt','manifest.json','census.json','edge_certificates.jsonl',
        'audit.json','square24.py','regular.out','regular24.out','pds_2797.json','pds_2848.json','verify_pds.g',
        'pds_check.out','verify_bundle.py','rebuild_bundle.py','uniform6.jsonl','compare595.jsonl','difference_orbits.jsonl','koval_verify.log']
    assert len((HERE/'koval_verify.log').read_text().splitlines())==25
    with ZipFile(HERE/'srg1296_595_spread_completions.zip','w') as z:
        for n in (6,10):
            p=next(p for p in sources() if p.name==f'srg1296_560_244_240_class{n}.txt')
            seed=[int(x) for x in kernel(p).ravel().nonzero()[0]]
            z.writestr(ZipInfo(f'seeds/I.{n}.json',(2026,11,20,0,0,0)),json.dumps(seed)+'\n',compress_type=ZIP_DEFLATED)
        for f in files:
            z.writestr(ZipInfo(f,(2026,11,20,0,0,0)),(HERE/f).read_bytes(),compress_type=ZIP_DEFLATED)
    print(sha256((HERE/'srg1296_595_spread_completions.zip').read_bytes()).hexdigest())
if __name__=='__main__':write()
