"""Package every certified isomorphism class as an individually named graph6 file."""
from hashlib import sha256
from pathlib import Path
from tempfile import TemporaryDirectory
from zipfile import ZipFile, ZipInfo, ZIP_DEFLATED
import json

from construct import build_cut, graph6
from family import GAUGES
from scan import ROOT

def put(z,name,data):
    info=ZipInfo(name,(1980,1,1,0,0,0))
    info.compress_type=ZIP_DEFLATED
    z.writestr(info,data)

def main():
    manifest=[]
    def record(z,name,raw,digest,param,detail):
        assert sha256(raw).hexdigest()==digest
        put(z,name,raw)
        manifest.append(dict(file=name,sha256=digest,parameters=param,**detail))
    dest=ROOT/'non_cayley_68.zip'
    with TemporaryDirectory(dir=ROOT) as tmp:
        path=Path(tmp)/dest.name
        with ZipFile(path,'w') as z:
            for k in (385,420,560,592):
                cert=json.loads((ROOT/f'switched_{k}.json').read_text())
                name=f'{k}/hyperplane.g6'
                record(z,name,(ROOT/cert['graph']).read_bytes(),cert['sha256'],cert['parameters'],
                       dict(source=cert['source'],normal=cert['normal'],negative_frequencies=cert['negative_frequencies']))
                if k not in GAUGES:continue
                spec=json.loads((ROOT/f'family_{k}.json').read_text())
                source=ROOT.parent/'connsets'/spec['source']
                for line in (ROOT/f'sample_{k}.jsonl').read_text().splitlines():
                    r=json.loads(line)
                    mask=int(r['fibre_mask'],16)
                    param,a=build_cut(source,[bool(mask>>j&1) for j in range(81)],GAUGES[k])
                    name=f"{k}/sample_{r['index']:02d}.g6"
                    record(z,name,graph6(a),r['sha256'],param,
                           dict(source=source.name,fibre_mask=r['fibre_mask'],
                                negative_frequencies=spec['negative_frequencies'],
                                four_cliques_through_vertex=r['four_cliques_through_vertex']))
            assert len(manifest)==68
            put(z,'manifest.json',(json.dumps(manifest,indent=2)+'\n').encode())
        with ZipFile(path,'r') as z:
            assert z.testzip() is None and len(z.namelist())==69
            assert all(sha256(z.read(r['file'])).hexdigest()==r['sha256'] for r in manifest)
        path.replace(dest)
    print(json.dumps(dict(archive=dest.name,graphs=len(manifest),sha256=sha256(dest.read_bytes()).hexdigest())))

if __name__=='__main__':main()
