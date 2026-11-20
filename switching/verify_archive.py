"""Cross-check each archive member with Kovaľ's unmodified independent verifier."""
from concurrent.futures import ProcessPoolExecutor
from hashlib import sha256
from importlib.util import module_from_spec,spec_from_file_location
from pathlib import Path
from zipfile import ZipFile
import json

ROOT=Path(__file__).resolve().parent
import os
# A. Kovaľ's verify_srg.py (not included here); set KOVAL_VERIFY_SRG to its path
SOURCE=Path(os.environ.get('KOVAL_VERIFY_SRG', ROOT.parents[1]/'koval/new_srgs_twoorbit_data/verify_srg.py'))
spec=spec_from_file_location('koval_verify',SOURCE)
verifier=module_from_spec(spec)
spec.loader.exec_module(verifier)

def check_member(r):
    with ZipFile(ROOT/'non_cayley_68.zip') as z:
        data=z.read(r['file'])
    assert sha256(data).hexdigest()==r['sha256']
    n,adj=verifier.decode(data.decode())
    report=verifier.check(n,adj)
    expected='OK srg('+','.join(map(str,r['parameters']))+');'
    assert report.startswith(expected),(r['file'],report)
    return dict(file=r['file'],report=report)

if __name__=='__main__':
    with ZipFile(ROOT/'non_cayley_68.zip') as z:manifest=json.loads(z.read('manifest.json'))
    with ProcessPoolExecutor(max_workers=6) as pool:
        for result in pool.map(check_member,manifest):print(json.dumps(result),flush=True)
