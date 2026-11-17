"""Use phase-retrieval uniqueness to enforce full signed spectral templates.

After unsigned difference consistency, the local code catalogue has at
most one cross-set for each absolute Fourier vector. Hence a shared-frame
solution is translation invariant on the binary fibres. Its remaining
signed spectra must come from the 140 or 870 integral Walsh subsets.
"""
from unsigned_domains import *


def main():
    started=time.monotonic();D=json.loads((HERE/'source.json').read_text());U=json.loads((HERE/'unsigned_domains.json').read_text())
    diffcycles=U['difference_cycles'];freqcycles=U['frequency_cycles'];dual=U['action']['dual'];inv=np.argsort(dual)
    H=np.array([[(-1)**((i&j).bit_count()%2) for j in range(16)] for i in range(16)],dtype=np.int16)
    full={4:[],8:[]}
    for size in [4,8]:
        for S in itertools.combinations(range(16),size):
            xor=0
            for a in S:xor^=a
            if xor:continue
            v=H[1:,S].sum(1);assert np.all(v%4==0)
            full[size].append((v//4).tolist())
    assert len(full[4])==140 and len(full[8])==870
    templates=[]
    for fo in freqcycles:
        p=fo[0]
        if p in D['matching_projective_channels']:
            d=D['binary_cap'][D['matching_projective_channels'].index(p)]
            patterns=[[2*int(a==d) for a in range(1,16)]]
        else:patterns=full[4 if D['diagonal_fourier'][p]==-6 else 8]
        templates.append(sorted({tuple(sum((pattern[d-1]+2)*5**k for k,d in enumerate(dc)) for dc in diffcycles) for pattern in patterns}))
    outputs=[]
    arrays={size:np.load(HERE/f'local{size}.npz') for size in [27,36]}
    for b,dc in enumerate(diffcycles):
        base=np.load(HERE/f'bucket{b}.npz');size=U['buckets'][b]['cross_size']
        # This is the uniqueness condition needed for the Cayley reduction.
        assert len(base['codes'])==len(base['local_indices'])
        assert len(np.unique(base['option_indices']))==len(base['option_indices'])
        indices=base['local_indices'];C=arrays[size]['coefficients'][indices]
        keep=np.ones(len(C),dtype=bool);codes=np.empty((len(C),8),dtype=np.int16)
        for f,fo in enumerate(freqcycles):
            p=fo[0];ids=[]
            for _ in range(5):ids.append(int(p));p=inv[p]
            codes[:,f]=(C[:,ids].astype(np.int16)+2)@np.array([5**i for i in range(5)],dtype=np.int16)
            keep &= np.isin(codes[:,f],[t[b] for t in templates[f]])
        selected=np.flatnonzero(keep)
        np.savez_compressed(HERE/f'signed_bucket{b}.npz',codes=codes[keep],local_indices=indices[keep])
        outputs.append(dict(index=b,differences=dc,cross_size=size,signed_options=len(selected)))
    out=dict(scope='Shared-frame path after verified local absolute-spectrum uniqueness',
             absolute_spectrum_injective=True,rank4_signed_templates=140,rank8_signed_templates=870,
             buckets=outputs,frequency_cycles=freqcycles,difference_cycles=diffcycles,
             frequency_templates=[[list(t) for t in ts] for ts in templates],
             seconds=time.monotonic()-started)
    (HERE/'signed_domains.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='frequency_templates'},indent=2),flush=True)


if __name__=='__main__':main()
