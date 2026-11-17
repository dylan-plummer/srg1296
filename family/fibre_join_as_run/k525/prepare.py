from pathlib import Path
import json
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]

def main():
    d=json.loads((ROOT/'work/reprioritise_20261115/c5_projective_reduction.json').read_text())
    P=d['points'];matching=d['matching_projective_channels'];lines=['40']
    for p,diagonal in zip(P,d['diagonal_fourier']):
        hits=[int(sum(a*b for a,b in zip(p,P[j]))%3==0) for j in matching]
        lines.append(' '.join(map(str,p+hits+[diagonal])))
    (HERE/'points.txt').write_text('\n'.join(lines)+'\n')
    (HERE/'source.json').write_text(json.dumps(d)+'\n')
    print('Prepared forty projective points and five matching-plane constraints.')

if __name__=='__main__':main()
