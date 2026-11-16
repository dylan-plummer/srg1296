"""Bounded degree-490 transfer using one size-36 and two size-27 seeds."""
from join455 import *


def local_arrays(D,file_index):
    masks={27:[],36:[]}
    for line in (HERE/f'local_sets_{file_index}.tsv').open():
        size,mask=line.split();masks[int(size)].append(int(mask,16))
    planes=np.array([sum(int(a)<<j for j,a in enumerate(row)) for row in D['N']],dtype=np.uint64)
    result={}
    for size,values in masks.items():
        values=np.array(values,dtype=np.uint64);C=np.empty((len(values),40),dtype=np.int8)
        for first in range(0,len(values),10000):
            hits=np.bitwise_count(values[first:first+10000,None]&planes).astype(np.int16)
            shift=4 if size==27 else 6;assert np.all((hits-shift)%3==0);C[first:first+10000]=(hits-shift)//3
        result[size]=(values,C)
    return result


def join(arrays,templates,budget,deadline):
    A,B,C=[a['codes'] for a in arrays];abit=[bitsets(A[:,f]) for f in range(8)];bbit=[bitsets(B[:,f]) for f in range(8)]
    allowed=[];left_by_cap=[];right_by_cap=[]
    def union(values,index):return functools.reduce(int.__or__,(index.get(v,0) for v in values),0)
    for f,ts in enumerate(templates):
        pairs=collections.defaultdict(set);left=collections.defaultdict(set);right=collections.defaultdict(set)
        for a,b,c in ts:pairs[int(a),int(c)].add(int(b));left[int(c)].add(int(a));right[int(c)].add(int(b))
        allowed.append({key:union(vals,bbit[f]) for key,vals in pairs.items()})
        left_by_cap.append({key:union(vals,abit[f]) for key,vals in left.items()})
        right_by_cap.append({key:union(vals,bbit[f]) for key,vals in right.items()})
    energy_bit=bitsets(arrays[1]['energy']);tested=0
    for ci,cap in enumerate(C):
        liveA=(1<<len(A))-1;liveB=(1<<len(B))-1
        for f in range(8):liveA &=left_by_cap[f].get(int(cap[f]),0);liveB &=right_by_cap[f].get(int(cap[f]),0)
        if not liveA or not liveB:continue
        order=sorted(range(8),key=lambda f:right_by_cap[f].get(int(cap[f]),0).bit_count())
        for ai in indices(liveA):
            tested+=1;alive=liveB &energy_bit.get(budget-int(arrays[0]['energy'][ai])-int(arrays[2]['energy'][ci]),0)
            for f in order:
                alive &=allowed[f].get((int(A[ai,f]),int(cap[f])),0)
                if not alive:break
            if alive:return dict(found=[ai,next(indices(alive)),ci],complete_join=False,pairs_tested=tested)
        if time.monotonic()>deadline:return dict(found=None,complete_join=False,pairs_tested=tested,stopped_on_time_limit=True)
    return dict(found=None,complete_join=True,pairs_tested=tested)


def main():
    started=time.monotonic();results=[]
    for mask,file_index in [(29,1),(51,2)]:
        D,frequency,e,templates,z=prepare(mask);local=local_arrays(D,file_index)
        for cycle in [c for c in frequency if e[c[0]]==4]:
            arrays=[];before=[]
            for b,size in enumerate([36,27,27]):
                masks,C=local[size];target=np.array([2,0,0,0,0]) if b==2 else np.zeros(5,dtype=np.int16)
                selected=np.flatnonzero(np.all(C[:,cycle]==target,axis=1));before.append(len(selected))
                cc=C[selected];codes=np.column_stack([(cc[:,orbit]+2)@np.array([5**i for i in range(5)]) for orbit in frequency])
                keep=np.ones(len(cc),dtype=bool)
                for f,ts in enumerate(templates):keep &=np.isin(codes[:,f],ts[:,b])
                arrays.append(dict(masks=masks[selected][keep],codes=codes[keep],energy=(cc[keep][:,e==4].astype(np.int16)**2).sum(1)))
            case=dict(fibre_mask=mask,matching_cycle=cycle,local_counts_before_frame_filter=before,domain_sizes=[len(a['masks']) for a in arrays])
            print(json.dumps(case),flush=True);join_start=time.monotonic();case.update(join(arrays,templates,16,join_start+20));case['join_seconds']=time.monotonic()-join_start
            case['seconds_since_start']=time.monotonic()-started
            if case['found']:case['certificate']=save_found(D,frequency,e,arrays,case['found'],z,dict(fibre_mask=mask,matching_cycle=cycle),[1296,490,190,182],mask)
            results.append(case);(HERE/'join490_results.json').write_text(json.dumps(results,indent=2)+'\n');print(json.dumps(case),flush=True)
            if case['found'] or case.get('stopped_on_time_limit') or time.monotonic()-started>45:return


if __name__=='__main__':main()
