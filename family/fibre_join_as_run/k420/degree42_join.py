"""Weaker, fixed-point-free C3 join for SRG(1296,420,144,132).

There is one 42-element zero fiber, one developed 18-element seed and
four developed 27-element seeds. No C3 x C3 symmetry is imposed.
"""
import argparse
import collections
import itertools
import json
import time
from pathlib import Path

import numpy as np
from c3_join import centralizer, support_index, bits

HERE = Path(__file__).resolve().parent


def cycles(perm):
    unseen = set(range(len(perm)))
    out = []
    while unseen:
        i = min(unseen)
        cyc = [i]
        while perm[cyc[-1]] != i:
            cyc.append(perm[cyc[-1]])
        unseen -= set(cyc)
        out.append(cyc)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--case', choices=['22', '31'], default='31')
    ap.add_argument('--seconds', type=float, default=45)
    ap.add_argument('--nodes', type=int, default=10000)
    ap.add_argument('--geometric18', action='store_true')
    ap.add_argument('--matching27', action='store_true')
    ap.add_argument('--rank8-matching', action='store_true')
    ap.add_argument('--two-caps', action='store_true')
    args = ap.parse_args()
    start = time.monotonic()
    V = np.array(list(itertools.product(range(3), repeat=4)), dtype=np.int16)
    vd = {tuple(v): i for i, v in enumerate(V)}
    P = np.array([v for v in V if np.any(v) and next(int(x) for x in v if x) == 1], dtype=np.int16)
    pd = {tuple(p): i for i, p in enumerate(P)}
    def normalize(v):
        return tuple(v*next(int(x) for x in v if x) % 3)
    N = P @ P.T % 3 == 0
    char = np.where(N, 2, -1).astype(np.int16)
    plane_masks = np.array([sum(int(x) << j for j, x in enumerate(row)) for row in N], dtype=np.uint64)
    H = np.array([[(-1)**((a & b).bit_count() % 2) for b in range(16)] for a in range(16)], dtype=np.int16)
    binary = list(itertools.product(range(2), repeat=4))
    bi = {v: i for i, v in enumerate(binary)}
    coarse = [bi[(v[1], (v[0]+v[1]) % 2, v[3], (v[2]+v[3]) % 2)] for v in binary]
    quotient_sizes = np.array([42] + [18 if b % 4 == 0 else 27 for b in range(1, 16)], dtype=np.int16)
    qev = H @ quotient_sizes
    assert qev[0] == 420 and np.all(np.isin(qev[1:], [-12, 24]))
    cc = [c for c in cycles(coarse) if c != [0]]
    cc.sort(key=lambda c: (quotient_sizes[c[0]], c))
    for i in range(1,5):
        j = next(j for j,a in enumerate(cc[i]) if a % 4 == 1)
        cc[i] = cc[i][j:]+cc[i][:j]
    cc[1:] = sorted(cc[1:])
    assert len(cc) == 5 and [int(quotient_sizes[c[0]]) for c in cc] == [18,27,27,27,27]
    if args.case == '22':
        G = np.array([[1,0,1,0], [0,1,0,1], [0,0,1,0], [0,0,0,1]], dtype=np.int16)
    else:
        G = np.array([[1,1,0,0], [0,1,1,0], [0,0,1,0], [0,0,0,1]], dtype=np.int16)
    primal = [pd[normalize(G @ p % 3)] for p in P]
    dual = [pd[normalize(G.T @ p % 3)] for p in P]
    frequency = sorted(cycles(dual), key=lambda c: (len(c), c))
    po = [sum(1 << j for j in c) for c in cycles(primal)]
    zero_masks, zero_eigen = [], []
    for subset in range(1 << len(po)):
        mask = sum(m for j, m in enumerate(po) if subset >> j & 1)
        if mask.bit_count() != 21:
            continue
        e = np.array([3*(mask & int(pm)).bit_count()-21 for pm in plane_masks], dtype=np.int16)
        if np.all(np.isin(e, [-12,-3,6,15,24])):
            zero_masks.append(mask)
            zero_eigen.append(e)
    zero_masks = np.array(zero_masks, dtype=np.uint64)
    zero_eigen = np.array(zero_eigen, dtype=np.int16).reshape(-1,40)
    masks18 = np.array([int(line.split()[1],16) for line in (HERE/'local18_0.tsv').open()], dtype=np.uint64)
    assert len(masks18) == 9880 and len(np.unique(masks18)) == 9880
    hits18 = np.bitwise_count(masks18[:,None] & plane_masks).astype(np.int16)
    assert np.all((hits18-3) % 3 == 0)
    C18 = ((hits18-3)//3).astype(np.int8)
    if args.geometric18:
        keep = np.any(C18 == 2, axis=1)
        masks18, C18 = masks18[keep], C18[keep]
        assert len(masks18) == 520
    z = np.load(HERE/'local27.npz')
    masks27, C27 = z['masks'], z['coefficients']
    patterns = []
    for rank in (0,4,8,12,16):
        for S in itertools.combinations(range(16), rank):
            xor = 0
            for x in S:
                xor ^= x
            if xor:
                continue
            w = H[:,S].sum(axis=1) if rank else np.zeros(16,dtype=np.int16)
            assert np.all(w % 4 == 0)
            w = w//4
            assert np.all(np.isin(H @ w, [0,4]))
            if args.rank8_matching and rank == 8 and not (np.count_nonzero(w[1:]) == 1 and np.max(w[1:]) == 2):
                continue
            patterns.append(w)
    relations = []
    for fo in frequency:
        rows = []
        for w in patterns:
            if len(fo) == 1 and any(len({int(w[a]) for a in cyc}) != 1 for cyc in cc):
                continue
            rows.append([int(w[0])] + [sum((int(w[a])+2)*5**j for j,a in enumerate(cyc[:len(fo)])) for cyc in cc])
        relations.append(np.array(sorted(set(map(tuple, rows))),dtype=np.int16))
    def encode(C):
        return np.column_stack([(C[:,fo].astype(np.int16)+2) @ np.array([5**j for j in range(len(fo))],dtype=np.int16)
                                for fo in frequency])
    masks = [zero_masks, masks18] + [masks27]*4
    codes = [np.column_stack([(zero_eigen[:,fo[0]]+12)//9 for fo in frequency]), encode(C18)] + [encode(C27)]*4
    # Directly reconstruct random block models, independently of the code
    # transport, and compare every represented row and its C3 transports.
    transport_checks = 0
    if len(zero_masks):
        for trial in range(3):
            choice = [(trial*37+11) % len(ms) for ms in masks]
            T = np.zeros((16,81), dtype=np.int16)
            def put(a, mask, zero):
                T[a,0] = int(zero)
                for j,p in enumerate(P):
                    if mask >> j & 1:
                        T[a,vd[tuple(p)]] = T[a,vd[tuple(-p % 3)]] = 1
            put(0, int(masks[0][choice[0]]), False)
            for v,cyc in enumerate(cc,1):
                mask = int(masks[v][choice[v]])
                for a in cyc:
                    put(a,mask,v>1)
                    mask = sum(1 << primal[j] for j in range(40) if mask >> j & 1)
            assert np.array_equal(T.sum(axis=1),quotient_sizes)
            F = T[:,0,None]+T[:,[vd[tuple(p)] for p in P]] @ char.T
            for f,fo in enumerate(frequency):
                b = fo[0]
                direct = [(int(F[0,b])+12)//9] + [sum((int(F[a,b])//9+2)*5**j for j,a in enumerate(cyc[:len(fo)])) for cyc in cc]
                assert direct == [int(codes[v][choice[v],f]) for v in range(6)]
                ok = np.all(np.isin(H @ (F[:,b]+12*np.eye(16,dtype=np.int16)[0]),[0,36]))
                if args.rank8_matching and F[0,b] == 6:
                    ok = ok and np.count_nonzero(F[1:,b]) == 1 and np.max(F[1:,b]) == 18
                assert ok == (tuple(direct) in set(map(tuple,relations[f])))
                for q in fo:
                    assert np.array_equal(F[np.array(coarse),q], F[:,dual[q]])
                transport_checks += 1
    # Join all four fixed fine characters in one compatibility relation.
    nf = sum(len(fo)==1 for fo in frequency)
    assert nf == 4
    combined = np.zeros((1,6),dtype=np.int16)
    super_counts = []
    for f in range(nf):
        combined = (combined[:,None,:]+relations[f][None,:,:]*5**f).reshape(-1,6)
        keep = np.ones(len(combined),dtype=bool)
        for v in range(6):
            prefix = codes[v][:,:f+1] @ np.array([5**j for j in range(f+1)],dtype=np.int16)
            keep &= np.isin(combined[:,v],np.unique(prefix))
        combined = np.unique(combined[keep],axis=0)
        super_counts.append(len(combined))
    relations = [combined]+relations[nf:]
    def merge(c):
        return np.column_stack([c[:,:nf] @ np.array([5**j for j in range(nf)],dtype=np.int16),c[:,nf:]])
    codes = [merge(c) for c in codes[:2]]+[merge(codes[2])]*4
    supports = [[support_index(c[:,f]) for f in range(len(relations))] for c in codes[:3]]
    supports += [supports[2]]*3
    rsup = [[support_index(R[:,v]) for v in range(6)] for R in relations]
    rfull = [(1 << len(R))-1 for R in relations]
    domains = [(1 << len(c))-1 for c in codes]
    if args.matching27:
        domains[2] = int.from_bytes(np.packbits(np.any(C27 == 2,axis=1),bitorder='little').tobytes(),'little')
    if args.two_caps:
        assert args.rank8_matching and args.geometric18
        capmask = int.from_bytes(np.packbits(np.any(C27 == 2,axis=1),bitorder='little').tobytes(),'little')
        domains[2] = domains[3] = capmask
        domains[4] &= ~capmask
        domains[5] &= ~capmask
    # Normalize T0 under the full common fine centralizer. Other seed
    # domains are complete and transported by this same coordinate change.
    group = centralizer(G)
    permutations = [[pd[normalize(K @ p % 3)] for p in P] for K in group]
    lookup0 = {int(mask):i for i,mask in enumerate(zero_masks)}
    remaining = set(range(len(zero_masks)))
    zero_reps = []
    orbit_sizes = []
    while remaining:
        i = min(remaining)
        zero_reps.append(i)
        mask = int(zero_masks[i])
        orbit = {lookup0[sum(1 << pm[j] for j in range(40) if mask >> j & 1)] for pm in permutations}
        remaining -= orbit
        orbit_sizes.append(len(orbit))
    stats = dict(parameters=[1296,420,144,132],zero_fibre_degree=42,case=args.case,
                 geometric18=args.geometric18,matching27=args.matching27,
                 rank8_matching=args.rank8_matching,two_caps=args.two_caps,
                 coarse=coarse,coarse_cycles=cc,frequency_cycles=frequency,
                 initial_domains=[len(c) for c in codes],fixed_character_join_rows=super_counts,
                 relation_rows=[len(R) for R in relations],fine_centralizer_order=len(group),
                 zero_fibre_orbit_sizes=orbit_sizes,transport_controls=transport_checks,
                 zero_fibre_orbit_spectra=[dict(collections.Counter(map(int,zero_eigen[i]))) for i in zero_reps],
                 preprocessing_seconds=time.monotonic()-start)
    print(json.dumps(stats),flush=True)
    deadline = time.monotonic()+args.seconds
    nodes = 0
    root_info = None
    def propagate(ds):
        active = rfull.copy()
        change = True
        rounds = 0
        while change:
            change = False
            rounds += 1
            for f,R in enumerate(relations):
                alive = active[f]
                for v in range(6):
                    allowed = 0
                    for c,rmask in rsup[f][v].items():
                        if ds[v] & supports[v][f].get(c,0):
                            allowed |= rmask
                    alive &= allowed
                    if not alive:
                        return None,None,dict(rounds=rounds,empty_relation=f)
                active[f] = alive
                for v in range(6):
                    allowed = 0
                    for c,rmask in rsup[f][v].items():
                        if rmask & alive:
                            allowed |= supports[v][f].get(c,0)
                    new = ds[v] & allowed
                    if not new:
                        return None,None,dict(rounds=rounds,empty_variable=v,relation=f)
                    if new != ds[v]:
                        ds[v] = new
                        change = True
        return ds,active,dict(rounds=rounds,domains=[d.bit_count() for d in ds],relations=[r.bit_count() for r in active])
    class Limit(Exception):
        pass
    def search(ds,depth=0):
        nonlocal nodes
        if time.monotonic()>deadline or nodes>=args.nodes:
            raise Limit
        nodes += 1
        ds,active,info = propagate(ds)
        if ds is None:
            return None
        sizes = [(d.bit_count(),v) for v,d in enumerate(ds) if d.bit_count()>1]
        if not sizes:
            return [next(bits(d)) for d in ds]
        choices = [(r.bit_count(),f) for f,r in enumerate(active) if r.bit_count()>1]
        if not choices:
            return [next(bits(d)) for d in ds]
        count,v = min(sizes)
        rcount,f = min(choices)
        if count <= min(128,rcount):
            for i in bits(ds[v]):
                child = ds.copy()
                child[v] = 1 << i
                found = search(child,depth+1)
                if found:
                    return found
        else:
            for ri in bits(active[f]):
                row = relations[f][ri]
                child = [ds[v] & supports[v][f].get(int(row[v]),0) for v in range(6)]
                assert child != ds
                found = search(child,depth+1)
                if found:
                    return found
        return None
    found = None
    done = 0
    pair_counts = []
    try:
        root,_,root_info = propagate(domains.copy())
        print(json.dumps(dict(stage='root_propagated',**root_info)),flush=True)
        if root is not None:
            for zi in zero_reps:
                if not (root[0] >> zi & 1):
                    continue
                ds = root.copy()
                ds[0] = 1 << zi
                ds,_,info = propagate(ds)
                if ds is None:
                    continue
                zmask = int(zero_masks[zi])
                stab = [pm for pm in permutations if sum(1 << pm[j] for j in range(40) if zmask >> j & 1) == zmask]
                live18 = set(bits(ds[1]))
                lookup18 = {int(masks18[i]):i for i in live18}
                remaining18 = live18.copy()
                reps18 = []
                while remaining18:
                    i = min(remaining18)
                    reps18.append(i)
                    mask = int(masks18[i])
                    orbit = {lookup18[sum(1 << pm[j] for j in range(40) if mask >> j & 1)] for pm in stab}
                    assert orbit <= live18
                    remaining18 -= orbit
                pair_counts.append(dict(zero_index=zi,stabilizer=len(stab),raw18=len(live18),representatives18=len(reps18)))
                print(json.dumps(dict(stage='pair_normalized',**pair_counts[-1])),flush=True)
                for si in reps18:
                    child = ds.copy()
                    child[1] = 1 << si
                    found = search(child)
                    done += 1
                    if found:
                        break
                if found:
                    break
        status = 'FOUND' if found else 'INFEASIBLE_IN_SPECIFIED_MODEL'
    except Limit:
        status = 'UNKNOWN_LIMIT'
    stats.update(status=status,root_propagation=root_info,nodes=nodes,fixed_pairs_completed=done,pair_normalization=pair_counts,
                 search_seconds=time.monotonic()-start-stats['preprocessing_seconds'],total_seconds=time.monotonic()-start)
    if found:
        chosen_masks = [int(ms[i]) for ms,i in zip(masks,found)]
        T = np.zeros((16,81),dtype=np.int16)
        def put(a,mask,zero):
            T[a,0] = int(zero)
            for j,p in enumerate(P):
                if mask >> j & 1:
                    T[a,vd[tuple(p)]] = T[a,vd[tuple(-p % 3)]] = 1
        put(0,chosen_masks[0],False)
        for v,cyc in enumerate(cc,1):
            mask = chosen_masks[v]
            for a in cyc:
                put(a,mask,v>1)
                mask = sum(1 << primal[j] for j in range(40) if mask >> j & 1)
        F = T[:,0,None]+T[:,[vd[tuple(p)] for p in P]] @ char.T
        assert np.array_equal(T.sum(axis=1),quotient_sizes)
        assert np.all(np.isin(H @ F,[-12,24]))
        np.savez_compressed(HERE/f'FOUND_degree42_{args.case}.npz',T=T,G=G,P=P,
                            coarse=np.array(coarse),selected_masks=np.array(chosen_masks,dtype=np.uint64))
        stats['selected_indices'] = found
        stats['selected_masks'] = [hex(x) for x in chosen_masks]
    suffix = '_geometric18' if args.geometric18 else ''
    if args.matching27:
        suffix += '_matching27'
    if args.rank8_matching:
        suffix += '_rank8matching'
    if args.two_caps:
        suffix += '_twocaps'
    (HERE/f'degree42_{args.case}{suffix}_result.json').write_text(json.dumps(stats,indent=2)+'\n')
    print(json.dumps(stats),flush=True)


if __name__ == '__main__':
    main()
