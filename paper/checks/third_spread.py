#!/usr/bin/env python3
r"""Kovaľ's 'third spread' question for our LS_16(36) graphs: for a (1296,560,244,240) PDS D, find subgroups
U of order 36 with U & D = {} and test whether D+U* (LS_17), D+U*+V* (LS_18), D+U*+V*+W* (LS_19) are PDSs,
for pairwise orthogonal U, V, W (pairwise intersections {0})."""
import os, glob, itertools
import numpy as np
from spread_structure import load, Us, is_pds, CONN, SHAPE
zero = np.zeros(SHAPE, bool); zero[(0,) * 8] = True
LS = {17: (19, -17), 18: (18, -18), 19: (17, -19)}          # LS_m(36): eigenvalues 36-m, -m
for fn in sorted(glob.glob(os.path.join(CONN, "srg1296_560_*_class*.txt")), key=lambda f: int(f.split("_class")[1][:-4])):
    prm, D = load(fn)
    free = [U for U, a, b in Us if not (U & D).any()]
    if not free:
        continue
    one = [U for U in free if is_pds(D | (U & ~zero), *LS[17])]
    orth = lambda X, Y: (X & Y).sum() == 1
    two = [(U, V) for U, V in itertools.combinations(one, 2) if orth(U, V) and is_pds(D | ((U | V) & ~zero), *LS[18])]
    three = [(U, V, W) for U, V, W in itertools.combinations(one, 3)
             if orth(U, V) and orth(U, W) and orth(V, W) and is_pds(D | ((U | V | W) & ~zero), *LS[19])]
    print(f"{os.path.basename(fn)[8:-4]}: free subgroups {len(free)}, D+U* is LS17: {len(one)}, "
          f"orthogonal pairs giving LS18: {len(two)}, orthogonal triples giving LS19: {len(three)}")
