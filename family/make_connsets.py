#!/usr/bin/env python3
"""Write the new connection sets of the paper's Section 6 into ../connsets:
 * IV.1 = srg(1296,525,216,210), V.1 = srg(1296,490,190,182), VI.1 = srg(1296,455,166,156): the imported sets
   (import_c5.py), mapped to their canonical closed form (canon.py);
 * the dual sets D* = {(c,d) != 0 : chi_(c,d)(D) = theta} of every connection set whose dual is not isomorphic to
   the original (list DUALS below, determined by family_structure.py), where chi_(c,d)(a,b) = (-1)^Tr(ca) w^Tr(db) and
   theta is the eigenvalue of multiplicity k.
usage: make_connsets.py RAW525 RAW490 RAW455"""
import sys, os
import numpy as np
from fields import *
from canon import canonical, rebuild
from duality import dual_set

HERE = os.path.dirname(os.path.abspath(__file__)); CONN = os.path.join(HERE, '..', 'connsets')
names = {525: ('525_216_210', 'IV.1'), 490: ('490_190_182', 'V.1'), 455: ('455_166_156', 'VI.1')}
for path in sys.argv[1:4]:
    params, D = read_connset(path); k = params[1]
    key, t, L = canonical(D)
    E = {(Z4 if i is None else A16[i], Z4 if j is None else B81[j]) for i, j in L}
    write_connset(os.path.join(CONN, f'srg1296_{names[k][0]}_class1.txt'), params, E,
                  f'graph {names[k][1]}: C5 fibre-join construction (16 Nov 2026), canonical form J0={key[0]} I0={key[1]}\n'
                  f'invariant under (alpha^12, beta^16) and -1; closed form in family/canon.py')
    print(names[k][1], 'J0', key[0], 'I0', key[1], 'S', [len(s) for s in key[2]])
DUALS = [('560_244_240', 5, 11), ('560_244_240', 7, 12), ('560_244_240', 8, 13), ('560_244_240', 9, 14),
         ('592_268_272', 2, 5), ('592_268_272', 3, 6), ('592_268_272', 4, 7), ('525_216_210', 1, 2), ('490_190_182', 1, 2)]
roman = {'560_244_240': 'I', '592_268_272': 'II', '525_216_210': 'IV', '490_190_182': 'V'}
for p, src, dst in DUALS:
    params, D = read_connset(os.path.join(CONN, f'srg1296_{p}_class{src}.txt'))
    S = dual_set(D, params[1])
    write_connset(os.path.join(CONN, f'srg1296_{p}_class{dst}.txt'), params, S,
                  f'graph {roman[p]}.{dst}: the dual set D* of graph {roman[p]}.{src} (family/duality.py)')
    print(f'{roman[p]}.{dst} = dual of {roman[p]}.{src}')
