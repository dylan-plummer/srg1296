Two Cayley srg(1296,595,274,272), 20 Nov 2026

Files (one graph6 line each, n = 1296, terminated by a single newline):

  srg1296_595_274_272_cayley_A.g6   sha256 660d83d273ac78a744182378dcca015c0568708fea84896acc3dca890dbc3a13
  srg1296_595_274_272_cayley_B.g6   sha256 17928169a8fa034b9b6b2433639e04733e76f035dde69bc19c35f729ec523d1e

Parameters: (v,k,lambda,mu) = (1296,595,274,272) for both graphs (not complements; degree 595).

Checks, re-run 20 Nov 2026 with our own graph6 decoder:
  - simple, 595-regular, exact integer identity A^2 = 595 I + 274 A + 272 (J - I - A);
  - separately, all 839,160 unordered vertex pairs by neighbour-set intersection (274 for adjacent, 272 for
    non-adjacent pairs): no exceptions;
  - |Aut| = 20736 for each (nauty), vertex-transitive (one orbit);
  - A and B are not isomorphic to each other or to any of the four classes srg1296_595_274_272_class1..4 in
    github.com/dylan-plummer/srg1296 (commit 8b151dd): six pairwise different nauty-labelg canonical forms.

Group and labelling. Both are Cayley graphs on G = Z2^4 x Z3^4 = (F16,+) x (F81,+). Vertex i (0-based) is the
element (a, x) with a = floor(i/81), x = i mod 81; a is read as 4 bits (coefficients of X^0..X^3 in
F16 = F2[X]/(X^4+X^3+1)), x as 4 base-3 digits (coefficients of Y^0..Y^3 in F81 = F3[Y]/(Y^4+Y^3+2)), lowest
coefficient = lowest digit. In this labelling all translations of G are automorphisms (checked for the 8
generators), so G acts regularly, and the connection set is N(0) (closed under negation). Slice sizes
|N(0) ∩ ({a} x F81)| = 28 for a = 0, 45 for three nonzero a, 36 for the other twelve.

How found (20 Oct 2026): a CP-SAT search over connection sets that are unions of orbits of a group K of semilinear
maps of F16 x F81 (your K(S, sigma) restriction of the rational PDS model). For both graphs N(0) is invariant under
multiplication by the elements of order 4 of F81* (in the second coordinate) and under a field automorphism:
  A: (x, y) -> (x^2, y^3)
  B: (x, y) -> (x^2, y^27)
(re-checked 20 Nov with our own field arithmetic). Each solution was rebuilt as an explicit Cayley graph and verified
exactly before being kept.
