# Verified construction: SRG(1296,525,216,210)

16 November 2026. **FOUND AND INDEPENDENTLY VERIFIED.** The requested graph
exists. All deliverables and an exact three-seed construction are in
[result/README.md](result/README.md). The complement SRG(1296,770,454,462)
was also fully verified. Brouwer's live table still marks the degree-525 row
open on the date of this result.

## What made the construction possible

The fixed-matching geometric local family had an exact positive-character
energy lower bound of 840 against the available 800. Increasing a solver
budget inside that family could not succeed. We enlarged the local domains
to all 13-/18-point projective moment-code cross-sets satisfying the retained
matching conditions, allowing additional matching channels.

`local_domains.cpp` enumerated the local domains completely using a 20+20
meet-in-the-middle join: 264 size-27 sets for each of five fixed matching
choices, and 580615 admissible size-36 sets. Local enumeration took 1.69675
seconds. `unsigned_domains.py` and `signed_domains.py` imposed the shared
Sylvester frame and exact Fourier projection conditions. The surviving
absolute local spectra were injective, removing a separate ambiguity in
the signs. Binary difference consistency is essential: the earlier scalar
13-old/3-new orbit aggregate in the planning note was not a valid global
assignment in this frame.

The final three domains had sizes 336665, 336665, 252.
`join_domains.py` found choices [292431,160395,5] after 2910 tested pairs in
3.2212539929896593 seconds, including indexing but excluding preceding domain
precomputation and final matrix checks. The three selected local sets have
projective masks 0x931645b94a, 0x559742520f and 0x381ff1, with sizes 36, 36
and 27 respectively. Their explicit points and simultaneous order-five
development appear in the standalone certificate.

## Verification and deliverables

The search reconstruction checked all 1296 degrees and all 839160 unordered
vertex pairs. The separately written standard-library verifier read the saved
matrix and checked those conditions again, checked the full partial difference
set identity directly in (Z/6Z)^4, and matched the resulting connection set
to every matrix row. It also verified the complement and both graph6 files.
NetworkX independently decoded the primary graph6 and matched every adjacency.
A standalone rebuild in a clean output directory reproduced every data-file
hash. None of these checks depends on accepting a solver status or the
correctness of the search-space reduction.

Primary adjacency SHA256:
`40fa23bf70f8a020bc90539be56cf4c490f480cb7b83a4c5a05c2220bb980bb9`.

The builder and independent verifier need only standard Python, and take
approximately one and three seconds here. Their certificate has three seeds;
there is no need to repeat the search to reconstruct the graph.

## Scope and next work

The requested parameter row is constructed. Stop treating it as a pending
construction attempt. This does not construct OA(36,15), classify all such
graphs, or compute the full automorphism group. The initial CP-SAT UNKNOWN
pilot and infeasible additional self-duality ansatz remain records of those
specific earlier tests; neither conflicts with this successful construction.

The useful next scientific work is to explain the seed geometry and determine
which parts of the construction vary or transfer to other parameters. Public
table updates, messages to collaborators, publication, and broader searches
have not been undertaken as part of this deliverable.
