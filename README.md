# Strongly regular graphs on 1296 vertices from partial difference sets in F₁₆ × F₈₁

Code and data for the paper *New strongly regular graphs on 1296 vertices from partial difference sets in
F₁₆ × F₈₁* (Dylan Plummer, draft, October 2026), in `paper/`.

## Results

| parameters | type | status before this work | graphs here (non-isomorphic) | \|Aut\| | ω |
|---|---|---|---|---|---|
| (1296, 560, 244, 240) | LS₁₆(36) | open ("OA(36,16)?") | 10 | 103680, 51840 (×3), 25920 (×2), 20736 (×4) | 12, 14 or 15 |
| (1296, 592, 268, 272) | NL₁₆(36) | open | 4 | 20736 each | 12 or 13 |
| (1296, 595, 274, 272) | LS₁₇(36) | realised independently by Kovaľ, by graphs that are not vertex-transitive | 4 | 20736 each | 24 |

Their complements (1296, 735, 414, 420), (1296, 703, 382, 380) and (1296, 700, 376, 380) are also realised.
Spectra: 560¹ 20⁵⁶⁰ (−16)⁷³⁵; 592¹ 16⁷⁰³ (−20)⁵⁹²; 595¹ 19⁵⁹⁵ (−17)⁷⁰⁰.

- All 18 graphs are vertex-transitive Cayley graphs, so their connection sets are abelian partial difference sets.
- None of the Latin-square-type graphs is the graph of an orthogonal array.
- Counts are lower bounds; the search is not exhaustive.
- Brouwer's table (checked 12 October 2026) lists (1296, 560, 244, 240) as realised, citing this work, and
  (1296, 592, 268, 272) as open.

## Construction

All graphs are Cayley graphs Cay(G, D) on G = (F₁₆, +) × (F₈₁, +) ≅ Z₂⁴ × Z₃⁴, with
F₁₆ = F₂[x]/(x⁴+x³+1), α = x, and F₈₁ = F₃[y]/(y⁴+y³+2), β = y (both primitive). In each case D = −D, 0 ∉ D, and
D is a union of orbits of a small group of field multiplications together with a Galois map:

- Graph I.1, srg(1296, 560, 244, 240): D is invariant under (a, b) ↦ (α³a, β⁴b) and (a, b) ↦ (a², b²⁷).
- Graphs II.1 and III.1, srg(1296, 592, 268, 272) and srg(1296, 595, 274, 272): D is invariant under the 4th roots
  of unity (a, b) ↦ (a, β²⁰b) and the Frobenius map (a, b) ↦ (a², b³).

The paper prints all three connection sets in closed form (Constructions 4.1–4.3).
Every connection set is in `connsets/*.txt` as discrete-log pairs: a line "i j" means (αⁱ, βʲ), and "-" stands for 0.

Structure (Section 5 of the paper):
- In graph I each coset {a} × F₈₁ induces the Brouwer–Haemers graph srg(81, 20, 1, 6).
- Graphs I.6 and I.10 arise from an srg(1296, 665, 340, 342) Cayley graph by deleting three pairwise orthogonal
  spreads of Delsarte cliques, a configuration asked for in Kovaľ's multi-orbit paper (Section 4.10).
- The graphs III arise from such a graph by deleting two orthogonal spreads.

## Layout

| path | content |
|---|---|
| `connsets/` | one connection set per isomorphism class (18 files: I.1–I.10, II.1–II.4, III.1–III.4) |
| `graphs/` | the 18 graphs in graph6 format, built by `build_from_logs.py` |
| `srg1296_*_classes.g6` | the same graphs collected per parameter set; line n is class n |
| `build_from_logs.py` | standalone builder and exhaustive verifier for any file in `connsets/` (plain Python 3, no packages) |
| `paper/` | LaTeX source (`srg1296.tex`, `refs.bib`), compiled PDF, figures |
| `paper/check_constructions.py` | rebuilds Constructions 4.1–4.3 exactly as printed and verifies them (plain Python 3) |
| `paper/make_figures.py` | regenerates `paper/figures/` from `connsets/` (numpy, matplotlib) |
| `paper/checks/` | computations behind the paper's claims: exact clique numbers (`clique_numbers.py`, output in `clique_numbers.out`), stabilisers and the known-construction checks (`known_constructions.py`), spread structure (`spread_structure.py`, `third_spread.py`) |
| `search/` | the search code, the target parameter sets, and a snapshot of all search logs; see `SEARCH.md` |
| `SEARCH.md` | how the search works, how it extends Kovaľ's, and its coverage so far (the search is still running) |
| `srg1296_visualization.html` | interactive visualisation of all 18 graphs (adjacency matrices, connection sets, ring layout) |

## Verification

```sh
python3 build_from_logs.py connsets/srg1296_592_268_272_class1.txt out.g6   # any connection set: VERIFIED
python3 paper/check_constructions.py                                       # Constructions 4.1-4.3 as printed
```

Both programs use plain Python 3 with no packages and no code shared with the search. They build Cay(G, D) with
their own arithmetic in F₁₆ and F₈₁, then check D = −D, 0 ∉ D, the degree, and the common-neighbour count for all
839,160 vertex pairs. All 18 graphs were also checked with Kovaľ's independent verifier (`verify_srg.py` in the
data of his paper *New strongly regular graphs on open parameter sets of Brouwer's table*,
doi:10.5281/zenodo.23210829). Isomorphism classes and automorphism groups were computed with nauty (pynauty).

The scripts in `paper/checks/` need numpy; `clique_numbers.py` also needs python-igraph.

## How the graphs were found

The graphs were found by an exact CP-SAT search for partial difference sets in products of finite fields that are
invariant under a subgroup of F_{q₁}* × F_{q₂}* plus a Galois map. It generalises Kovaľ's search on F₄ × F₈₁, which
produced srg(324, 136, 58, 56). See `SEARCH.md`.
