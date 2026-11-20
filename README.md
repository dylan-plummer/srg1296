# Strongly regular graphs on 1296 vertices for every Latin square parameter set

Code and data for the paper *Strongly regular graphs on 1296 vertices for every Latin square parameter set*
(Dylan Plummer, draft, November 2026), in `paper/`.

## Results

**Every Latin square parameter set LS_m(36), 1 ≤ m ≤ 35, is now realized by a strongly regular graph.** The open cases
were m = 11, …, 17 and their complements. This repository realizes m = 11, …, 16 (and vertex-transitive m = 17),
all as Cayley graphs of Z₂⁴ × Z₃⁴. The rest: m ≤ 10 from orthogonal arrays; m = 17 from Kovaľ and from our graphs III;
m = 18, 19 from graph I.6 plus two or three free subgroups (Lemma 5.2); m ≥ 20 by complementation.

| parameters | type | status before this work | graphs here (non-isomorphic) | \|Aut\| | ω |
|---|---|---|---|---|---|
| (1296, 385, 124, 110) | LS₁₁(36) | open ("OA(36,11)?") | 2 (VIII.1, VIII.2) | 7776 each | 18, 13 |
| (1296, 420, 144, 132) | LS₁₂(36) | open ("OA(36,12)?") | 3 (VII.1–VII.3) | 7776 each | 18, 36, 36 |
| (1296, 455, 166, 156) | LS₁₃(36) | open ("OA(36,13)?") | 1 (VI.1) | 12960 | 24 |
| (1296, 490, 190, 182) | LS₁₄(36) | open ("OA(36,14)?") | 2 (V.1, V.2) | 12960 each | 18, 12 |
| (1296, 525, 216, 210) | LS₁₅(36) | open ("OA(36,15)?") | 2 (IV.1, IV.2) | 12960 each | 12 |
| (1296, 560, 244, 240) | LS₁₆(36) | open ("OA(36,16)?") | 14 (I.1–I.14) | 103680, 51840 (×3), 25920 (×3), 20736 (×7) | 12–15 |
| (1296, 592, 268, 272) | NL₁₆(36) | open | 7 (II.1–II.7) | 20736 each | 12 or 13 |
| (1296, 595, 274, 272) | LS₁₇(36) | realized independently by Kovaľ, by graphs that are not vertex-transitive | 4 (III.1–III.4) | 20736 each | 24 |

Their complements are also realized.

- All 35 graphs are vertex-transitive Cayley graphs, so their connection sets are abelian partial difference sets.
- None of the Latin-square-type graphs is the graph of an orthogonal array. VII.2 and VII.3 contain one spread of
  Delsarte cliques, but some of their edges lie in no 36-clique.
- Counts are lower bounds; the searches are not exhaustive.
- **Two further Cayley graphs with parameters (1296, 595, 274, 272)**, found by A. Kovaľ on 20 October 2026 in the
  same two models as III.1–III.4 and provided by him, are in [`graphs/koval/`](graphs/koval/) (his README included).
  With III.1–III.4 they form six pairwise non-isomorphic vertex-transitive graphs (paper, Remark 4.9).
- **68 graphs that are not vertex-transitive** (1 each at k = 385, 420 and 33 each at k = 560, 592) come from a
  switching of the binary fibers of our Cayley graphs; see [`switching/`](switching/) and Section 6.4 of the paper.
  They were found by GPT 6.2 Astra (OpenAI) after an audit of the paper.
- **25 spread completions with parameters (1296, 595, 274, 272)**: exactly 25 classes contain a spanning copy of one of
  I.1–I.14 (Lemma 6.4: spanning LS₁₇ ⊃ LS₁₆ ⇔ partition into 36-cocliques). Three are abelian Cayley graphs, one is a
  Cayley graph only of non-abelian groups, one is vertex-transitive but not Cayley, and 20 are not vertex-transitive;
  see [`spread_completions/`](spread_completions/) and Section 6.5 of the paper.
- Status of the rows: Brouwer's table (checked 14 October 2026) lists (1296, 560, 244, 240) and (1296, 592, 268, 272)
  as realized, citing this work. The rows with k = 525, 490, 455, 420, 385 were open when this work began.

**Shared structure** (Section 5 of the paper, `family/`):
- a forced binary skeleton (planes, ovoids, grids in PG(3,2), by m mod 8);
- a local alphabet of nine ring-set types for all 35 graphs;
- two symmetry types: order 5 for m = 13–16, order 3 for m = 11, 12;
- links by duality and by spreads, e.g. LS₁₁ + one spread = LS₁₂.

The **structure theorem** (Theorem 2.7 of the paper, `family/fibre_data.py`) covers every LS_m(36) PDS in
Z₂⁴ × Z₃⁴ with m ≤ 18. In any such PDS, the elements in F₁₆ × {0} and the positive characters trivial on F₈₁ form
ovoids, grids, planes or their complements in PG(3,2), according to m mod 8; |T₀| and all ring sizes are forced.

The graphs I, IV, V, VI form the **order-five family** (Section 5): each F₈₁-coset induces the cyclotomic NL_g(9) graph,
with g = 18 − m. Under ρ = (α¹², β¹⁶), LS_m(36) occurs exactly for m = 5, 6, 13, 14, 15, 16; m = 7 is excluded
exhaustively.

## Construction

All graphs are Cayley graphs Cay(G, D) on G = (F₁₆, +) × (F₈₁, +) ≅ Z₂⁴ × Z₃⁴, with
F₁₆ = F₂[x]/(x⁴+x³+1), α = x, and F₈₁ = F₃[y]/(y⁴+y³+2), β = y (both primitive). In each case D = −D and 0 ∉ D.
For graphs I–VI, D is a union of orbits of a small group of field multiplications, sometimes together with a Galois map:

- Graph I.1, srg(1296, 560, 244, 240): D is invariant under (a, b) ↦ (α³a, β⁴b) and (a, b) ↦ (a², b²⁷).
- Graphs II.1 and III.1, srg(1296, 592, 268, 272) and srg(1296, 595, 274, 272): D is invariant under the 4th roots
  of unity (a, b) ↦ (a, β²⁰b) and the Frobenius map (a, b) ↦ (a², b³).

- Graphs IV.1, V.1, VI.1 (k = 525, 490, 455): D is invariant under ρ = (α¹², β¹⁶), of order 5, and −1 (Construction 4.4).
- Graphs VII.1 (k = 420) and VIII.1 (k = 385): D is invariant under (α⁵a, u(b)), u unipotent of order 3 (not a field multiplication) (Constructions 4.5, 4.6). VII.2, VII.3 are D(VIII.1) ∪ U* and D(VIII.2) ∪ U* for a free subgroup U of order 36.
- The duals I.11–I.14, II.5–II.7, IV.2, V.2 and VIII.2 are the sets D* of Proposition 2.5.

The paper prints the connection sets of I.1, II.1, III.1, IV.1, V.1, VI.1, VII.1 and VIII.1 (Constructions 4.1–4.6).
Every connection set is in `connsets/*.txt` as discrete-log pairs: a line "i j" means (αⁱ, βʲ), and "-" stands for 0.

Structure (Section 6 of the paper; the fiber structure and the order-five family are in Sections 2.3 and 5 and in `family/`):
- In graph I each coset {a} × F₈₁ induces the Brouwer–Haemers graph srg(81, 20, 1, 6).
- Graphs I.6 and I.10 arise from an srg(1296, 665, 340, 342) Cayley graph by deleting three pairwise orthogonal
  spreads of Delsarte cliques. That Cayley graph is the graph of a symmetric Bush-type Hadamard matrix, so this answers
  the question of Kovaľ's multi-orbit paper (Section 4.10) positively for n = 36 (Proposition 6.3).
- The graphs III arise from such a graph by deleting two orthogonal spreads (graphs Γ_K of Bush-type Hadamard matrices).
- Adding pairwise orthogonal subgroups of order 36 that are disjoint from D always gives a Latin square type PDS:
  graph I.6 is an LS₂₁(36) Cayley graph with five orthogonal spreads deleted (the maximum in this group), I.10 an
  LS₂₀(36) graph with four, and the graphs III LS₂₁(36) graphs with four (Proposition 6.2, Figure 4).

## Layout

| path | content |
|---|---|
| `connsets/` | one connection set per isomorphism class (35 files: I.1–I.14, II.1–II.7, III.1–III.4, IV.1–IV.2, V.1–V.2, VI.1, VII.1–VII.3, VIII.1–VIII.2) |
| `graphs/` | the 35 graphs in graph6 format, built by `build_from_logs.py` |
| `srg1296_*_classes.g6` | the same graphs collected per parameter set; line n is class n |
| `build_from_logs.py` | standalone builder and exhaustive verifier for any file in `connsets/` (plain Python 3, no packages) |
| `paper/` | LaTeX source (`srg1296.tex`, `refs.bib`), compiled PDF, figures |
| `paper/check_constructions.py` | rebuilds Constructions 4.1–4.6 exactly as printed and verifies them (plain Python 3; output in `check_constructions.out`) |
| `paper/make_figures.py` | regenerates `paper/figures/rings.pdf` (Figure 3) from `connsets/` (numpy, matplotlib) |
| `paper/make_fig_orthogonality.py` | regenerates `paper/figures/orthogonality.pdf` (Figure 4(b)) from `paper/checks/spread_lattice.json` |
| `paper/make_fig_family.py` | regenerates `paper/figures/zero_rings.pdf` (Figure 1: the zero ring for m = 11, …, 17) from `connsets/` (numpy, matplotlib) |
| `paper/checks/` | computations behind the paper's claims: exact clique numbers (`clique_numbers.py`, output in `clique_numbers.out`), stabilizers and the known-construction checks (`known_constructions.py`), spread structure (`spread_structure.py`, `third_spread.py`), and the independent check of Propositions 5.2 and 5.3 (`spread_lattice.py`: all subgroups of order 36, every orthogonal set via the character criterion, and the Bush-type Hadamard matrices; output in `spread_lattice.out` and `.json`) |
| `family/` | Sections 2.3–5 of the paper: the structure theorem, the order-five family, duals, imports of graphs IV–VII, the fiber-join code as run, the LS₁₁/LS₁₂ searches; see `family/README.md` |
| `CHECKS.md`, `verify_all.sh` | every claim of the paper with the program and output that support it; a script that reruns them all |
| `search/` | the search code, the target parameter sets, and a snapshot of all search logs; see `SEARCH.md` |
| `SEARCH.md` | how the search works, how it extends Kovaľ's, and its coverage |
| `srg1296_visualization.html` | interactive visualization of the 18 graphs I.1–I.10, II.1–II.4, III.1–III.4 (adjacency matrices, connection sets, ring layout) |

## Verification

```sh
python3 build_from_logs.py connsets/srg1296_592_268_272_class1.txt out.g6   # any connection set: VERIFIED
python3 paper/check_constructions.py                                       # Constructions 4.1-4.6 as printed
./verify_all.sh                                                            # everything in CHECKS.md
```

Both programs use plain Python 3 with no packages and no code shared with the search. They build Cay(G, D) with
their own arithmetic in F₁₆ and F₈₁, then check D = −D, 0 ∉ D, the degree, and the common-neighbor count for all
839,160 vertex pairs. All 35 graphs were also checked with Kovaľ's independent verifier (`verify_srg.py` in the
data of his paper *New strongly regular graphs on open parameter sets of Brouwer's table*,
doi:10.5281/zenodo.23210829). Isomorphism classes and automorphism groups were computed with nauty (pynauty) and
bliss (igraph). Graphs IV.1–VIII.1 were also verified in their original (Z/6)⁴ coordinates (`family/fibre_join_as_run/`).

The scripts in `paper/checks/` need numpy; `clique_numbers.py` also needs python-igraph.

## How the graphs were found

Graphs IV–VIII were found by the fiber decomposition of Section 3.4 (`family/fibre_join_as_run/`); the duals by
Proposition 2.5. The graphs I–III were found by an exact CP-SAT search for partial difference sets in products of finite fields that are
invariant under a subgroup of F_{q₁}* × F_{q₂}* plus a Galois map. It generalizes Kovaľ's search on F₄ × F₈₁, which
produced srg(324, 136, 58, 56). See `SEARCH.md`.
