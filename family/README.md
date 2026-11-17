# Fibre structure, shared structure, and graphs IV–VIII (Sections 2.3, 3.4, 4, 5 of the paper)

Notation follows the paper:
- ρ = (α¹², β¹⁶) is multiplication by elements of order 5 in both fields.
- T_a = {b : (a, b) ∈ D} are the ring sets and U_0 = {a : (a, 0) ∈ D}.
- S₁ is the set of characters φ of F₁₆ with (φ × 1)(D) = r.

See `../CHECKS.md` for the map from claims to programs.

| file | purpose | needs |
|---|---|---|
| `fields.py` | field arithmetic, connset I/O (same conventions as `../build_from_logs.py`) | — |
| `duality.py` | exact character values via the trace forms; the dual set D* of Proposition 2.5 | — |
| `fibre_data.py` → `fibre_data.out` | Theorem 2.7 and Table 1: forced S₁, U₀, \|T₀\|, ring sizes for every LS_m(36), m ≤ 18, with geometric identification (planes, ovoids, grids) | — |
| `family_theory.py` → `family_theory.out` | finite checks in the proof of Proposition 5.4 (ρ-model): trace identities, Gauss periods of index 8, NL_g(9) unions, Walsh sums of μ₅-cosets | — |
| `import_c5.py`, `import_c3.py` | map the (Z/6)⁴ connection sets of the fibre-join runs to F₁₆ × F₈₁ (IV–VI with R ↦ ρ; VII, VIII with their binary C3 ↦ α⁵) | — |
| `shared_structure.py` → `shared_structure.out` | Section 5: local types L1–L9 of all ring sets (Table 8) with geometric identification, zero-ring spectra, \|U_b\| and \|S_ψ\| distributions, matching | — |
| `spread_inside.py` → `spread_inside.out` | subgroups of order 36 inside / disjoint from every connection set | — |
| `spread_link.py` → `spread_link.out` | VII.2, VII.3 = D(VIII.i) ∪ U*: PDS check and non-isomorphism | numpy, pynauty |
| `oa_test.py` → `oa_test.out` | VII.2, VII.3 have ω = 36 but edges in no 36-clique (not OA graphs) | numpy, igraph |
| `sigma_cpsat.py` | CP-SAT in the model ⟨σ⟩ with forced data fixed (Section 7: LS₈–LS₁₀, all UNKNOWN after 30 min) | ortools |
| `canon.py` | closed form D(J₀, I₀, S) of Construction 4.4, canonical under the normalizer of ⟨ρ⟩ | — |
| `make_connsets.py` | writes IV.1, V.1, VI.1 and the nine new duals into `../connsets` | — |
| `family_structure.py` → `family_structure.out` | all 35 graphs: ring and column graphs, ring sizes, Lemma 2.3 divisibility, \|S_ψ\|, \|Aut\| (bliss), self-duality | numpy, igraph |
| `nauty_check.py` → `nauty_check.out` | \|Aut\|, pairwise non-isomorphism and duality pairings with nauty | numpy, pynauty |
| `rho_cpsat.py` | CP-SAT in the model ⟨ρ⟩, optionally with T₀, U₀ and ring sizes fixed as in Proposition 5.4 | ortools |
| `low_range/` | Theorem 5.5 for m = 5, 6, 7: `exhaustive_low.py` (solver-free; outputs `exhaustive_m*.out`) and the PDSs found for m = 5, 6 | numpy |
| `search_m11_m12/` | Section 5.6: diag(F₄*, F₄*) models for LS₁₂ and LS₁₁, with logs | numpy |
| `fibre_join_as_run/` | the fibre-join code, logs and certificates that found IV.1 (`k525/`), V.1 and VI.1 (`k490_k455/`), VII.1 (`k420/`) and VIII.1 (`k385/`), copied as run | see below |
| `verification_new_graphs.log` | `../build_from_logs.py` and Kovaľ's `verify_srg.py` on all new graphs | — |

## Results in one table

| graph | parameters | T₀ | U₀ (type) | ring sizes | dual |
|---|---|---|---|---|---|
| I.1–I.14 | LS₁₆(36) | Brouwer–Haemers NL₂(9) | ∅ | 36 ×15 | I.5↔I.11, I.7↔I.12, I.8↔I.13, I.9↔I.14; others self-dual |
| IV.1, IV.2 | LS₁₅(36) | NL₃(9), 3 index-8 classes | ovoid (Clebsch) | 27 ×5, 36 ×10 | IV.1 ↔ IV.2 |
| V.1, V.2 | LS₁₄(36) | NL₄(9), 4 classes | complement of ovoid | 27 ×10, 36 ×5 | V.1 ↔ V.2 |
| VI.1 | LS₁₃(36) | NL₅(9), 5 classes | all of F₁₆* | 27 ×15 | self-dual |
| VIII.1, VIII.2 | LS₁₁(36) | 34-set, not a PDS (−11, −2, 7, 16) | grid | 18 ×6, 27 ×9 | VIII.1 ↔ VIII.2 |
| VII.1–VII.3 | LS₁₂(36) | 42-set, not a PDS (−3, 6, 15) | complement of a plane | 18 ×3, 27 ×12 | VII.1 self-dual; VII.2 ↔ VII.3 |
| III.1–III.4 | LS₁₇(36) | 28-set, not a PDS | plane minus 0 (F₄*) | 36 ×12, 45 ×3 | self-dual |

The U₀ column and the ring sizes are those forced by Theorem 2.7 for every Cayley graph on Z₂⁴ × Z₃⁴ with these parameters.

## Reproduce

```sh
python3 fibre_data.py > fibre_data.out
python3 family_theory.py > family_theory.out
python3 family_structure.py > family_structure.out       # about 2 minutes
python3 nauty_check.py > nauty_check.out                 # about 1 minute
python3 low_range/exhaustive_low.py 7 0,4 0              # about 15 minutes: RESULT m=7: NO solution (exhaustive)
```

The imports were:

```sh
python3 import_c5.py RAW_Z6_CONNECTIONS.txt 525 216 210 raw525.txt   # likewise 490, 455
python3 make_connsets.py raw525.txt raw490.txt raw455.txt
python3 import_c3.py RAW_Z6_CONNECTIONS_420.txt ../connsets/srg1296_420_144_132_class1.txt
```

The raw (Z/6)⁴ connection sets are `connections_Z6_4.txt` in the result folders of the fibre-join runs.

## Notes on `fibre_join_as_run/`

This folder is a record, not a package. The scripts contain the absolute paths of the working tree where they ran
(`work/construct525_20261116`, `work/transfer525_20261116`, `work/construct420_20261116`, `work/transfer420_20261116`). The large intermediate files
(local catalogues, `.npz` domains) are not included. Each can be regenerated in seconds by the `.cpp` and `.py` files there.

The certificates (`construction.json`, `FOUND455_certificate.json`, `FOUND490_certificate.json`, `k420/construction.json`)
give the seeds explicitly. The joins used optional restrictions on the positive Walsh patterns, which are not necessary
conditions, and they stopped at the first solution.

`k420/dual_parity_gate.py` is the first version of Lemma 2.6. It was applied to LS₁₂ and checked on IV–VI; Theorem 2.7
supersedes it.

## Notes on `search_m11_m12/`

```sh
g++ -O2 -o local_catalogue local_catalogue.cpp && ./local_catalogue cat_input.txt cat   # cat_0.tsv: all 9-divisible 13-/18-point sets
python3 pairs.py && python3 dedupe.py                                                 # pairs_cls.txt: 89 classes of (g1, g2)
python3 lsdiag.py pairs_cls.txt 0 89                                                  # LS_12, data (42,12,12,42)
M_LS=11 T0PTS=17 E1=9,0 E2=9,0 python3 lsdiag.py pairs_cls.txt 0 89                   # LS_11, data (34,9,9,34)
```

`MAXTRIP` (default 2·10⁵) skips classes with more triples. The logs record every class:
- **LS₁₂:** all 89 classes ran, with 0 solutions. VII.1–VII.3 are therefore not invariant under such a group.
- **LS₁₁:** 86 classes ran, with 0 solutions.
  - Two classes with 1.57·10⁶ triples were rerun with `MAXTRIP=2000000`; their lines appear twice in the log, first as skipped, then completed.
  - Three classes with about 9.9·10⁷ triples each were not run.

By Theorem 2.7, taking S₁ to be the 3×3 grid (the orbit of size 9) is no restriction for LS₁₁; only the symmetry is. VIII.1 was later found in the smaller model ⟨σ⟩.

`test_engine.py` checks the engine's row test on IV.1, V.1, VI.1 and the dual-action plumbing. It reads the (Z/6)⁴
connection sets from the working tree.
