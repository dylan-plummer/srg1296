# Claims and their evidence

Each claim of the paper (numbering of `paper/srg1296.tex`, draft of 16 November 2026) is listed with the program that
checks it and the stored output. `verify_all.sh` reruns everything.

Unless noted, programs are plain Python 3. Those marked † need numpy, igraph, pynauty or ortools
(`~/.conda/envs/srg`).

## The graphs (Theorems 1.1–1.3, 4.7; Table 6)

| claim | program | output |
|---|---|---|
| the 35 connection sets in `connsets/` are PDSs with the stated parameters (all 839,160 pairs) | `build_from_logs.py` | `family/verification_new_graphs.log` (new graphs); the 18 original graphs as in the October release |
| Constructions 4.1–4.6 as printed give these PDSs | `paper/check_constructions.py` | `paper/check_constructions.out` |
| independent check of all graph6 files | Kovaľ's `verify_srg.py` (not in this repository) | `family/verification_new_graphs.log` |
| \|Aut\|, pairwise non-isomorphism, the duality pairings of Proposition 2.5 | `family/nauty_check.py` † | `family/nauty_check.out` |
| the same with bliss; ring and column graphs; divisibility of Lemma 2.3 on every graph | `family/family_structure.py` † | `family/family_structure.out` |
| \|Stab(D)\| (Table 6) | `paper/checks/known_constructions.py` † | `paper/checks/known_constructions.out` |
| clique numbers (Section 6.1) | `paper/checks/clique_numbers.py` † | `paper/checks/clique_numbers.out` |
| free subgroups, spreads and the Bush-type matrices (Propositions 6.2, 6.3) | `paper/checks/spread_lattice.py` † | `paper/checks/spread_lattice.out`, `.json` |
| IV.1, V.1, VI.1, VII.1, VIII.1 equal the sets found in (Z/6)⁴; independent checks there | `family/import_c5.py`, `family/import_c3.py`; `fibre_join_as_run/k420/indep420.py`, `k385/indep385.py` | `fibre_join_as_run/*/` logs |
| VII.2, VII.3 = D(VIII.1) ∪ U*, D(VIII.2) ∪ U*; not isomorphic to VII.1 or each other | `family/spread_link.py` † | `family/spread_link.out` |
| VII.2, VII.3 are not OA graphs although ω = 36 (edges in no 36-clique) | `family/oa_test.py` † | `family/oa_test.out` |

## Structure (Section 2.3, Theorems 1.4, 1.5, Section 5)

| claim | program | output |
|---|---|---|
| Table 7 (skeleton), Table 8 (the nine local types L1–L9, with geometric identification of L1, L6, L8), zero-ring spectra, \|S_ψ\| vs \|U_b\|, cospectral zero rings of dual pairs | `family/shared_structure.py` | `family/shared_structure.out` |
| subgroups of order 36 inside / disjoint from every D (Section 5.4, Proposition 6.2) | `family/spread_inside.py` | `family/spread_inside.out` |
| Lemma 5.1 (Σc, Σc²); Lemma 5.2 (adding a free subgroup gives a PDS); Lemma 5.3 (spreads commute with duality) | proofs; consistent with `shared_structure.out`, `spread_link.out` | |
| Theorem 2.7 / Table 1: the binary fibre data of every LS_m(36) PDS, m ≤ 18 (enumeration of all 2¹⁵ sets, geometric identification) | `family/fibre_data.py` | `family/fibre_data.out` |
| Proposition 5.4: trace identities in the ρ-model; Gauss periods of index 8; NL_g(9) unions; ring sizes | `family/family_theory.py` | `family/family_theory.out` |
| Theorem 5.5, m = 5, 6: ρ-invariant PDSs exist | `family/rho_cpsat.py` †, `family/low_range/exhaustive_low.py` † | `family/low_range/*.txt` (verified by `build_from_logs.py`) |
| Theorem 5.5, m = 7: none exists (exhaustive) | `family/low_range/exhaustive_low.py 7 0,4 0` † | `family/low_range/exhaustive_m7.out` |
| … and by CP-SAT, two normalizations | `family/rho_cpsat.py 245 64 42 … 0,4 0 18` and `… 1,5 2 18` † | `family/low_range/solver_logs/rho_245*.log` (INFEASIBLE in 849 s and 125 s) |
| Section 5.6: LS_11 / LS_12 searches in the diag(F₄*, F₄*) models | `family/search_m11_m12/lsdiag.py` † | `family/search_m11_m12/*.jsonl` |

## Searches (Section 3)

| claim | where |
|---|---|
| CP-SAT field-model search and its coverage | `search/`, `SEARCH.md` |
| fibre-join constructions of IV.1, V.1, VI.1, VII.1, VIII.1, as run | `family/fibre_join_as_run/` |
| Section 7: σ-model attempts for LS_8, LS_9, LS_10 (CP-SAT UNKNOWN after 30 min each) | `family/sigma_cpsat.py` †; logs in `family/low_range/solver_logs/` |
| Section 3.3 / 5.5: CP-SAT timings quoted in the text | `family/low_range/solver_logs/` (see its README) |

## Not covered here

- The "Relation to known constructions" arguments of Section 6.4 are from the literature, apart from the cyclotomy
  check in `known_constructions.py`. The product-construction bullet relies on Tables 1.1 and 1.2 of Jedwab–Li,
  Combinatorial Theory 1 (2021) #18.
- The duality statement of Proposition 2.5 is cited (Delsarte), and is checked on every graph by difference counting.
