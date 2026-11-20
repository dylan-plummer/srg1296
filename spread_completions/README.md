# Completions of the LS₁₆(36) Cayley graphs by spreads of cocliques

These graphs were found by GPT 6.2 Astra (OpenAI), 20 November 2026; see the paper, Section 6.5.

**Result.** There are exactly 25 isomorphism classes of srg(1296, 595, 274, 272) that contain a spanning copy of one of
the fourteen degree-560 Cayley graphs I.1–I.14 in `../connsets`. By symmetry type:
- 3 are Cayley graphs of Z₂⁴ × Z₃⁴ (uniform partitions, Cay(G, D ∪ U*));
- 1 is vertex-transitive and a Cayley graph only of non-abelian groups: SmallGroup(1296, 2797) and (1296, 2848);
- 1 is vertex-transitive and not a Cayley graph;
- 20 are not vertex-transitive.

None is the graph of an orthogonal array. All 25 are pairwise non-isomorphic, and none is isomorphic to III.1–III.4,
to Kovaľ's two Cayley graphs in `../graphs/koval`, or to Kovaľ's earlier multi-orbit graphs.

## Why spreads of cocliques

Let A be an LS_m(n) graph and P the projector onto its eigenvalue n − m, so that A = nP − mI + (m/n)J. Suppose B
is an LS_{m+d}(n) graph that contains every edge of A, with projector Q, and put D = B − A. Then:
- tr((Q − P)²) = d(n − 1) = rank Q − rank P;
- hence im P ⊆ im Q, and Q − P is a projector, so (D + dI)² = n(D + dI) + d(d − 1)J.

For d = 1, D is a disjoint union of n cliques K_n, and their vertex sets are cocliques of A of size n, which is the
Hoffman bound. Conversely, adding the edges inside a partition of the vertices into n-cocliques of A gives an
LS_{m+1}(n) graph.

So the spanning completions of A are exactly its partitions into 36-cocliques.

## The census

- **I.10.** Its 36-cocliques are the 288 cosets of eight subgroups. The exact cover problem on these cosets has
  1,124 solutions: 8 use one direction, 524 two and 592 three. Aut(I.10), of order 20,736, has 24 orbits on them, and
  the 24 graphs are pairwise non-isomorphic. Two of them are uniform (Cayley).
- **I.6.** It has five 36-cocliques through 0, all subgroups, meeting pairwise only in 0. So only uniform partitions
  exist, and they give one class.
- **The other twelve graphs** have no 36-coclique.

**The two exceptional vertex-transitive classes.**
- **class_18** has |Aut| = 5,184. Its regular subgroups are of types SmallGroup(1296, 2797) and (1296, 2848), both
  non-abelian. `pds_2797.json` and `pds_2848.json` give the connection sets as exponent tuples relative to
  Pcgs(SmallGroup(1296, ID)), and `verify_pds.g` checks all 354,025 differences in GAP.
- **class_24** has |Aut| = 2,592 and is vertex-transitive. Let g be the automorphism (a, y) ↦ (a ⊕ 1, Ly) with
  L² = −I. Then g² fixes 16 points, and every square lies in each index-2 subgroup. So no subgroup acts regularly, and
  the graph is not a Cayley graph. `square24.py` checks this on the file.

**Positive control.** The same construction applied to Kovaľ's graph B gives `mixed630_B.g6`, an
srg(1296, 630, 306, 306) with four vertex orbits. No census was made at degree 630.

## Files

| file | contents |
|---|---|
| `srg1296_595_spread_completions.zip` | the 25 graph6 files (`classes/class_01..24.g6`, `uniform6_1.g6`), construction data, the PDS exponent lists, and verifiers |
| `verify_bundle.py` | standard-library verifier. It checks all 20,979,000 vertex pairs, the four-clique totals, and a non-OA certificate per graph: an edge whose 274 common neighbors are colored with 33 colors |
| `rebuild_bundle.py` | rebuilds all 25 graph6 files byte-for-byte from the two degree-560 seeds and the partitions |
| `survey.py`, `cocliques.py`, `difference_orbits.py`, `spread_census.py`, `uniform6.py`, `audit.py` | coclique enumeration (two independent orders), the exact covers, the orbit census, and an OR-Tools cross-check |
| `compare.py`, `regular.py`, `regular.g`, `stabilizer24.py`, `square24.py`, `duals.py`, `edge_witnesses.py` | pair checks, canonical forms, regular subgroups (GAP), the non-Cayley witness, Fourier duals, non-OA witnesses |
| `build.py`, `mixed595.g6`, `mixed630_B.g6` | an explicit two-direction example, and the degree-630 control |
| `*.json`, `*.jsonl`, `*.out`, `*.log` | outputs of the runs above, including Kovaľ's verifier on all 25 graphs (`koval_verify.log`) |

To check the archive:

    python verify_bundle.py srg1296_595_spread_completions.zip
    python rebuild_bundle.py srg1296_595_spread_completions.zip
    gap -q -b verify_pds.g        # in the extracted archive
