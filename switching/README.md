# Non-vertex-transitive graphs by switching binary fibers

These graphs were found by GPT 6.2 Astra (OpenAI) after an audit of the paper. See Section 6.4 of the paper.

We give 68 pairwise non-isomorphic strongly regular graphs on 1296 vertices for four parameter sets that this
repository already realizes by Cayley graphs. Each of the 68 has at least two different numbers of 4-cliques through a
vertex. So none is vertex-transitive, and none is isomorphic to a Cayley graph, in particular to any graph in `graphs/`.

| parameters | new classes certified here |
|---|---|
| (1296, 385, 124, 110) | 1 |
| (1296, 420, 144, 132) | 1 |
| (1296, 560, 244, 240) | 33 |
| (1296, 592, 268, 272) | 33 |

## Sources and coordinates

The input is the class 1 connection set of the given degree in `../connsets`. The coordinates are those of
`../family/fields.py`. A vertex is a pair (x, y) in F₂⁴ × F₃⁴:
- the four binary coordinates of x are read as a little-endian integer x = 0..15;
- the ternary 4-tuples y are ordered lexicographically, numbered from 0;
- vertex 16i + x has binary part x and ternary tuple number i.

## The switch

Let H_ab = (−1)^{a·b} be the unnormalized binary Walsh matrix of order 16. Choose a frequency set M not containing 0,
and put e_a = −1 for a in M and e_a = 1 otherwise. Then Q = H diag(e) H / 16 is a symmetric orthogonal involution that
fixes the all-ones vector. For M = all nonzero frequencies it is the Godsil–McKay matrix (2/16) J − I.

Choose a set S of the 81 fibers F₁₆ × {y}. Let U be block diagonal, with Q on the fibers in S and I elsewhere. The new
adjacency matrix is U A U. Every block of A is invariant under binary translations, so it commutes with Q. Hence only
the blocks that cross the cut can change.

For a ternary difference d, let b_d(x) = 1 if (x, d) is in the connection set. A crossing block has the new kernel

    b'_d = ((b_d H) * e) H / 16      (* = coordinatewise product).

The exact test is that every numerator is 0 or 16. If it passes:
- the diagonal blocks are unchanged, so the new diagonal is zero;
- symmetry is automatic;
- conjugation by U fixes I and J, so it preserves A² = (k − μ) I + (λ − μ) A + μ J and the degree.

### Hyperplane examples

The four examples below use S = {y : y₄ = 0} (one-based coordinates). The last two columns give the number of
4-cliques through vertex 0 (on the hyperplane) and through vertex 16 (off it).

| degree | M | \|Aut\| | vertex 0 | vertex 16 |
|---|---|---|---|---|
| 385 | 5, 9 | 864 | 338402 | 338646 |
| 420 | 3, 11 | 864 | 529872 | 530616 |
| 560 | 6, 7, 10, 12, 14 | 1728 | 2431560 | 2433372 |
| 592 | 3, 6, 10, 12, 13, 14 | 864 | 3167996 | 3166426 |

In every row the two clique counts occur at 432 and 864 vertices respectively. bliss independently finds
automorphism orbits of these sizes.

### Arbitrary subsets of the 81 fibers

Take M = {6, 7, 9, 10, 12, 14} for degree 560 and M = {3, 6, 10, 12, 13, 14} for degree 592. Then all 80 nonzero
ternary differences pass the test, so every subset S is allowed.

Exactly 52 differences change, and they generate F₃⁴. So the graph on the 81 fibers that links pairs at these
differences is connected. Two choices of S give the same labelled adjacency matrix exactly when they differ by a
constant. Hence there are exactly 2⁸⁰ distinct labelled graphs for each of the two parameter sets.

`sample_560.jsonl` and `sample_592.jsonl` give 32 cut masks per degree. Each record contains:
- the full local 4-clique histogram;
- the SHA-256 of the graph6 file;
- an independent census of all vertex pairs.

The 32 histograms differ from each other and from the hyperplane histogram. This certifies at least 33 isomorphism
classes per degree. The sample with index 0 in each row has full automorphism group of order 16: the binary
translations, acting with the 81 fibers as orbits. We make no claim about how many of the other labelled graphs are
isomorphic.

## Files

| file | contents |
|---|---|
| `non_cayley_68.zip` | the 68 graph6 files and a manifest with checksums |
| `construct.py`, `family.py`, `scan.py` | the construction (numpy, python-igraph; srg environment) |
| `verify_graph6.py` | standalone verifier, plain Python. It decodes the graph6 bits, checks the degree and every unordered pair, checks the binary translations, and counts triangles in one neighborhood per fiber |
| `verify_switch.py` | independent check. It rebuilds every row from the source graphs in `../graphs` and the manifest by direct character sums, and checks the Boolean test and the 52 changed differences |
| `verify_archive.py` | runs A. Kovaľ's `verify_srg.py` (not included; set `KOVAL_VERIFY_SRG` to its path) on the archive; its 68 reports are in `koval_verification.jsonl` |
| `census.py`, `export_collection.py` | reconcile the certificates and graph hashes, and build the archive |
| `*.jsonl`, `*.json` | outputs of the runs above |

To reproduce:

    python scan.py > scan.jsonl
    python construct.py 385 1        # likewise 420, 560, 592
    python family.py 560 32 > sample_560.jsonl
    python family.py 592 32 > sample_592.jsonl
    python verify_graph6.py <graph6 files>
    python census.py
    python verify_switch.py
