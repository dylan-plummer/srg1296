# The degree-24 fiber template is impossible

Let D be an inverse-closed connection set for SRG(1296,420,144,132) on
A x B = F2^4 x F3^4. Write T_a={v:(a,v) in D}. Let S consist of the
nonzero binary characters alpha for which (alpha,0)(D)=24. For a!=0,

    |T_a| = (9/4)(12 + W_S(a)).

Because T_a=-T_a and B has odd order, zero belongs to T_a exactly when
|T_a| is odd. Thus U={a:(a,0) in D} is determined by the parity of these
fifteen sizes. The graph has no loops, so 0 is absent from U.

The character-dual connection set D* consists of characters attaining 24.
It has size 420, is inverse-closed, and has the same SRG parameters:
Fourier inversion gives, for a nonprincipal character indexed by g,

    g(D*) = 36*1_D(g) - 12.

Consequently the dual's positive binary quotient set is U, whereas its
binary zero-fiber connection set is S. Applying the same size/parity rule
to U must therefore recover S. This is a necessary involution condition
on the binary quotient, not an assumption that the graph is self-dual.

If |T_0|=24, then S is an affine plane avoiding zero and has size four.
Its cross sizes are 36 once, 18 twice and 27 twelve times. Hence U is the
complement of a linear plane, of size twelve. The size/parity rule applied
to that U gives the complement of a linear plane again, of size twelve.
It cannot recover the original four-element S: contradiction.

All 105 labeled degree-24 quotient choices fail this test. The 35 degree-
42 choices survive: S is the complement of a linear plane, and the parity
map applied twice recovers it. Thus the degree-24 template is excluded
on this group, while the degree-42 template and the SRG row remain open.

The standard-library [checker](dual_parity_gate.py) enumerates all 140
quotients and tests the same rule on the verified degree-455,490,525
connection sets as positive controls. See [results](dual_parity_gate.json).
This supersedes the degree-24 recommendation in the preceding reassessment.
