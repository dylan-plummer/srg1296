#!/usr/bin/env python3
"""Figure: the zero ring R_0 = {0} x F_81 of one graph for each m = 11,...,17 (Section 5), i.e. Cay(F_81, T_0).
(0,0) in the center, its neighbors in R_0 on the inner circle and its non-neighbors on the outer circle (at angle
2 pi j/80 for beta^j); vertices colored by their number of common neighbors with (0,0) inside R_0.
usage: python3 make_fig_family.py   (writes figures/zero_rings.pdf; needs numpy, matplotlib)"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
HERE = os.path.dirname(os.path.abspath(__file__)); CONN = os.path.join(HERE, "..", "connsets"); OUT = os.path.join(HERE, "figures")
plt.rcParams.update({"font.family": "serif", "font.size": 7, "pdf.fonttype": 42, "mathtext.fontset": "cm"})
def powers(p, poly, n):
    cur = [1] + [0] * (n - 1); out = []
    for _ in range(p ** n - 1):
        out.append(tuple(cur)); top = cur[-1]; cur = [0] + cur[:-1]
        cur = [(c - top * poly[i]) % p for i, c in enumerate(cur)]
    return out
B81 = powers(3, [2, 0, 0, 1, 1], 4); Z = (0, 0, 0, 0)
F3 = [Z] + B81; idx = {b: i for i, b in enumerate(F3)}; logB = {b: j for j, b in enumerate(B81)}
sub = np.array([[idx[tuple((p - q) % 3 for p, q in zip(x, y))] for y in F3] for x in F3])
PANELS = [("385_124_110", 1, "VIII.1", 11, "grid"), ("420_144_132", 1, "VII.1", 12, "complement of a plane"),
          ("455_166_156", 1, "VI.1", 13, "all of $\\mathbb{F}_{16}^*$"), ("490_190_182", 1, "V.1", 14, "complement of an ovoid"),
          ("525_216_210", 1, "IV.1", 15, "ovoid"), ("560_244_240", 1, "I.1", 16, "empty"), ("595_274_272", 1, "III.1", 17, "plane minus $0$")]
BLUE, ORANGE, AQUA, INK, GREY = "#2a78d6", "#eb6834", "#1baf7a", "#1b1f27", "#9aa1ad"
fig, axes = plt.subplots(2, 4, figsize=(6.9, 3.95), gridspec_kw={"wspace": 0.04, "hspace": 0.32})
axes = axes.ravel(); axes[-1].axis("off")
axes[-1].text(0.02, 0.55, "center: $(0,0)$\ninner circle: its neighbors in $R_0$\nouter circle: its non-neighbors in $R_0$\ncolor: number of common\nneighbors with $(0,0)$ in $R_0$", transform=axes[-1].transAxes, fontsize=6.5, va="center", linespacing=1.4)
for ax, (p, c, name, m, u0) in zip(axes, PANELS):
    T0 = set()
    for l in open(os.path.join(CONN, f"srg1296_{p}_class{c}.txt")):
        if l.startswith("#") or len(l.split()) != 2: continue
        i, j = l.split()
        if i == "-" and j != "-": T0.add(B81[int(j) % 80])
    inT = np.zeros(81, bool); inT[[idx[b] for b in T0]] = True
    R = inT[sub]; np.fill_diagonal(R, False); C = R.astype(int) @ R.astype(int)
    inn = [x for x in range(1, 81) if R[0, x]]; out = [x for x in range(1, 81) if not R[0, x]]
    lamv = sorted({int(C[0, x]) for x in inn}); muv = sorted({int(C[0, x]) for x in out})
    srg = len({int(C[x, y]) for x in range(81) for y in range(81) if R[x, y]}) == 1 and \
          len({int(C[x, y]) for x in range(81) for y in range(81) if x != y and not R[x, y]}) == 1
    zp = {0: np.zeros(2)}
    for x in inn: th = 2 * np.pi * logB[F3[x]] / 80; zp[x] = 1.15 * np.array([np.cos(th), np.sin(th)])
    for x in out: th = 2 * np.pi * logB[F3[x]] / 80; zp[x] = 1.75 * np.array([np.cos(th), np.sin(th)])
    E = [(zp[x], zp[y]) for x in range(1, 81) for y in range(x + 1, 81) if R[x, y]]
    ax.add_collection(LineCollection(E, colors=GREY, linewidths=0.15, alpha=0.18))
    ax.add_collection(LineCollection([(zp[0], zp[x]) for x in inn], colors=BLUE, linewidths=0.4, alpha=0.55))
    lam_col = dict(zip(lamv, [BLUE, ORANGE, AQUA, INK])); mu_col = dict(zip(muv, (["#c3c8d0", GREY, INK, "#5b616e"] if len(muv) > 1 else [INK])))
    for x in inn: ax.scatter(*zp[x], s=5, c=lam_col[int(C[0, x])], lw=0, zorder=4)
    for x in out: ax.scatter(*zp[x], s=5, c=mu_col[int(C[0, x])], lw=0, zorder=4)
    ax.scatter([0], [0], s=14, c=INK, lw=0, zorder=5)
    k0 = len(inn)
    print(name, 'inner', lamv, 'outer', muv, 'srg', srg)
    desc = (f"$\\mathrm{{srg}}(81,{k0},{lamv[0]},{muv[0]})$" if srg else f"{k0}-regular, not srg")
    ax.set_title(f"$\\mathrm{{LS}}_{{{m}}}(36)$, {name}\n{desc}", fontsize=7, pad=2)
    ax.text(0, -2.15, f"$U_0$: {u0}", ha="center", fontsize=6)
    ax.set_xlim(-1.95, 1.95); ax.set_ylim(-2.35, 1.95); ax.set_aspect("equal"); ax.axis("off")
fig.savefig(os.path.join(OUT, "zero_rings.pdf"), bbox_inches="tight", pad_inches=0.03)
print("written", os.path.join(OUT, "zero_rings.pdf"))
