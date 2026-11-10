#!/usr/bin/env python3
"""Figure: orthogonality graphs of the subgroups of order 36 disjoint from D (graphs I.6, I.10, III.1, III.2).
Data: checks/spread_lattice.json (written by checks/spread_lattice.py). Output: OUT/orthogonality.pdf
usage: make_fig_orthogonality.py OUTDIR"""
import os, sys, json, itertools, math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "figures")
D = json.load(open(os.path.join(HERE, "checks", "spread_lattice.json")))
plt.rcParams.update({"font.family": "serif", "font.size": 8, "pdf.fonttype": 42, "mathtext.fontset": "cm"})
BLUE, GREY, INK, MUTED = "#2a78d6", "#b9bec7", "#1b1f27", "#5b616e"

PANELS = [("srg1296_560_244_240_class6", "I.6", 16), ("srg1296_560_244_240_class10", "I.10", 16),
          ("srg1296_595_274_272_class1", "III.1, III.3", 17), ("srg1296_595_274_272_class2", "III.2, III.4", 17)]


def max_clique(n, E):
    S = {tuple(sorted(e)) for e in E}
    return max((T for j in range(1, n + 1) for T in itertools.combinations(range(n), j)
                if all((a, b) in S for a, b in itertools.combinations(T, 2))), key=len)


fig, axes = plt.subplots(1, 4, figsize=(6.3, 1.85))
for ax, (key, name, m) in zip(axes, PANELS):
    r = D[key]; n = r["free_subgroups"]; E = [tuple(e) for e in r["orth_edges"]]
    C = max_clique(n, E); rest = [v for v in range(n) if v not in C]
    pos = {}
    for t, v in enumerate(C):                      # maximal orthogonal set on the inner polygon
        a = math.pi / 2 + 2 * math.pi * t / len(C); pos[v] = (0.62 * math.cos(a), 0.62 * math.sin(a))
    for t, v in enumerate(rest):                   # the other subgroups on an outer ring
        a = math.pi / 2 + math.pi / len(C) + 2 * math.pi * t / max(1, len(rest)); pos[v] = (1.12 * math.cos(a), 1.12 * math.sin(a))
    inC = set(C)
    for a, b in E:
        hot = a in inC and b in inC
        ax.plot(*zip(pos[a], pos[b]), color=BLUE if hot else GREY, lw=1.6 if hot else 0.8, zorder=1 if not hot else 2,
                solid_capstyle="round")
    for v in range(n):
        ax.scatter(*pos[v], s=46, color=BLUE if v in inC else "white", edgecolors=BLUE if v in inC else MUTED,
                   linewidths=1.0, zorder=3)
    j = len(C)
    ax.set_title(f"{name}", fontsize=8.5, color=INK, pad=3)
    ax.text(0, -1.38, f"{n} subgroups, {j} orthogonal\n"
                      rf"$\mathrm{{LS}}_{{{m}}}(36)\to\mathrm{{LS}}_{{{m + j}}}(36)$",
            ha="center", va="top", fontsize=7, color=INK, linespacing=1.5)
    ax.set_xlim(-1.3, 1.3); ax.set_ylim(-2.05, 1.3); ax.set_aspect("equal"); ax.axis("off")
fig.subplots_adjust(left=0.01, right=0.99, top=0.92, bottom=0.0, wspace=0.05)
os.makedirs(OUT, exist_ok=True)
fig.savefig(os.path.join(OUT, "orthogonality.pdf"), bbox_inches="tight", pad_inches=0.02)
fig.savefig(os.path.join(OUT, "orthogonality.png"), dpi=200, bbox_inches="tight", pad_inches=0.02)
print("written", OUT)
