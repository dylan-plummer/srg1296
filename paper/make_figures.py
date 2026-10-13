#!/usr/bin/env python3
"""Figures for the paper, rebuilt from ../connsets/*.txt with independent field arithmetic.
usage: python3 make_figures.py   (writes figures/rings.pdf)"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
CONN = os.path.join(HERE, "..", "connsets")
OUT = os.path.join(HERE, "figures")
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.family": "serif", "font.size": 8, "axes.linewidth": 0.5,
                     "xtick.major.width": 0.5, "ytick.major.width": 0.5, "pdf.fonttype": 42})


def powers(p, poly, n):
    cur = [1] + [0] * (n - 1); out = []
    for _ in range(p ** n - 1):
        out.append(tuple(cur)); top = cur[-1]; cur = [0] + cur[:-1]
        cur = [(c - top * poly[i]) % p for i, c in enumerate(cur)]
    return out


A16 = powers(2, [1, 0, 0, 1, 1], 4)          # alpha^i, F_16 = F_2[x]/(x^4+x^3+1)
B81 = powers(3, [2, 0, 0, 1, 1], 4)          # beta^j,  F_81 = F_3[y]/(y^4+y^3+2)
encA = lambda t: sum(c << s for s, c in enumerate(t))
encB = lambda t: sum(c * 3 ** s for s, c in enumerate(t))
logA = np.full(16, -1); logB = np.full(81, -1)
for i, t in enumerate(A16): logA[encA(t)] = i
for j, t in enumerate(B81): logB[encB(t)] = j
SUBB = np.array([[sum(((x // 3 ** s) % 3 - (y // 3 ** s) % 3) % 3 * 3 ** s for s in range(4)) for y in range(81)] for x in range(81)])


def load(name):
    L = [l.split() for l in open(os.path.join(CONN, name + ".txt")) if l.strip() and not l.startswith("#")]
    prm = tuple(map(int, L[0])); el = [(-1 if i == "-" else int(i), -1 if j == "-" else int(j)) for i, j in L[1:]]
    ids = [(0 if i < 0 else encA(A16[i])) * 81 + (0 if j < 0 else encB(B81[j])) for i, j in el]
    inD = np.zeros(1296, bool); inD[ids] = True
    u = np.arange(1296)
    diff = ((u[:, None] // 81) ^ (u[None, :] // 81)) * 81 + SUBB[(u % 81)[:, None], (u % 81)[None, :]]
    A = inD[diff] & ~np.eye(1296, dtype=bool)
    A2 = A.astype(np.int64) @ A.astype(np.int64)
    k, lam, mu = prm[1:]
    assert (A.sum(1) == k).all() and (A2[A] == lam).all() and (A2[~A & ~np.eye(1296, dtype=bool)] == mu).all()
    return prm, el, inD, A


G = {"560": load("srg1296_560_244_240_class1"), "592": load("srg1296_592_268_272_class1"),
     "595": load("srg1296_595_274_272_class1")}
NAME = {"560": r"$\mathrm{srg}(1296,560,244,240)$", "592": r"$\mathrm{srg}(1296,592,268,272)$",
        "595": r"$\mathrm{srg}(1296,595,274,272)$"}
AQUA = "#1baf7a"

BLUE, ORANGE, INK = "#2a78d6", "#eb6834", "#1b1f27"

BLUE, ORANGE, INK = "#2a78d6", "#eb6834", "#1b1f27"

# ---- Figure: ring layout with edges, and a zoom into the ring R_0 ----------------------------
from matplotlib.collections import LineCollection
from matplotlib.patches import ConnectionPatch, Rectangle
GREY = "#9aa1ad"
def ring_pos():
    P = np.zeros((1296, 2))
    for v in range(1296):
        a, b = divmod(v, 81); la, lb = logA[a], logB[b]
        cx, cy = (0, 0) if la < 0 else (0.71 * np.cos(2 * np.pi * la / 15 + np.pi / 2), 0.71 * np.sin(2 * np.pi * la / 15 + np.pi / 2))
        P[v] = (cx, cy) if lb < 0 else (cx + 0.135 * np.cos(2 * np.pi * lb / 80), cy + 0.135 * np.sin(2 * np.pi * lb / 80))
    return P
P = ring_pos()
ring_of = np.arange(1296) // 81
zpos = np.array([(0.0, 0.0) if logB[b] < 0 else (np.cos(2 * np.pi * logB[b] / 80), np.sin(2 * np.pi * logB[b] / 80)) for b in range(81)])
fig, ax = plt.subplots(3, 3, figsize=(6.9, 7.4), gridspec_kw={"wspace": 0.05, "hspace": 0.28})
v1 = encB(B81[0])                                   # the vertex (0, beta^0)
for row, (key, tag) in enumerate([("560", "I"), ("592", "II"), ("595", "III")]):
    A = G[key][3]; k, lam, mu = G[key][0][1:]
    u = 0; nb = np.where(A[u])[0]
    # (left) the edges at (0,0)
    a0 = ax[row, 0]
    inside = nb[ring_of[nb] == 0]; outside = nb[ring_of[nb] != 0]
    a0.add_collection(LineCollection([(P[u], P[x]) for x in outside], colors=BLUE, linewidths=0.25, alpha=0.35))
    a0.add_collection(LineCollection([(P[u], P[x]) for x in inside], colors=ORANGE, linewidths=0.6, alpha=0.9))
    a0.scatter(P[:, 0], P[:, 1], s=0.4, c="#c3c8d0", lw=0)
    a0.scatter(P[nb, 0], P[nb, 1], s=1.4, c=INK, lw=0)
    a0.scatter([P[u, 0]], [P[u, 1]], s=12, c=ORANGE, lw=0, zorder=4)
    a0.set_title(f"$\\Gamma_{{{tag}}}$: the {k} edges at $(0,0)$\n{len(inside)} in its ring, {len(outside)} to others", fontsize=7)
    # (middle) a non-adjacent pair and its common neighbours; frame the ring R_0
    w = next(x for x in range(1296) if ring_of[x] == 6 and x != u and not A[u, x])
    com = np.where(A[u] & A[w])[0]
    a1 = ax[row, 1]
    a1.add_collection(LineCollection([(P[u], P[x]) for x in com], colors=BLUE, linewidths=0.25, alpha=0.4))
    a1.add_collection(LineCollection([(P[w], P[x]) for x in com], colors=ORANGE, linewidths=0.25, alpha=0.4))
    a1.scatter(P[:, 0], P[:, 1], s=0.4, c="#c3c8d0", lw=0)
    a1.scatter(P[com, 0], P[com, 1], s=1.4, c=INK, lw=0)
    a1.scatter([P[u, 0]], [P[u, 1]], s=12, c=BLUE, lw=0, zorder=4)
    a1.scatter([P[w, 0]], [P[w, 1]], s=12, c=ORANGE, lw=0, zorder=4)
    a1.set_title(f"$\\Gamma_{{{tag}}}$: a non-adjacent pair,\n$\\mu={mu}$ common neighbours", fontsize=7)
    r0 = 0.16
    a1.add_patch(Rectangle((-r0, -r0), 2 * r0, 2 * r0, fill=False, ec=INK, lw=0.6, zorder=5))
    for a in (a0, a1):
        a.set_xlim(-0.9, 0.9); a.set_ylim(-0.9, 0.9); a.set_aspect("equal"); a.axis("off")
    # (right) zoom: the subgraph induced on R_0, layered from its centre (0,0)
    a2 = ax[row, 2]
    R = A[:81, :81]; C = R.astype(int) @ R.astype(int)
    inn = [x for x in range(81) if R[0, x]]; out = [x for x in range(81) if x and not R[0, x]]
    zp = {0: np.zeros(2)}
    if tag == "I":                               # windmill: pair x with its unique common neighbour -x
        prs = []
        for x in inn:
            y = [z for z in inn if R[x, z]]; assert len(y) == 1
            if x < y[0]: prs.append((x, y[0]))
        for t, (x, y) in enumerate(prs):
            th = 2 * np.pi * t / len(prs) + np.pi / 2
            zp[x] = np.array([np.cos(th - 0.15), np.sin(th - 0.15)]); zp[y] = np.array([np.cos(th + 0.15), np.sin(th + 0.15)])
    else:
        for x in inn:
            th = 2 * np.pi * logB[x] / 80; zp[x] = 1.2 * np.array([np.cos(th), np.sin(th)])
    for x in out:
        th = 2 * np.pi * logB[x] / 80; zp[x] = 1.75 * np.array([np.cos(th), np.sin(th)])
    E = [(zp[x], zp[y]) for x in range(81) for y in range(x + 1, 81) if R[x, y] and x and y]
    a2.add_collection(LineCollection(E, colors=GREY, linewidths=0.2, alpha=0.25 if tag == "I" else 0.1))
    a2.add_collection(LineCollection([(zp[0], zp[x]) for x in inn], colors=BLUE, linewidths=0.5, alpha=0.6))
    lamv = sorted(set(int(C[0, x]) for x in inn)); muv = sorted(set(int(C[0, x]) for x in out))
    lam_col = {l: c_ for l, c_ in zip(lamv, [BLUE, ORANGE, AQUA])}
    mu_col = {m: c_ for m, c_ in zip(muv, ([INK] if len(muv) == 1 else ["#c3c8d0", GREY, INK]))}
    if tag == "I":
        a2.add_collection(LineCollection([(zp[x], zp[y]) for x, y in prs], colors=ORANGE, linewidths=1.3))
        lam_col = {lamv[0]: ORANGE}
    for x in inn:
        a2.scatter([zp[x][0]], [zp[x][1]], s=10 if tag == "I" else 7, c=lam_col[int(C[0, x])], lw=0, zorder=4)
    for x in out:
        a2.scatter([zp[x][0]], [zp[x][1]], s=7, c=mu_col[int(C[0, x])], lw=0, zorder=4)
    a2.scatter([0], [0], s=18, c=INK, lw=0, zorder=5)
    if tag == "I":
        a2.set_title("zoom into $R_0$: the Brouwer–Haemers\ngraph $\\mathrm{srg}(81,20,1,6)$", fontsize=7)
        a2.text(0, -2.2, "inner: $\\lambda=1$ (triangles), outer: $\\mu=6$", ha="center", fontsize=6)
    else:
        a2.set_title(f"zoom into $R_0$: {len(inn)}-regular,\nnot strongly regular", fontsize=7)
        a2.text(0, -2.2, "inner: " + " / ".join(map(str, lamv)) + ", outer: " + " / ".join(map(str, muv)), ha="center", fontsize=6)
    a2.set_xlim(-2.0, 2.0); a2.set_ylim(-2.35, 2.0); a2.set_aspect("equal"); a2.axis("off")
    for (xa_, ya_), (xb_, yb_) in [((r0, r0), (-1.85, 1.85)), ((r0, -r0), (-1.85, -1.85))]:
        fig.add_artist(ConnectionPatch(xyA=(xa_, ya_), coordsA=a1.transData, xyB=(xb_, yb_), coordsB=a2.transData,
                                       color=INK, lw=0.4, ls=":"))
fig.savefig(os.path.join(OUT, "rings.pdf"), bbox_inches="tight", pad_inches=0.04); plt.close(fig)
print("figures written to", OUT)
