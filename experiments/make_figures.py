"""Figure generation.  Every figure answers a specific mathematical question.

fig_convex_packings.png  -- what does an optimal packing look like? (Thm 1)
fig_three_functions.png  -- f(n), q(n) and a(n-2) all coincide (Thm 1, Thm 2)
fig_configurations.png   -- nu(P) for many configurations vs the convex line (Conj. D)
fig_extremal_gap.png     -- convex position vs the best general configuration (Thm. C)
"""

from __future__ import annotations

import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon

from planetri.convex import (
    convex_nu_closed_form,
    convex_nu_dp_tables,
    convex_optimal_packing,
    convex_packing_avoids_edge,
)
from planetri.geom import PointSet, is_in_convex_position
from planetri.solver import nu_exact, optimal_packing
from planetri.trees import alpha_max_subcubic

RESULTS = os.path.join(os.path.dirname(__file__), "results")
FIGURES = os.path.join(os.path.dirname(__file__), "..", "figures")

INK = "#1a1a1a"
ACCENT = "#B03A2E"
BLUE = "#1F618D"
GREEN = "#1E8449"
GREY = "#7F8C8D"


def _ring(n: int, r: float = 1.0) -> np.ndarray:
    """Ring.
    
    Args:
        n:
        r (float):
    
    Returns:
        The computed result
    
    """
    a = 2 * math.pi * np.arange(n) / n
    return np.stack([r * np.cos(a), r * np.sin(a)], axis=1)


def fig_convex_packings() -> None:
    """Fig convex packings.
    
    """
    sizes = [9, 12, 15]
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 4.1))
    for ax, n in zip(axes, sizes):
        pts = _ring(n)
        pk = convex_optimal_packing(n)
        for i, t in enumerate(pk):
            tri = pts[list(t)]
            ax.add_patch(
                Polygon(
                    tri,
                    closed=True,
                    facecolor=plt.cm.Blues(0.18 + 0.5 * (i % 3) / 3.0),
                    edgecolor=BLUE,
                    linewidth=1.3,
                    alpha=0.85,
                    zorder=2,
                )
            )
        ax.plot(np.append(pts[:, 0], pts[0, 0]), np.append(pts[:, 1], pts[0, 1]),
                color=INK, linewidth=1.0, zorder=3)
        ax.scatter(pts[:, 0], pts[:, 1], s=22, color=INK, zorder=4)
        for k, (x, y) in enumerate(pts):
            ax.annotate(str(k), (x, y), textcoords="offset points", xytext=(0, 6),
                        ha="center", fontsize=7, color=GREY)
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_title(
            f"$n={n}$:  $|\\mathcal{{P}}|={len(pk)}=\\lfloor(2n-3)/3\\rfloor$",
            fontsize=11,
        )
    fig.suptitle(
        "Optimal edge-disjoint non-crossing triangle packings of convex $n$-gons",
        fontsize=13,
    )
    fig.tight_layout()
    fig.savefig(os.path.join(FIGURES, "fig_convex_packings.png"), dpi=190)
    plt.close(fig)


def fig_three_functions() -> None:
    """Fig three functions.
    
    """
    ns = np.arange(3, 121)
    p, q = convex_nu_dp_tables(121)
    f = np.array([p[n] for n in ns], dtype=float)
    qq = np.array([q[n] for n in ns], dtype=float)
    a = np.array([alpha_max_subcubic(n - 2) for n in ns], dtype=float)
    cf = np.array([convex_nu_closed_form(n) for n in ns], dtype=float)

    fig, ax = plt.subplots(figsize=(7.4, 5.0))
    ax.plot(ns, f, color=BLUE, lw=2.4, label="$f(n)$  (convex packings)")
    ax.plot(ns, a, color=GREEN, lw=1.8, ls="--",
            label="$a(n-2)$  (max indep. number of a subcubic tree)")
    ax.plot(ns, cf, color=ACCENT, lw=1.6, ls=":",
            label="$\\lfloor(2n-3)/3\\rfloor$")
    ax.plot(ns, qq, color=GREY, lw=1.8, ls="-.",
            label="$q(n)$  (a hull edge forbidden)")
    ax.plot(ns, (2 * ns - 4) / 3.0, color=GREY, lw=1.0, ls="--", alpha=0.6,
            label="$\\lfloor(2n-4)/3\\rfloor$")
    ax.set_xlabel("$n$")
    ax.set_ylabel("value")
    ax.set_title("Thms. 4.1 and 5.1: three quantities coincide", fontsize=12)
    ax.legend(fontsize=8.5, loc="upper left")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGURES, "fig_three_functions.png"), dpi=190)
    plt.close(fig)


def _families(n: int, rng) -> dict[str, np.ndarray]:
    """Families.
    
    Args:
        n:
        rng:
    
    Returns:
        The computed result
    
    """
    ang = 2 * math.pi * np.arange(n) / n
    ring = np.stack([np.cos(ang), np.sin(ang)], axis=1)
    out = {"convex": ring}
    out["convex + noise"] = ring + 0.06 * rng.normal(size=(n, 2))
    # an outer ring of n-1 points plus a single point in the middle
    outer = ring[: n - 1]
    out["ring + centre"] = np.vstack([outer, [[0.0, 0.0]]])
    # two concentric rings
    n_out = n // 2 + 1
    n_in = n - n_out
    ao = 2 * math.pi * np.arange(n_out) / n_out
    out_ring = np.stack([np.cos(ao), np.sin(ao)], axis=1)
    ai = 2 * math.pi * np.arange(n_in) / max(n_in, 1)
    r_in = 0.45 + 0.12 * rng.random(n_in)
    in_ring = np.stack([r_in * np.cos(ai), r_in * np.sin(ai)], axis=1)
    out["nested rings"] = np.vstack([out_ring, in_ring]) if n_in else out_ring
    out["uniform random"] = rng.normal(size=(n, 2))
    k = 3
    cent = rng.normal(size=(k, 2)) * 2.5
    idx = rng.integers(0, k, size=n)
    out["3 clusters"] = cent[idx] + rng.normal(size=(n, 2)) * 0.25
    return out


def fig_configurations() -> None:
    """Fig configurations.
    
    """
    rng = np.random.default_rng(7)
    ns = list(range(4, 13))
    series: dict[str, list[float]] = {}
    convex_line = [convex_nu_closed_form(n) for n in ns]
    for name in ("convex", "convex + noise", "ring + centre", "nested rings",
                 "uniform random", "3 clusters"):
        series[name] = []
    for n in ns:
        fams = _families(n, rng)
        for name in series:
            coords = fams[name]
            if name != "convex":
                # generic perturbation: the families are meant to represent
                # configurations in general position, and exact degeneracies
                # (e.g. a centre plus an antipodal pair) are not the object of
                # study
                coords = coords + 1e-3 * rng.normal(size=coords.shape)
            coords = np.round(coords, 5)
            ps = PointSet([tuple(r) for r in coords])
            if ps.n != n or not _gp(ps):
                series[name].append(float("nan"))
                continue
            series[name].append(nu_exact(ps))

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.6))
    styles = {
        "convex": (ACCENT, "-", 2.2),
        "convex + noise": ("#C0392B", "--", 1.3),
        "ring + centre": (BLUE, "-.", 1.3),
        "nested rings": (BLUE, ":", 1.3),
        "uniform random": (GREEN, "-", 1.5),
        "3 clusters": ("#7D3C98", "--", 1.5),
    }
    ax = axes[0]
    for name, ys in series.items():
        c, ls, lw = styles[name]
        ax.plot(ns, ys, color=c, ls=ls, lw=lw, marker="o", ms=3.5, label=name)
    ax.plot(ns, convex_line, color=INK, lw=2.4, ls="-", marker="s", ms=5,
            label="$\\lfloor(2n-3)/3\\rfloor$")
    ax.set_xlabel("$n$")
    ax.set_ylabel(r"$\nu(P)$")
    ax.set_title(r"$\nu(P)$ for six configuration families", fontsize=11)
    ax.legend(fontsize=7.5, ncol=2)
    ax.grid(alpha=0.25)

    ax = axes[1]
    for name, ys in series.items():
        c, ls, lw = styles[name]
        excess = [y - c0 for y, c0 in zip(ys, convex_line)]
        ax.plot(ns, excess, color=c, ls=ls, lw=lw, marker="o", ms=3.5, label=name)
    ax.axhline(0, color=INK, lw=2.0)
    ax.set_xlabel("$n$")
    ax.set_ylabel(r"$\nu(P)-\lfloor(2n-3)/3\rfloor$")
    ax.set_title("Conjecture D: the excess never goes below $0$", fontsize=11)
    ax.legend(fontsize=7.5, ncol=2)
    ax.grid(alpha=0.25)

    fig.suptitle("Convex position appears to be the worst configuration", fontsize=13)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGURES, "fig_configurations.png"), dpi=190)
    plt.close(fig)

    with open(os.path.join(RESULTS, "configuration_scan.csv"), "w") as fh:
        fh.write("n," + ",".join(series.keys()) + ",convex_value\n")
        for i, n in enumerate(ns):
            fh.write(str(n) + "," + ",".join(str(series[k][i]) for k in series)
                     + "," + str(convex_line[i]) + "\n")


def _gp(ps) -> bool:
    """Gp.
    
    Args:
        ps:
    
    Returns:
        bool: Result of type bool
    
    """
    from planetri.geom import orient

    n = ps.n
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                if orient(ps[i], ps[j], ps[k]) == 0:
                    return False
    return True


def fig_extremal_gap() -> None:
    """Fig extremal gap.
    
    """
    convex6 = [tuple(map(tuple, np.round(_ring(6), 4)))]
    nonconvex = [(-300, -125), (0, 100), (300, -125), (-75, 25), (75, -50), (0, 175)]

    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.3))
    for ax, coords, title, val in (
        (axes[0], convex6[0], "convex hexagon", convex_nu_closed_form(6)),
        (axes[1], nonconvex, "octahedron (Eulerian triangulation)", 4),
    ):
        arr = np.array(coords, dtype=float)
        ps = PointSet([tuple(r) for r in arr])
        pk = optimal_packing(ps)
        for i, t in enumerate(pk):
            tri = arr[list(t)]
            ax.add_patch(Polygon(tri, closed=True,
                                 facecolor=plt.cm.Oranges(0.2 + 0.5 * (i % 3) / 3.0),
                                 edgecolor=ACCENT, lw=1.4, alpha=0.9, zorder=2))
        edges = sorted({tuple(sorted(e)) for t in pk for e in
                        ((t[0], t[1]), (t[0], t[2]), (t[1], t[2]))})
        for a, b in edges:
            ax.plot(arr[[a, b], 0], arr[[a, b], 1], color=ACCENT, lw=0.8, alpha=0.4, zorder=1)
        ax.scatter(arr[:, 0], arr[:, 1], s=26, color=INK, zorder=4)
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_title(f"{title}\n$\\nu(P)={val}$", fontsize=11)

    fig.suptitle("Theorem C: the same $n=6$, two very different answers", fontsize=13)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGURES, "fig_extremal_gap.png"), dpi=190)
    plt.close(fig)


def main() -> None:
    """Entry point — parse arguments and run the main computation.
    
    """
    os.makedirs(FIGURES, exist_ok=True)
    os.makedirs(RESULTS, exist_ok=True)
    fig_convex_packings()
    fig_three_functions()
    fig_configurations()
    fig_extremal_gap()
    print("figures written to", os.path.abspath(FIGURES))


if __name__ == "__main__":
    main()