'''Experiment 1: convex-position exact values, by recurrence and by brute force.

Two independent computations of f(n) = max number of pairwise non-crossing,
edge-disjoint triangles spanned by n points in convex position:

(a) the two-parameter dynamic program ``convex_nu_dp``;
(b) exhaustive enumeration of every triangulation of the n-gon (C_{n-2} of
    them), computing for each triangulation the exact maximum number of
    edge-disjoint 3-cycles by integer programming, and taking the maximum.

(a) and (b) agreeing is strong evidence that the recurrence is right; (b) is
exponential but completely independent of (a).
'''  # noqa: D400,D401

from __future__ import annotations

import csv
import os
import sys
import time
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp

# Ensure the repository's ``src`` directory is on the import path.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

from planetri.convex import (
    CATALAN,
    convex_nu_closed_form,
    convex_nu_dp,
    enumerate_triangulations,
    triangulation_edges,
    triangulation_triangle_count,
)
from planetri.trees import alpha_max_subcubic

RESULTS = Path(__file__).with_name("results")
N_ENUMERATED = 10

__all__ = [
    "max_edge_disjoint_triangles",
    "main",
]


def max_edge_disjoint_triangles(diagonals: List[Tuple[int, int]], n: int) -> int:
    """Return the exact maximum number of edge‑disjoint triangles for a triangulation.

    The function builds a conflict graph where each vertex corresponds to a
    triangle of the triangulation and an edge indicates that the two triangles
    share a polygon edge.  The problem is then a maximum independent set on a
    bipartite conflict graph, which we solve via a mixed‑integer linear program
    (MILP) using :func:`scipy.optimize.milp`.

    Parameters
    ----------
    diagonals:
        List of interior diagonals that, together with the polygon boundary,
        define a triangulation of the convex *n*-gon.
    n:
        Number of vertices of the convex polygon.

    Returns
    -------
    int
        The size of a largest set of pairwise edge‑disjoint triangles.
    """
    # Compute all triangles of the triangulation.
    tris = triangulation_triangle_count(diagonals, n)
    m = len(tris)
    if m == 0:
        return 0

    # Map each edge to the list of triangle indices that contain it.
    owners: Dict[Tuple[int, int], List[int]] = {}
    for idx, t in enumerate(tris):
        for e in ((t[0], t[1]), (t[0], t[2]), (t[1], t[2])):
            owners.setdefault(e, []).append(idx)

    # Build the set of conflicting triangle pairs.
    conflict: set[Tuple[int, int]] = set()
    for lst in owners.values():
        for i, a in enumerate(lst):
            for b in lst[i + 1 :]:
                conflict.add((min(a, b), max(a, b)))

    # If there are no conflicts the answer is simply the number of triangles.
    if not conflict:
        return m

    # Construct the constraint matrix A where each row forbids selecting both
    # triangles of a conflicting pair.
    A = np.zeros((len(conflict), m), dtype=float)
    for row, (a, b) in enumerate(sorted(conflict)):
        A[row, a] = 1.0
        A[row, b] = 1.0

    # Linear constraints: each conflicting pair may contribute at most 1 to the
    # sum of selected triangles.
    cons = LinearConstraint(A, -np.inf, np.ones(len(conflict)))

    # Solve the MILP: maximise c·x with x ∈ {0,1}^m subject to the constraints.
    res = milp(
        c=-np.ones(m),  # maximise → minimise negative
        constraints=cons,
        integrality=np.ones(m),
        bounds=Bounds(0, 1),
    )

    if not res.success:
        # In the unlikely event that the solver fails, fall back to a safe
        # heuristic: return the trivial lower bound (0).
        return 0
    # ``res.fun`` is the minimized objective value (negative of the optimum).
    return int(round(-res.fun))


def main() -> int:
    """Run the experiment and write a CSV file with the results.

    Returns
    -------
    int
        Exit status compatible with ``sys.exit`` – ``0`` for success,
        ``1`` for mismatches.
    """
    RESULTS.mkdir(parents=True, exist_ok=True)
    dp = convex_nu_dp(400)

    rows: List[Dict[str, int]] = []
    for n in range(3, 31):
        rows.append(
            {
                "n": n,
                "f_dp": dp[n],
                "closed_form": convex_nu_closed_form(n),
                "alpha_max_subcubic_n_minus_2": alpha_max_subcubic(n - 2),
                "agree": int(dp[n] == convex_nu_closed_form(n)),
            }
        )

    # Exhaustive cross‑check for small n.
    brute: Dict[int, Tuple[int, int, float]] = {}
    for n in range(3, N_ENUMERATED + 1):
        t0 = time.time()
        best = 0
        count = 0
        for diagonals in enumerate_triangulations(n):
            count += 1
            val = max_edge_disjoint_triangles(diagonals, n)
            if val > best:
                best = val
        brute[n] = (best, count, time.time() - t0)
        print(
            f"n={n:2d}  brute={best:2d}  dp={dp[n]:2d}  closed_form={convex_nu_closed_form(n):2d}"
            f"  #triangulations={count:6d} (Catalan={CATALAN(n - 2):6d})  {time.time() - t0:6.1f}s",
            flush=True,
        )

    mismatches: List[int] = []
    for n, (best, _, _) in brute.items():
        if best != dp[n] or best != convex_nu_closed_form(n):
            mismatches.append(n)

    csv_path = RESULTS / "convex_exact.csv"
    with csv_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        for r in rows:
            writer.writerow(r)

    print()
    print(
        "DP matches closed form floor((2n-3)/3) for all 3<=n<=30 :",
        all(r["agree"] for r in rows),
    )
    print(
        f"brute force (n<={N_ENUMERATED}) matches DP                        : {len(mismatches) == 0}"
    )
    if mismatches:
        print("MISMATCHES at n =", sorted(set(mismatches)))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
