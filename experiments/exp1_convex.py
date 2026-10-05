"""Experiment 1: convex-position exact values, by recurrence and by brute force.

Two independent computations of f(n) = max number of pairwise non-crossing,
edge-disjoint triangles spanned by n points in convex position:

(a) the two-parameter dynamic program ``convex_nu_dp``;
(b) exhaustive enumeration of every triangulation of the n-gon (C_{n-2} of
    them), computing for each triangulation the exact maximum number of
    edge-disjoint 3-cycles by integer programming, and taking the maximum.

(a) and (b) agreeing is strong evidence that the recurrence is right; (b) is
exponential but completely independent of (a).
"""

from __future__ import annotations

import csv
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp

from planetri.convex import (
    CATALAN,
    convex_nu_closed_form,
    convex_nu_dp,
    enumerate_triangulations,
    triangulation_edges,
    triangulation_triangle_count,
)
from planetri.trees import alpha_max_subcubic

RESULTS = os.path.join(os.path.dirname(__file__), "results")
N_ENUMERATED = 10


def max_edge_disjoint_triangles(diagonals, n: int) -> int:
    """Exact max # pairwise edge-disjoint 3-cycles of one triangulation."""
    tris = triangulation_triangle_count(diagonals, n)
    m = len(tris)
    if m == 0:
        return 0
    owners: dict[tuple[int, int], list[int]] = {}
    for idx, t in enumerate(tris):
        for e in ((t[0], t[1]), (t[0], t[2]), (t[1], t[2])):
            owners.setdefault(e, []).append(idx)
    conflict = set()
    for _, lst in owners.items():
        for i, item in enumerate(lst):
            for j in range(i + 1, len(lst)):
                conflict.add((min(item, lst[j]), max(item, lst[j])))
    A = np.zeros((len(conflict), m))
    for row, (a, b) in enumerate(sorted(conflict)):
        A[row, a] = 1.0
        A[row, b] = 1.0
    cons = LinearConstraint(A, -np.inf, np.ones(len(conflict))) if conflict else None
    res = milp(c=-np.ones(m), constraints=cons, integrality=np.ones(m), bounds=Bounds(0, 1))
    return int(round(-res.fun))


def main() -> None:
    """Entry point — parse arguments and run the main computation.
    
    Returns:
        int: Result of type int
    
    """
    os.makedirs(RESULTS, exist_ok=True)
    dp = convex_nu_dp(400)

    rows = []
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

    # Exhaustive cross-check.
    brute = {}
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

    mismatches = []
    for n, (best, count, dt) in brute.items():
        if best != dp[n]:
            mismatches.append(n)
        if best != convex_nu_closed_form(n):
            mismatches.append(n)

    path = os.path.join(RESULTS, "convex_exact.csv")
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(r)

    print()
    print("DP matches closed form floor((2n-3)/3) for all 3<=n<=30 :",
          all(r["agree"] for r in rows))
    print(f"brute force (n<={N_ENUMERATED}) matches DP                        : {...}")
    if mismatches:
        print("MISMATCHES at n =", sorted(set(mismatches)))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())