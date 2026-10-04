"""Exact maximisation of non-crossing, edge-disjoint triangle packings.

A *packing* of a point set ``P`` is a family of 3-element subsets ("triangles")
such that

1. no two triangles share an edge (two common vertices), and
2. no edge of one triangle crosses an edge of another triangle.

Condition 2 is exactly what makes the union of the triangles a plane straight
line graph; condition 1 is a packing condition on the edges.

:func:`nu_exact` computes the optimum by maximum-weight independent set in the
conflict graph, solved as a 0/1 integer program with HiGHS through
``scipy.optimize.milp``.  The result is exact for integral input.
"""

from __future__ import annotations

import itertools
from typing import Iterable, Sequence

import numpy as np
from scipy.optimize import LinearConstraint, Bounds, milp

from .geom import Point, PointSet, segments_cross

__all__ = [
    "Triangle",
    "all_triangles",
    "conflicting",
    "conflict_pairs",
    "is_valid_packing",
    "nu_exact",
    "optimal_packing",
    "greedy_lower_bound",
]


class Triangle(tuple):
    """A triple of vertex indices, stored in increasing order."""

    __slots__ = ()

    def __new__(cls, a: int, b: int, c: int):
        return super().__new__(cls, sorted((a, b, c)))

    @property
    def a(self) -> int:
        return self[0]

    @property
    def b(self) -> int:
        return self[1]

    @property
    def c(self) -> int:
        return self[2]

    def edges(self) -> tuple[tuple[int, int], ...]:
        a, b, c = self
        return ((a, b), (a, c), (b, c))

    def __repr__(self) -> str:  # pragma: no cover - debugging helper
        return f"Triangle({self[0]}, {self[1]}, {self[2]})"


def all_triangles(n: int) -> list[Triangle]:
    """All ``C(n,3)`` triangles on vertex set ``{0,...,n-1}``."""
    return [Triangle(*t) for t in itertools.combinations(range(n), 3)]


def conflicting(t1: Triangle, t2: Triangle, pts: Sequence[Point]) -> bool:
    """Do ``t1`` and ``t2`` violate condition 1 or condition 2?"""
    if len(set(t1) & set(t2)) >= 2:
        return True
    p1, p2, p3 = pts[t1[0]], pts[t1[1]], pts[t1[2]]
    q1, q2, q3 = pts[t2[0]], pts[t2[1]], pts[t2[2]]
    for a, b in ((p1, p2), (p1, p3), (p2, p3)):
        for c, d in ((q1, q2), (q1, q3), (q2, q3)):
            if segments_cross(a, b, c, d):
                return True
    return False


def conflict_pairs(triangles: Sequence[Triangle], pts: Sequence[Point]) -> list[tuple[int, int]]:
    """Index pairs of conflicting triangles (each pair listed once)."""
    out: list[tuple[int, int]] = []
    for i in range(len(triangles)):
        ti = triangles[i]
        for j in range(i + 1, len(triangles)):
            if conflicting(ti, triangles[j], pts):
                out.append((i, j))
    return out


def is_valid_packing(packing: Iterable[Triangle], pts: Sequence[Point]) -> bool:
    """Verify the two packing conditions directly (used by the test suite)."""
    packing = list(packing)
    for a, b in itertools.combinations(range(len(packing)), 2):
        if conflicting(packing[a], packing[b], pts):
            return False
    return True


def nu_exact(point_set: PointSet, *, time_limit: float | None = None) -> int:
    """Exact value of the maximum non-crossing edge-disjoint triangle packing.

    Parameters
    ----------
    point_set:
        Points in general position (no three collinear).
    time_limit:
        Optional solver time limit in seconds.  If the solver hits the limit the
        *best incumbent* is returned; the function then no longer guarantees an
        exact answer and sets the attribute ``last_solve_was_exact``.
    """
    n = point_set.n
    if n < 3:
        solver.nu_exact.last_solve_was_exact = True
        return 0
    triangles = all_triangles(n)
    pairs = conflict_pairs(triangles, point_set.points)
    m = len(triangles)
    # Variables x_t in {0,1}; maximise sum x_t subject to x_s + x_t <= 1.
    cost = -np.ones(m)
    integrality = np.ones(m)
    bounds = Bounds(0, 1)
    if pairs:
        A = np.zeros((len(pairs), m))
        for row, (s, t) in enumerate(pairs):
            A[row, s] = 1.0
            A[row, t] = 1.0
        constraints = LinearConstraint(A, -np.inf, np.ones(len(pairs)))
    else:
        constraints = None
    options = {"presolve": True}
    if time_limit is not None:
        options["time_limit"] = float(time_limit)
    res = milp(
        c=cost,
        constraints=constraints,
        integrality=integrality,
        bounds=bounds,
        options=options,
    )
    nu_exact.last_solve_was_exact = bool(res.status == 0)
    return int(round(-res.fun)) if res.fun is not None else 0


nu_exact.last_solve_was_exact = True  # type: ignore[attr-defined]


def optimal_packing(point_set: PointSet, *, time_limit: float | None = None) -> list[Triangle]:
    """An optimal packing attaining :func:`nu_exact`."""
    n = point_set.n
    if n < 3:
        return []
    triangles = all_triangles(n)
    pairs = conflict_pairs(triangles, point_set.points)
    m = len(triangles)
    cost = -np.ones(m)
    integrality = np.ones(m)
    bounds = Bounds(0, 1)
    if pairs:
        A = np.zeros((len(pairs), m))
        for row, (s, t) in enumerate(pairs):
            A[row, s] = 1.0
            A[row, t] = 1.0
        constraints = LinearConstraint(A, -np.inf, np.ones(len(pairs)))
    else:
        constraints = None
    options = {"presolve": True}
    if time_limit is not None:
        options["time_limit"] = float(time_limit)
    res = milp(
        c=cost,
        constraints=constraints,
        integrality=integrality,
        bounds=bounds,
        options=options,
    )
    x = np.round(res.x).astype(bool)
    return [triangles[i] for i in range(m) if x[i]]


def greedy_lower_bound(point_set: PointSet, order: Sequence[int] | None = None) -> list[Triangle]:
    """A fast, always-valid packing (used as a sanity net and in stress tests).

    Triangles are considered in the given vertex order (``order`` defaults to
    lexicographic order) and accepted greedily if compatible with everything
    accepted so far.
    """
    n = point_set.n
    triangles = all_triangles(n)
    if order is not None:
        triangles = sorted(triangles, key=lambda t: tuple(order[i] for i in t))
    chosen: list[Triangle] = []
    for t in triangles:
        if all(not conflicting(t, u, point_set.points) for u in chosen):
            chosen.append(t)
    return chosen