"""Test suite for the `planetri` package.

The tests are written to catch *mathematical* errors rather than to cover code:
they check exact closed forms against recurrences, validate constructed
packings geometrically (not merely combinatorially), and verify the structural
lemmas that the proofs rely on.
"""

from __future__ import annotations

import math

import pytest

from planetri.convex import (
    CATALAN,
    convex_nu_closed_form,
    convex_nu_dp,
    convex_nu_dp_tables,
    convex_optimal_packing,
    convex_packing_avoids_edge,
    enumerate_triangulations,
    triangulation_edges,
    triangulation_triangle_count,
)
from planetri.geom import Point, PointSet, convex_hull_order, is_in_convex_position, orient, segments_cross
from planetri.solver import Triangle, conflicting, is_valid_packing, nu_exact
from planetri.trees import alpha_max_subcubic, alpha_tree, weak_dual


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def regular_convex(n: int) -> PointSet:
    pts = []
    for i in range(n):
        a = 2 * math.pi * i / n
        pts.append((round(1000 * math.cos(a), 6), round(1000 * math.sin(a), 6)))
    return PointSet(pts)


def edge_set(tri) -> set[tuple[int, int]]:
    return {tuple(sorted(e)) for e in tri.edges()}


# --------------------------------------------------------------------------- #
# predicates
# --------------------------------------------------------------------------- #
def test_orient_basic():
    assert orient(Point(0, 0), Point(1, 0), Point(0, 1)) == 1
    assert orient(Point(0, 0), Point(0, 1), Point(1, 0)) == -1
    assert orient(Point(0, 0), Point(1, 1), Point(2, 2)) == 0


def test_shared_endpoint_is_not_a_crossing():
    a, b, c = Point(0, 0), Point(2, 0), Point(1, 1)
    assert segments_cross(a, b, b, c) is False
    assert segments_cross(a, b, a, c) is False


def test_proper_crossing_detected():
    # unit square diagonals cross; opposite sides do not
    p = [Point(0, 0), Point(1, 0), Point(1, 1), Point(0, 1)]
    assert segments_cross(p[0], p[2], p[1], p[3]) is True
    assert segments_cross(p[0], p[1], p[2], p[3]) is False


def test_convex_hull_detection():
    assert is_in_convex_position(regular_convex(7))
    sq = PointSet([(0, 0), (2, 0), (2, 2), (0, 2), (1, 1)])
    assert not is_in_convex_position(sq)
    assert len(convex_hull_order(sq)) == 4


def test_pointset_exact_decimal_arithmetic():
    """Decimals must be handled exactly, not as binary doubles.

    With naive float arithmetic ``0.1 + 0.2 != 0.3``; here the three points are
    exactly collinear, and the predicate must report that.
    """
    ps = PointSet([(0, 0), (0.1, 0.2), (0.3, 0.6)])
    assert orient(ps[0], ps[1], ps[2]) == 0
    assert 0.1 + 0.2 != 0.3  # the naive float computation that we avoid


# --------------------------------------------------------------------------- #
# the convex recurrence and its closed form  (Theorem 1)
# --------------------------------------------------------------------------- #
def test_dp_matches_closed_form():
    tab = convex_nu_dp(400)
    for m in range(3, 401):
        assert tab[m] == convex_nu_closed_form(m), m


def test_q_table_matches_closed_form():
    _p, q = convex_nu_dp_tables(400)
    for m in range(2, 401):
        assert q[m] == (2 * m - 4) // 3, m


def test_constructed_packing_is_optimal_and_valid():
    for n in range(3, 60):
        ps = regular_convex(n)
        pk = convex_optimal_packing(n)
        assert len(pk) == convex_nu_closed_form(n), n
        assert is_valid_packing([Triangle(*t) for t in pk], ps.points), n


def test_constructed_avoiding_packing():
    for n in range(4, 60):
        ps = regular_convex(n)
        pk = convex_packing_avoids_edge(n)
        assert is_valid_packing([Triangle(*t) for t in pk], ps.points), n
        used = {tuple(sorted(e)) for t in pk for e in ((t[0], t[1]), (t[0], t[2]), (t[1], t[2]))}
        assert (0, n - 1) not in used, n


def test_two_triangles_in_a_quadrilateral_always_share_an_edge():
    """sanity check of the base case n = 4: nu = 1"""
    assert convex_nu_closed_form(4) == 1


# --------------------------------------------------------------------------- #
# triangulations of the convex polygon
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("n", range(3, 11))
def test_triangulation_enumeration_is_catalan(n):
    trs = list(enumerate_triangulations(n))
    assert len(trs) == CATALAN(n - 2)
    assert len(set(trs)) == len(trs)
    for t in trs:
        assert len(t) == n - 3


def test_triangulation_edge_normalisation():
    n = 5
    diag = [(1, 3), (3, 4)]
    edges = triangulation_edges(diag, n)
    for (a, b) in edges:
        assert a < b


@pytest.mark.parametrize("n", range(3, 10))
def test_every_three_cycle_of_a_triangulation_is_a_face(n):
    """Lemma L1 -- the crux of the weak-dual reduction."""
    for diagonals in enumerate_triangulations(n):
        tris = set(triangulation_triangle_count(diagonals, n))
        faces, adj = weak_dual(diagonals, n)
        assert tris == {tuple(f) for f in faces}
        # the weak dual really is a tree on n-2 vertices of degree <= 3
        assert len(faces) == n - 2
        assert all(len(a) <= 3 for a in adj)
        # the dual adjacency graph is a tree: n-2 vertices, n-3 edges
        assert sum(len(a) for a in adj) == 2 * (n - 3)


@pytest.mark.parametrize("n", range(3, 10))
def test_triangulation_max_packing_equals_max_weak_dual_alpha(n):
    """f(n) = max over triangulations of alpha(weak dual)."""
    best_alpha = 0
    for diagonals in enumerate_triangulations(n):
        faces, adj = weak_dual(diagonals, n)
        best_alpha = max(best_alpha, alpha_tree({i: adj[i] for i in range(len(faces))}))
    assert best_alpha == convex_nu_closed_form(n)


# --------------------------------------------------------------------------- #
# subcubic trees  (Theorem 2)
# --------------------------------------------------------------------------- #
def test_alpha_max_subcubic_closed_form():
    for N in range(1, 41):
        assert alpha_max_subcubic(N) == (2 * N + 1) // 3, N


def test_alpha_max_subcubic_matches_convex_extremal():
    """f(n) = a(n-2) exactly."""
    for n in range(3, 41):
        assert alpha_max_subcubic(n - 2) == convex_nu_closed_form(n), n


def test_alpha_max_subcubic_lower_bound_random_search():
    """The DP value must dominate a wide random search over subcubic trees."""
    import random

    rng = random.Random(12345)
    for N in range(4, 30):
        best = 0
        for _ in range(3000):
            adj = {i: [] for i in range(N)}
            ok = True
            for v in range(1, N):
                cands = [u for u in range(N) if u != v and len(adj[u]) < 3 and len(adj[v]) < 3]
                if not cands:
                    ok = False
                    break
                u = rng.choice(cands)
                adj[u].append(v)
                adj[v].append(u)
            if ok:
                best = max(best, alpha_tree(adj))
        assert best <= alpha_max_subcubic(N), (N, best, alpha_max_subcubic(N))


def test_alpha_tree_on_known_trees():
    # path on 6 vertices
    path = {i: ([i - 1] if i else []) + ([i + 1] if i < 5 else []) for i in range(6)}
    assert alpha_tree(path) == 3
    # star on 5 vertices is not subcubic, but alpha_tree must still be right
    star = {i: ([0] if i else [1, 2, 3, 4]) for i in range(5)}
    assert alpha_tree(star) == 4
    # the subdivided claw K_{1,3} (4 vertices): alpha = 3
    claw = {0: [1, 2, 3], 1: [0], 2: [0], 3: [0]}
    assert alpha_tree(claw) == 3
    assert alpha_max_subcubic(4) == 3


# --------------------------------------------------------------------------- #
# the exact solver on small instances
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("n", range(3, 10))
def test_nu_exact_agrees_with_convex_formula(n):
    ps = regular_convex(n)
    assert nu_exact(ps) == convex_nu_closed_form(n)


def test_nu_exact_on_a_non_convex_configuration_beats_convex():
    """The octahedron (an Eulerian triangulation) drawn with straight lines.

    This 6-point set supports 4 = n-2 edge-disjoint non-crossing triangles,
    whereas convex position supports only 3.
    """
    ps = PointSet([(-300, -125), (0, 100), (300, -125), (-75, 25), (75, -50), (0, 175)])
    assert not is_in_convex_position(ps)
    assert nu_exact(ps) == 4 == 6 - 2
    assert convex_nu_closed_form(6) == 3


def test_conflicting_is_symmetric_and_detects_shared_edges():
    ps = regular_convex(6)
    t1, t2 = Triangle(0, 1, 2), Triangle(0, 1, 3)
    assert conflicting(t1, t2, ps.points)
    assert conflicting(t2, t1, ps.points)
    assert not conflicting(Triangle(0, 1, 2), Triangle(2, 3, 4), ps.points)


def test_nu_is_at_most_n_minus_two():
    for n in range(3, 9):
        assert nu_exact(regular_convex(n)) <= n - 2