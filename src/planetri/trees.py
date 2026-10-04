"""Weak duals of triangulations and independence numbers of subcubic trees.

For a point set ``P`` in general position, a *triangulation* ``T`` is a maximal
plane straight line graph on ``P``.  Its number of bounded triangular faces is

    |F(T)| = 2|P| - h - 2,      h = number of convex-hull vertices of P,

and the *weak dual* ``D(T)`` is the tree whose vertices are the faces of ``T``
and whose edges join faces sharing an edge.  ``D(T)`` is a tree of maximum
degree at most 3.

Two facts make the weak dual the right combinatorial object:

* If ``P`` is in convex position then every 3-cycle of ``T`` bounds a face of
  ``T``.  Indeed no vertex of ``P`` lies strictly inside a triangle spanned by
  three vertices, so a 3-cycle of ``T`` cannot contain another edge of ``T`` in
  its interior.  Hence for convex ``P`` the maximum number of pairwise
  edge-disjoint triangles equals exactly the independence number of ``D(T)``.
* In general position a 3-cycle of ``T`` may contain other points, so the
  number of edge-disjoint 3-cycles of ``T`` can exceed ``alpha(D(T))``.

The extremal quantity appearing in the convex analysis is

    a(N) = max { alpha(T) : T a tree with N vertices and maximum degree <= 3 }.

:func:`alpha_max_subcubic` computes ``a(N)`` exactly by dynamic programming over
rooted subcubic trees, keeping the Pareto frontier of the
"(root excluded, root included)" pair because the two optima are not attained
simultaneously.
"""

from __future__ import annotations

import functools
from typing import Dict, Iterable, List, Sequence, Set, Tuple

__all__ = [
    "weak_dual",
    "alpha_tree",
    "alpha_max_subcubic",
    "max_subcubic_tree",
    "IndependenceNumber",
]


# --------------------------------------------------------------------------- #
# Independence number of explicit trees
# --------------------------------------------------------------------------- #
def alpha_tree(adj: Dict[int, Sequence[int]]) -> int:
    """Independence number of a tree given by adjacency lists (iterative)."""
    if not adj:
        return 0
    nodes = sorted(adj)
    parent: Dict[int, int] = {nodes[0]: -1}
    order = [nodes[0]]
    stack = [nodes[0]]
    while stack:
        u = stack.pop()
        for w in adj.get(u, ()):
            if w not in parent:
                parent[w] = u
                order.append(w)
                stack.append(w)
    inc = {u: 1 for u in nodes}
    exc = {u: 0 for u in nodes}
    for u in reversed(order):
        p = parent[u]
        if p >= 0:
            inc[p] += exc[u]
            exc[p] += max(inc[u], exc[u])
    root = nodes[0]
    return max(inc[root], exc[root])


# --------------------------------------------------------------------------- #
# a(N) = max independence number over subcubic trees on N vertices
# --------------------------------------------------------------------------- #
def _pareto_merge(frontier_a: Iterable[Tuple[int, int]],
                  frontier_b: Iterable[Tuple[int, int]]) -> Set[Tuple[int, int]]:
    """Pareto frontier of the union of sums of two sets of (out, in) pairs."""
    out: Set[Tuple[int, int]] = set()
    for (o1, i1) in frontier_a:
        for (o2, i2) in frontier_b:
            out.add((o1 + o2, i1 + i2))
    # keep only non-dominated pairs (max out, max in)
    result: Set[Tuple[int, int]] = set()
    for p in out:
        if not any(q[0] >= p[0] and q[1] >= p[1] and q != p for q in out):
            result.add(p)
    return result


@functools.lru_cache(maxsize=None)
def alpha_max_subcubic(N: int) -> int:
    """Exact value of ``a(N)`` for ``N >= 0``.

    DP over rooted subcubic trees.  A non-root vertex has at most two children
    and the root has at most three.  For a subtree rooted at ``r`` on ``k``
    vertices write

    * ``in(k, s)``  = best independence number with ``r`` **included**,
    * ``out(k, s)`` = best independence number with ``r`` **excluded**,

    where ``r`` has at most ``s`` children.  Then

    * ``in(k, s)  = 1 + sum over children of out(k_i, 2)``   (children forced out),
    * ``out(k, s) = sum over children of max(in(k_i,2), out(k_i,2))``.

    Both are maxima over compositions of ``k-1`` into ``s`` parts, and since the
    per-child contributions are additive and independent, each maximum is the
    sum of the per-child maxima (a plain max-plus convolution, no Pareto
    bookkeeping required).
    """
    if N <= 0:
        return 0
    NEG = float("-inf")
    # knapsack tables: K_F[k][s] = max sum over a composition of k into s parts of F(part)
    memo: Dict[Tuple[str, int, int], int] = {}

    def knap(which: str, k: int, s: int) -> float:
        """Max over splits of ``k`` vertices into at most ``s`` child subtrees."""
        if k < 0:
            return NEG
        if k == 0:
            return 0.0
        if s == 0:
            return NEG
        key = (which, k, s)
        if key in memo:
            return memo[key]
        best = NEG
        for t in range(0, k + 1):
            f = child_value(which, t)
            rest = knap(which, k - t, s - 1)
            if rest == NEG:
                continue
            best = max(best, f + rest)
        memo[key] = best
        return best

    def child_value(which: str, t: int) -> float:
        """Contribution of a child subtree of size ``t`` (its root has <=2 children).

        ``t == 0`` means "no child in this slot".
        """
        if t == 0:
            return 0.0
        o = subtree_out(t)
        if which == "any":
            return max(o, subtree_in(t))
        return o

    def subtree_in(k: int) -> float:
        if k <= 0:
            return NEG
        if k == 1:
            return 1.0
        return 1.0 + knap("out", k - 1, 2)

    def subtree_out(k: int) -> float:
        if k <= 0:
            return 0.0
        if k == 1:
            return 0.0
        return knap("any", k - 1, 2)

    root_in = 1.0 + knap("out", N - 1, 3)
    root_out = knap("any", N - 1, 3)
    return int(round(max(root_in, root_out)))


def _compositions(total: int, k: int) -> Iterable[Tuple[int, ...]]:
    """All k-tuples of non-negative integers summing to ``total``."""
    if k == 1:
        if total >= 0:
            yield (total,)
        return
    if total < 0:
        return
    for first in range(total + 1):
        for rest in _compositions(total - first, k - 1):
            yield (first,) + rest


def max_subcubic_tree(N: int) -> Tuple[List[Tuple[int, int]], int]:
    """An explicit subcubic tree on ``N`` vertices attaining ``a(N)``.

    Construction: walk along a spine and attach one or two leaves to every
    spine vertex, keeping the total degree at most 3.  The number of leaves is
    maximised, which is what drives the independence number.
    """
    if N <= 0:
        return [], 0
    edges: List[Tuple[int, int]] = []
    if N == 1:
        return edges, 1
    nxt = 1
    spine = [0]
    prev: int | None = None
    while nxt < N:
        # decide how many leaves the current spine vertex may carry
        #   degree so far = 1 (to previous) + (0 or 1 for the next spine link)
        if prev is not None:
            edges.append((prev, nxt))
        prev = nxt
        spine.append(nxt)
        nxt += 1
        # leaves
        if nxt < N:
            # one leaf, keeping room for the spine to continue
            edges.append((spine[-1], nxt))
            nxt += 1
            if nxt < N:
                edges.append((spine[-1], nxt))
                nxt += 1
    val = alpha_tree(_adj(N, edges))
    return edges, val


def _adj(N: int, edges: Sequence[Tuple[int, int]]) -> Dict[int, List[int]]:
    adj: Dict[int, List[int]] = {i: [] for i in range(N)}
    for a, b in edges:
        adj[a].append(b)
        adj[b].append(a)
    return adj


# --------------------------------------------------------------------------- #
# Weak duals
# --------------------------------------------------------------------------- #
def weak_dual(diagonals: Sequence[tuple[int, int]], n: int) -> Tuple[List[Tuple[int, int, int]], List[List[int]]]:
    """Weak dual of the triangulation of the convex ``n``-gon.

    Returns ``(faces, adjacency)`` where ``faces`` lists the triangles as index
    triples and ``adjacency[i]`` lists the indices of faces sharing an edge
    with face ``i``.
    """
    from .convex import triangulation_edges, triangulation_triangle_count

    faces = triangulation_triangle_count(diagonals, n)
    edge_faces: Dict[Tuple[int, int], List[int]] = {}
    for idx, (a, b, c) in enumerate(faces):
        for e in ((a, b), (a, c), (b, c)):
            edge_faces.setdefault(tuple(sorted(e)), []).append(idx)
    adjacency: List[List[int]] = [[] for _ in faces]
    for key, owners in edge_faces.items():
        if len(owners) == 2:
            u, v = owners
            adjacency[u].append(v)
            adjacency[v].append(u)
    return faces, adjacency


class IndependenceNumber:
    """Convenience wrapper computing ``alpha(T)`` for explicit trees."""

    @staticmethod
    def of_edges(N: int, edges: Sequence[tuple[int, int]]) -> int:
        return alpha_tree(_adj(N, edges))