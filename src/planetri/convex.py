"""The convex-position extremal function and its exact recurrence.

Let ``f(m)`` denote the maximum number of pairwise non-crossing, edge-disjoint
triangles spanned by ``m`` points in convex position.  Every non-crossing edge
set on ``m`` convexly positioned points extends to a triangulation of the
``m``-gon, so ``f(m)`` equals

    max over triangulations T of the convex m-gon of
        (max number of pairwise edge-disjoint 3-cycles of T).

Triangulation structure gives a two-parameter dynamic program.  For a
sub-polygon ``[i, j]`` (vertices ``v_i, ..., v_j`` in cyclic order) with root
chord ``(v_i, v_j)`` let

    P(i, j) = max # edge-disjoint 3-cycles over triangulations of [i, j],
    Q(i, j) = the same, subject to the chord (v_i, v_j) being unused.

The root triangle of a triangulation of ``[i, j]`` is ``(v_i, v_k, v_j)``; the
two sub-triangulations share only the vertex ``v_k``.  Moreover every 3-cycle
of ``T(i, j)`` is either a 3-cycle of ``T(i, k)``, a 3-cycle of ``T(k, j)``, or
the root triangle itself, because the vertices of ``[i, j]`` are in convex
position and therefore no 3-cycle of ``T(i, j)`` can strictly contain another
edge of ``T(i, j)``.  Hence

    P(i, j) = max_{i < k < j} max( P(i,k) + P(k,j),  1 + Q(i,k) + Q(k,j) )
    Q(i, j) = max_{i < k < j} ( P(i,k) + P(k,j) )

(the root triangle uses ``(v_i, v_j)``, so it is excluded in ``Q``).  By cyclic
symmetry both depend only on the number of vertices.
"""

from __future__ import annotations

from typing import Iterator, Sequence
import functools


__all__ = [
    "convex_nu_dp",
    "convex_nu_table",
    "convex_nu_closed_form",
    "enumerate_triangulations",
    "triangulation_edges",
    "triangulation_triangle_count",
    "CATALAN",
]


def CATALAN(k: int) -> int:
    """The Catalan number ``C_k``."""
    num = 1
    den = 1
    for i in range(1, k + 1):
        num = num * 2 * (2 * i - 1)
        den = den * (i + 1)
    return num // den


@functools.lru_cache(maxsize=None)
def convex_nu_dp(m_max: int) -> tuple[int, ...]:
    """Return ``(p_0, p_1, ..., p_{m_max})`` with ``p_m = f(m)``.

    Also computes ``q_m = Q`` internally and asserts the closed form
    ``f(m) = floor((2m - 3) / 3)`` whenever the recurrence agrees with it.
    """
    p = [0] * (m_max + 1)
    q = [0] * (m_max + 1)
    for m in range(2, m_max + 1):
        best_p = 0
        best_q = 0
        # a + b = m + 1 with a, b >= 2
        for a in range(2, m):
            b = m + 1 - a
            best_q = max(best_q, p[a] + p[b])
            best_p = max(best_p, p[a] + p[b], 1 + q[a] + q[b])
        p[m] = best_p
        q[m] = best_q
    return tuple(p)


def convex_nu_table(m_max: int) -> list[tuple[int, int, int, int]]:
    """Rows ``(m, f(m), closed form, agreement)`` for ``2 <= m <= m_max``."""
    tab = convex_nu_dp(m_max)
    return [(m, tab[m], convex_nu_closed_form(m), tab[m] == convex_nu_closed_form(m))
            for m in range(2, m_max + 1)]


def convex_nu_closed_form(m: int) -> int:
    """``floor((2m - 3) / 3)``, the exact closed form of ``f(m)``."""
    return (2 * m - 3) // 3


def convex_nu_dp_tables(m_max: int) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Return ``(p, q)`` for ``0 <= m <= m_max`` where ``p_m = f(m)`` and

    ``q_m`` = the largest packing of a convex ``m``-gon whose *closing edge*
    ``(v_0, v_{m-1})`` is not used by any triangle.

    The two recurrences are

    ``p_m = max_{a+b=m+1} max( p_a + p_b , 1 + q_a + q_b )``
    ``q_m = max_{a+b=m+1} ( p_a + p_b )``

    with ``p_0 = p_1 = q_0 = q_1 = 0`` and ``p_2 = q_2 = 0``.
    """
    p = [0] * (m_max + 1)
    q = [0] * (m_max + 1)
    for m in range(2, m_max + 1):
        best_p = 0
        best_q = 0
        for a in range(2, m):
            b = m + 1 - a
            best_q = max(best_q, p[a] + p[b])
            best_p = max(best_p, p[a] + p[b], 1 + q[a] + q[b])
        p[m] = best_p
        q[m] = best_q
    return tuple(p), tuple(q)


@functools.lru_cache(maxsize=None)
def _build_packing(p, q, m: int, avoid: bool) -> list[tuple[int, int, int]]:
    """Backtrack the recurrence on the sub-polygon ``(v_0, ..., v_m)``.

    ``m`` is the largest local index, so the sub-polygon has ``m + 1`` vertices.
    Splitting at the root triangle ``(v_0, v_a, v_m)`` leaves two sub-polygons
    with ``a`` and ``m + 2 - a`` vertices.  When the root triangle is used, its
    edges ``(v_0, v_{a-1})`` and ``(v_a, v_m)`` become the closing edges of the
    two sub-polygons and must therefore be avoided by them; this is what the
    ``avoid`` flag records.

    If ``avoid`` is set on entry, the closing edge ``(v_0, v_m)`` must not occur
    in the output.
    """
    if m < 2:
        return []
    if m == 2:
        # A 3-gon has a single triangle and it uses all three edges, including
        # the closing edge, so nothing is available when the closing edge is
        # forbidden.
        return [] if avoid else [(0, 1, 2)]
    best_val, best = -1, (2, False)
    for a in range(2, m + 1):
        b = m + 2 - a
        v1 = p[a] + p[b]                      # root triangle not used
        if v1 > best_val:
            best_val, best = v1, (a, False)
        if not avoid:
            v2 = 1 + q[a] + q[b]              # root triangle used
            if v2 > best_val:
                best_val, best = v2, (a, True)
    a, use_root = best
    b = m + 2 - a
    k = a - 1                       # split position: root triangle is (v_0, v_k, v_m)
    left = _build_packing(p, q, k, use_root)
    right = [(x + k, y + k, z + k) for (x, y, z) in _build_packing(p, q, m - k, use_root)]
    root = [(0, k, m)] if use_root else []
    return left + root + right


def convex_optimal_packing(n: int) -> list[tuple[int, int, int]]:
    """An explicit optimal packing of the convex ``n``-gon (indices 0..n-1).

    Backtracks through the recurrence of :func:`convex_nu_dp_tables`, so the
    output has exactly ``f(n) = floor((2n-3)/3)`` triangles, no two sharing an
    edge, and no two having crossing edges.
    """
    p, q = convex_nu_dp_tables(max(n + 4, 8))
    return sorted(_build_packing(p, q, n - 1, False))


def convex_packing_avoids_edge(n: int) -> list[tuple[int, int, int]]:
    """Optimal packing of the convex ``n``-gon avoiding the closing edge
    ``(v_0, v_{n-1})``.

    Realising the ``q`` branch of the recurrence; see :func:`convex_optimal_packing`
    for the coding conventions.
    """
    p, q = convex_nu_dp_tables(max(n + 4, 8))
    return sorted(_build_packing(p, q, n - 1, True))


def triangulation_edges(diagonals: Sequence[tuple[int, int]], n: int) -> set[tuple[int, int]]:
    """Edge set of the triangulation of the ``n``-gon with the given diagonals.

    Every edge is stored as a normalised pair ``(min, max)``.
    """
    edges = {(i, (i + 1) % n) if i < (i + 1) % n else ((i + 1) % n, i) for i in range(n)}
    edges |= {(min(a, b), max(a, b)) for a, b in diagonals}
    return edges


def enumerate_triangulations(n: int) -> Iterator[frozenset[tuple[int, int]]]:
    """Yield every triangulation of the convex ``n``-gon as a set of diagonals.

    The number of triangulations is the Catalan number ``C_{n-2}``, so this is
    feasible by exhaustive enumeration for ``n`` up to about 12.
    """
    if n < 3:
        yield frozenset()
        return

    @functools.lru_cache(maxsize=None)
    def gen(lo: int, hi: int, acc: frozenset) -> Iterator[frozenset]:
        """Triangulate the sub-polygon with vertices ``lo, lo+1, ..., hi``.

        Splitting at ``k`` creates the root triangle ``(lo, k, hi)``; the chords
        ``(lo, k)`` and ``(k, hi)`` are diagonals of the original ``n``-gon
        exactly when their index gap is at least 2.  The top-level closing edge
        ``(0, n-1)`` is never emitted because it is a genuine polygon side.
        """
        if hi - lo < 2:
            yield acc
            return
        for k in range(lo + 1, hi):
            left = acc | (frozenset({(lo, k)}) if (k - lo) >= 2 else frozenset())
            right = left | (frozenset({(k, hi)}) if (hi - k) >= 2 else frozenset())
            for a in gen(lo, k, left):
                for b in gen(k, hi, a | (right - left)):
                    yield b

    yield from gen(0, n - 1, frozenset())


def triangulation_triangle_count(diagonals: Sequence[tuple[int, int]], n: int) -> list[tuple[int, int, int]]:
    """All 3-cycles of the triangulation with the given diagonals."""
    edges = triangulation_edges(diagonals, n)
    out = []
    for a in range(n):
        for b in range(a + 1, n):
            if (a, b) not in edges:
                continue
            for c in range(b + 1, n):
                if (a, c) in edges and (b, c) in edges:
                    out.append((a, b, c))
    return out