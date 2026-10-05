"""Exact geometric predicates for planar point sets.

Every routine works on exact integer coordinates.  The :class:`PointSet`
constructor accepts arbitrary numeric input (including floats and numpy scalars)
and converts it to exact integers by going through the shortest decimal
representation of each coordinate, so that e.g. ``0.1 + 0.2`` is exactly
``0.3`` rather than ``0.30000000000000004``.  Callers that care about exactness
should pass integers, ``Fraction`` or ``Decimal``.

The predicates provided are exactly the ones the mathematics needs:

* :func:`orient` --- the sign of the orientation determinant of three points;
* :func:`segments_cross` --- whether two closed segments meet *away from* their
  endpoints.  A shared endpoint is a touch, not a crossing, and is reported as
  ``False``: two triangles of a packing legitimately meet at a vertex, and
  treating that as a conflict collapses the whole problem;
* :func:`convex_hull_order` / :func:`is_in_convex_position`.
"""

from __future__ import annotations

from typing import Iterable, Sequence
import decimal
import fractions
import functools


Coord = tuple[int, int]

__all__ = [
    "Point",
    "PointSet",
    "cross",
    "orient",
    "segments_cross",
    "is_in_convex_position",
    "convex_hull_order",
]


@functools.lru_cache(maxsize=None)
def _to_exact_int(value) -> tuple[int, int]:
    """Return ``value`` as an exact pair ``(numerator, denominator)``.

    Floats and numpy scalars are routed through their shortest decimal
    representation, so decimal input keeps its decimal value exactly.
    """
    if isinstance(value, fractions.Fraction):
        return value.numerator, value.denominator
    if isinstance(value, decimal.Decimal):
        sign, digits, exp = value.as_tuple()
        mantissa = int("".join(map(str, digits)))
        if sign:
            mantissa = -mantissa
        if exp >= 0:
            return mantissa * 10**exp, 1
        return mantissa, 10 ** (-exp)
    if isinstance(value, int):
        return value, 1
    # numpy scalars and anything else exposing .item()
    item = getattr(value, "item", None)
    if callable(item):
        return _to_exact_int(item())
    if isinstance(value, float):
        sign, digits, exp = decimal.Decimal(repr(value)).as_tuple()
        mantissa = int("".join(map(str, digits)))
        if sign:
            mantissa = -mantissa
        if exp >= 0:
            return mantissa * 10**exp, 1
        return mantissa, 10 ** (-exp)
    raise TypeError(f"cannot interpret {value!r} as an exact coordinate")


def _gcd(a: int, b: int) -> int:
    """Gcd.
    
    Args:
        a:
        b:
    
    Returns:
        The computed result
    
    """
    while b:
        a, b = b, a % b
    return a


def _common_scale(raw: Sequence[tuple[tuple[int, int], tuple[int, int]]]) -> int:
    """Smallest positive integer clearing every denominator."""
    scale = 1
    for num_den in raw:
        for _, den in num_den:
            scale = scale * den // _gcd(scale, den)
    return scale


class Point(tuple):
    """An exact planar point with integer coordinates."""

    __slots__ = ()

    def __new__(cls, x, y):
        """New.
        
        Args:
            x:
            y:
        
        Returns:
            The computed result
        
        """
        rx, dx = _to_exact_int(x)
        ry, dy = _to_exact_int(y)
        return super().__new__(cls, (rx * dy, ry * dx))

    @property
    def x(self) -> int:
        """X.
        
        Returns:
            The computed result
        
        """
        return self[0]

    @property
    def y(self) -> int:
        """Y.
        
        Returns:
            The computed result
        
        """
        return self[1]

    def as_float(self) -> tuple[float, float]:
        """As float.
        
        Returns:
            tuple: Result of type tuple
        
        """
        return (float(self[0]), float(self[1]))

    def __repr__(self) -> str:  # pragma: no cover - debugging helper
        """Repr.
        
        Returns:
            The computed result
        
        """
        return f"Point({self[0]}, {self[1]})"


def cross(o: Point, a: Point, b: Point) -> int:
    """Twice the signed area of triangle ``(o, a, b)``, exactly."""
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def orient(o: Point, a: Point, b: Point) -> int:
    """Sign of the orientation determinant of ``(o, a, b)``: -1, 0 or +1."""
    d = cross(o, a, b)
    return (d > 0) - (d < 0)


def segments_cross(p: Point, q: Point, r: Point, s: Point) -> bool:
    """Do the closed segments ``pq`` and ``rs`` meet away from their endpoints?

    A single shared endpoint is a *touch*, not a crossing, and is reported as
    ``False``.  Under the general-position hypothesis (no three of the four
    points collinear) the test is the classical orientation test.  For robustness
    the function also reports ``True`` when the two segments are collinear and
    overlap in a segment of positive length.
    """
    if len({p, q, r, s}) < 4:
        # A common endpoint: the segments touch but do not cross.
        return False
    d1 = cross(p, q, r)
    d2 = cross(p, q, s)
    d3 = cross(r, s, p)
    d4 = cross(r, s, q)
    if d1 == 0 and d2 == 0:
        lo1, hi1 = min(p, q), max(p, q)
        lo2, hi2 = min(r, s), max(r, s)
        lo, hi = max(lo1, lo2), min(hi1, hi2)
        return lo < hi
    return (
        d1 != 0
        and d2 != 0
        and d3 != 0
        and d4 != 0
        and ((d1 > 0) != (d2 > 0))
        and ((d3 > 0) != (d4 > 0))
    )


class PointSet:
    """A finite set of exact points, indexed by position."""

    __slots__ = ("points", "n")

    def __init__(self, coords: Iterable[Sequence]):
        """Init.
        
        Args:
            coords:
        
        """
        raw = [(_to_exact_int(c[0]), _to_exact_int(c[1])) for c in coords]
        scale = _common_scale(raw)
        pts = tuple(
            Point(xn * scale // xd, yn * scale // yd)
            for (xn, xd), (yn, yd) in raw
        )
        seen: set[Coord] = set()
        uniq: list[Point] = []
        for p in pts:
            if p not in seen:
                seen.add(p)
                uniq.append(p)
        self.points: tuple[Point, ...] = tuple(uniq)
        self.n = len(self.points)

    def __len__(self) -> int:
        """Len.
        
        Returns:
            The computed result
        
        """
        return self.n

    def __getitem__(self, i: int) -> Point:
        """Getitem.
        
        Args:
            i:
        
        Returns:
            The computed result
        
        """
        return self.points[i]

    def __iter__(self):
        """Iter.
        
        Returns:
            The computed result
        
        """
        return iter(self.points)

    def __repr__(self) -> str:  # pragma: no cover - debugging helper
        """Repr.
        
        Returns:
            The computed result
        
        """
        return f"PointSet(n={self.n}, {self.points})"

    def crossings(self) -> set[tuple[tuple[int, int], tuple[int, int]]]:
        """All unordered pairs of vertex pairs whose segments cross."""
        n = self.n
        chords = [((i, j), self.points[i], self.points[j])
                  for i in range(n) for j in range(i + 1, n)]
        out: set[tuple[tuple[int, int], tuple[int, int]]] = set()
        for a, item in enumerate(chords):
            (i, j), pi, pj = item
            for b in range(a + 1, len(chords)):
                (k, l), pk, pl = chords[b]
                if len({i, j, k, l}) < 4:
                    continue
                if segments_cross(pi, pj, pk, pl):
                    key = ((i, j), (k, l)) if (i, j) < (k, l) else ((k, l), (i, j))
                    out.add(key)
        return out


def is_in_convex_position(ps: PointSet) -> bool:
    """Are all points of ``ps`` vertices of their convex hull?"""
    return len(convex_hull_order(ps)) == ps.n


def convex_hull_order(ps: PointSet) -> list[int]:
    """Indices of the points on the convex hull, counter-clockwise."""
    pts = list(ps.points)
    n = len(pts)
    if n == 0:
        return []
    order = sorted(range(n), key=lambda i: (pts[i][0], pts[i][1]))
    lower: list[int] = []
    for i in order:
        while len(lower) >= 2 and orient(pts[lower[-2]], pts[lower[-1]], pts[i]) <= 0:
            lower.pop()
        lower.append(i)
    upper: list[int] = []
    for i in reversed(order):
        while len(upper) >= 2 and orient(pts[upper[-2]], pts[upper[-1]], pts[i]) <= 0:
            upper.pop()
        upper.append(i)
    hull = lower[:-1] + upper[:-1]
    if hull:
        m = min(range(len(hull)), key=lambda t: hull[t])
        hull = hull[m:] + hull[:m]
    return hull