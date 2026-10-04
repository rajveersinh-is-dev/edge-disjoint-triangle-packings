"""planetri: edge-disjoint, non-crossing triangle packings of planar point sets.

The central object of study is

    nu(P) = the largest size of a family of triangles spanned by a finite point
            set P in general position in the plane such that no two triangles
            share an edge and no edge of one triangle crosses an edge of
            another.

The package provides

* exact predicates for convex position (:mod:`planetri.geom`),
* an exact maximiser for arbitrary point sets (:mod:`planetri.solver`),
* an exact recurrence for the convex-position extremal function
  (:mod:`planetri.convex`),
* weak-dual tree combinatorics (:mod:`planetri.trees`).
"""

from planetri.geom import Point, PointSet, cross, orient, segments_cross
from planetri.solver import (
    Triangle,
    all_triangles,
    conflicting,
    is_valid_packing,
    nu_exact,
    optimal_packing,
)
from planetri.convex import (
    convex_nu_dp,
    convex_nu_table,
    convex_nu_closed_form,
    enumerate_triangulations,
    triangulation_triangle_count,
)
from planetri.trees import (
    IndependenceNumber,
    alpha_max_subcubic,
    weak_dual,
)

__all__ = [
    "Point",
    "PointSet",
    "cross",
    "orient",
    "segments_cross",
    "Triangle",
    "all_triangles",
    "conflicting",
    "is_valid_packing",
    "nu_exact",
    "optimal_packing",
    "convex_nu_dp",
    "convex_nu_table",
    "convex_nu_closed_form",
    "enumerate_triangulations",
    "triangulation_triangle_count",
    "IndependenceNumber",
    "alpha_max_subcubic",
    "weak_dual",
]

__version__ = "1.0.0"