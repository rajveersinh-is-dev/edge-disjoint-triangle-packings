"""Experiment 4: the general-position side.

Claim under test (Theorem C):
    (i)   nu(P) <= n - 2 for every n-point set P in general position;
    (ii)  max over n-point sets P of nu(P) equals n-2 exactly when an
          *Eulerian* maximal planar graph on n vertices exists, i.e. for
          n in {3} u {n >= 6 : n != 7};
    (iii) the bridge: nu(P) = n-2 forces the union of the packing to be a
          maximal planar graph all of whose edges lie in the chosen triangles,
          and a maximal planar graph whose edges admit a partition into
          triangles is Eulerian (its cubic dual is bipartite).

Evidence collected:
  * exhaustive over *all* maximal planar graphs on n = 3..7 vertices:
    nu(M) by ILP, degree parity, and whether E(M) admits an exact triangle
    partition;
  * explicit Eulerian constructions for n <= 24, drawn with straight lines via
    networkx's planar layout and partitioned into triangles (verified
    geometrically, not just combinatorially);
  * random n-point sets for n <= 13, checking nu(P) <= n-2.
"""

from __future__ import annotations

import csv
import functools
import itertools
import os
import sys
import time

from planetri.geom import PointSet, orient, segments_cross
from planetri.solver import nu_exact
from scipy.optimize import Bounds, LinearConstraint, milp
import networkx as nx
import numpy as np


sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))



RESULTS = os.path.join(os.path.dirname(__file__), "results")
# The exhaustive scan of *all* maximal planar graphs is exponential; n = 7
# already takes a few minutes.  Override with PLANETRI_MAX_ENUM_N.
MAX_ENUM_N = int(os.environ.get("PLANETRI_MAX_ENUM_N", "7"))


# --------------------------------------------------------------------------- #
# exact triangle packing on an abstract graph
# --------------------------------------------------------------------------- #


def triangle_list(n: int, edges):
    """Triangle list.
    
    Args:
        n:
        edges:
    
    Returns:
        The computed result
    
    """
    es = {tuple(sorted(e)) for e in edges}
    return [
        t
        for t in itertools.combinations(range(n), 3)
        if (t[0], t[1]) in es and (t[0], t[2]) in es and (t[1], t[2]) in es
    ]


def max_edge_disjoint_triangles_of_graph(n: int, edges):
    """ILP: max # pairwise edge-disjoint triangles of G, with a witness."""
    tris = triangle_list(n, edges)
    m = len(tris)
    if m == 0:
        return 0, []
    owners: dict[tuple[int, int], list[int]] = {}
    for idx, t in enumerate(tris):
        for e in ((t[0], t[1]), (t[0], t[2]), (t[1], t[2])):
            owners.setdefault(e, []).append(idx)
    conflict = set()
    for lst in owners.values():
        for i, item in enumerate(lst):
            for j in range(i + 1, len(lst)):
                conflict.add((min(item, lst[j]), max(item, lst[j])))
    A = np.zeros((len(conflict), m))
    for row, (a, b) in enumerate(sorted(conflict)):
        A[row, a] = 1.0
        A[row, b] = 1.0
    cons = LinearConstraint(A, -np.inf, np.ones(len(conflict)))
    res = milp(c=-np.ones(m), constraints=cons, integrality=np.ones(m), bounds=Bounds(0, 1))
    x = np.round(res.x).astype(bool)
    return int(round(-res.fun)), [tris[i] for i in range(m) if x[i]]


def partition_edges_into_triangles(n: int, edges):
    """Exact cover of E(G) by triangles, or None if impossible."""
    es = {tuple(sorted(e)) for e in edges}
    tris = triangle_list(n, edges)
    trisedges = [
        {tuple(sorted(x)) for x in ((t[0], t[1]), (t[0], t[2]), (t[1], t[2]))} for t in tris
    ]

    @functools.lru_cache(maxsize=None)
    def cover(rem, used_tris):
        """Cover.
        
        Args:
            rem:
            used_tris:
        
        Returns:
            The computed result
        
        """
        if not rem:
            return used_tris
        e = min(rem)
        for i, te in enumerate(trisedges):
            if e in te and te <= rem:
                res = cover(rem - te, used_tris + [tris[i]])
                if res is not None:
                    return res
        return None

    return cover(es, [])


def all_maximal_planar_graphs(n: int):
    """Every labelled maximal planar graph on n vertices, by brute force."""
    total = 3 * n - 6
    allpairs = list(itertools.combinations(range(n), 2))
    for combo in itertools.combinations(allpairs, total):
        G = nx.Graph()
        G.add_nodes_from(range(n))
        G.add_edges_from(combo)
        if nx.check_planarity(G)[0]:
            yield G


def is_eulerian(G: nx.Graph) -> bool:
    """Is eulerian.
    
    Args:
        G:
    
    Returns:
        The computed result
    
    """
    return all(d % 2 == 0 for _, d in G.degree())


# --------------------------------------------------------------------------- #
# explicit Eulerian constructions
# --------------------------------------------------------------------------- #


def bipyramid_over_even_cycle(k: int) -> nx.Graph:
    """Cycle on 0..k-1 plus two apexes k, k+1 joined to every cycle vertex.

    Eulerian exactly when k is even (apex degree k, cycle vertices degree 4).
    """
    G = nx.Graph()
    G.add_nodes_from(range(k + 2))
    G.add_edges_from((i, (i + 1) % k) for i in range(k))
    for i in range(k):
        G.add_edge(k, i)
        G.add_edge(k + 1, i)
    return G


def octahedron() -> nx.Graph:
    """The octahedron K_{2,2,2} on 6 vertices; all degrees 4, so Eulerian."""
    G = nx.Graph()
    G.add_nodes_from(range(6))
    for u in range(6):
        for v in range(u + 1, 6):
            same = (u < 2 and v < 2) or (2 <= u < 4 and 2 <= v < 4) or (u >= 4 and v >= 4)
            if not same:
                G.add_edge(u, v)
    return G


def glue_along_face(G1: nx.Graph, f1, G2: nx.Graph, f2) -> nx.Graph:
    """Connected sum of two sphere triangulations along the faces ``f1``, ``f2``.

    Delete the open faces, then identify the two triangular boundaries with
    reversed orientation.  The result has ``|V1| + |V2| - 3`` vertices and is a
    triangulation of the sphere.  A vertex lying on the glued face has degree
    ``d1 + d2 - 4``, so Eulerian-ness of both factors is inherited.
    """
    off = G1.number_of_nodes()
    # orientation-reversing identification of the two face boundaries
    ident = {f2[(2 - i) % 3]: f1[i] for i in range(3)}
    fresh: dict[int, int] = {}
    for v in sorted(G2.nodes()):
        if v not in ident:
            fresh[v] = off + len(fresh)
    H = nx.Graph()
    H.add_nodes_from(range(off + len(fresh)))
    for a, b in G1.edges():
        H.add_edge(a, b)
    for a, b in G2.edges():
        H.add_edge(ident[a] if a in ident else fresh[a], ident[b] if b in ident else fresh[b])
    if not nx.check_planarity(H)[0]:
        raise AssertionError("connected sum is not planar")
    return H


OCTAHEDRON_FACE = (0, 2, 4)


def bipyramid_face(k: int) -> tuple[int, int, int]:
    """A face of the bipyramid over the cycle on 0..k-1 (apex k)."""
    return (k, 0, 1)


def eulerian_triangulation(n: int):
    """Explicit maximal planar Eulerian graph on n vertices, or None.

    * even n >= 6: the bipyramid over the even cycle C_{n-2};
    * odd n >= 9:   the connected sum of the octahedron (6 vertices) with the
      bipyramid over C_{n-5} (n-3 vertices), which has 6 + (n-3) - 3 = n
      vertices and even degrees throughout;
    * n = 3:        the single triangle.

    No Eulerian triangulation exists for n in {4, 5, 7}.
    """
    if n == 3:
        G = nx.Graph()
        G.add_edges_from([(0, 1), (1, 2), (0, 2)])
        return G
    if n in (4, 5, 7):
        return None
    if n % 2 == 0:
        return bipyramid_over_even_cycle(n - 2)
    k = n - 5  # even and >= 4
    bip = bipyramid_over_even_cycle(k)
    return glue_along_face(octahedron(), OCTAHEDRON_FACE, bip, bipyramid_face(k))


def geometric_check(G: nx.Graph, tris) -> bool:
    """Draw G with straight lines; verify the chosen triangles do not cross."""
    pos = nx.planar_layout(G)
    coords = [(round(1000 * pos[v][0], 4), round(1000 * pos[v][1], 4)) for v in sorted(G)]
    ps = PointSet(coords)
    for i, item in enumerate(tris):
        for j in range(i + 1, len(tris)):
            a, b = item, tris[j]
            if len(set(a) & set(b)) >= 2:
                return False
            for x, y in ((a[0], a[1]), (a[0], a[2]), (a[1], a[2])):
                for u, v in ((b[0], b[1]), (b[0], b[2]), (b[1], b[2])):
                    if len({x, y, u, v}) < 4:
                        continue
                    if segments_cross(ps[x], ps[y], ps[u], ps[v]):
                        return False
    return True


def general_position(ps: PointSet) -> bool:
    """General position.
    
    Args:
        ps:
    
    Returns:
        bool: Result of type bool
    
    """
    n = ps.n
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                if orient(ps[i], ps[j], ps[k]) == 0:
                    return False
    return True


# --------------------------------------------------------------------------- #


def main() -> int:
    """Entry point — parse arguments and run the main computation.
    
    Returns:
        int: Result of type int
    
    """
    os.makedirs(RESULTS, exist_ok=True)
    problems: list[str] = []

    print("=== Exhaustive scan of all maximal planar graphs ===")
    rows = []
    for n in range(3, MAX_ENUM_N + 1):
        t0 = time.time()
        count = 0
        best_nu = 0
        eulerian_nu: set[int] = set()
        non_eulerian_nu: set[int] = set()
        part_flags: set[tuple[int, int]] = set()
        for G in all_maximal_planar_graphs(n):
            count += 1
            edges = list(G.edges())
            nu, _ = max_edge_disjoint_triangles_of_graph(n, edges)
            eul = int(is_eulerian(G))
            best_nu = max(best_nu, nu)
            (eulerian_nu if eul else non_eulerian_nu).add(nu)
            part = partition_edges_into_triangles(n, edges)
            if part is not None:
                part_flags.add((1, eul))
        print(
            f"n={n}: {count:7d} graphs, max nu = {best_nu} (n-2 = {n-2}); "
            f"nu on Eulerian = {sorted(eulerian_nu)}; "
            f"nu on non-Eulerian = {sorted(non_eulerian_nu)}; "
            f"(partitioned, eulerian) = {sorted(part_flags)};  {time.time() - t0:.1f}s",
            flush=True,
        )
        rows.append(
            {
                "n": n,
                "num_maximal_planar_graphs": count,
                "max_nu": best_nu,
                "n_minus_2": n - 2,
                "nu_values_eulerian": ";".join(map(str, sorted(eulerian_nu))),
                "nu_values_non_eulerian": ";".join(map(str, sorted(non_eulerian_nu))),
                "reaches_n_minus_2": int(best_nu == n - 2),
            }
        )
        if best_nu > n - 2:
            problems.append(f"n={n}: nu={best_nu} exceeds n-2={n - 2}")
        if max(eulerian_nu, default=-1) == n - 2:
            pass
        if max(non_eulerian_nu, default=-1) == n - 2:
            problems.append(f"n={n}: a NON-Eulerian maximal planar graph reaches n-2")

    print()
    print("=== Explicit Eulerian constructions, straight-line verified ===")
    con_rows = []
    for n in range(3, 25):
        G = eulerian_triangulation(n)
        if G is None:
            con_rows.append(
                {"n": n, "constructed": 0, "nu": "", "eulerian": "", "straight_line_ok": "",
                 "n_minus_2": n - 2}
            )
            continue
        edges = list(G.edges())
        nu, witness = max_edge_disjoint_triangles_of_graph(n, edges)
        eul = int(is_eulerian(G))
        sl = int(geometric_check(G, witness))
        con_rows.append(
            {"n": n, "constructed": 1, "nu": nu, "eulerian": eul,
             "straight_line_ok": sl, "n_minus_2": n - 2}
        )
        ok = (nu == n - 2) and eul and sl
        if not ok:
            problems.append(f"n={n}: construction gave nu={nu}, eulerian={eul}, sl={sl}")
        print(f"  n={n:2d}  nu = {nu} (n-2 = {n-2}), Eulerian = {eul}, "
              f"straight-line packing verified = {sl}", flush=True)

    print()
    print("=== nu(P) <= n - 2 on random point sets ===")
    rng = np.random.default_rng(20240)
    violations = 0
    tested = 0
    for n in range(3, 14):
        for _ in range(12):
            coords = np.round(rng.normal(size=(n, 2)) * 1000, 4)
            ps = PointSet([tuple(row) for row in coords])
            if ps.n != n or not general_position(ps):
                continue
            v = nu_exact(ps)
            tested += 1
            if v > n - 2:
                violations += 1
                problems.append(f"nu(P)={v} > n-2={n-2}")
    print(f"  {tested} random sets tested, violations = {violations}")

    for name, data in (
        ("maximal_planar_scan.csv", rows),
        ("eulerian_constructions.csv", con_rows),
    ):
        with open(os.path.join(RESULTS, name), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(data[0].keys()))
            w.writeheader()
            w.writerows(data)

    print()
    if problems:
        print("PROBLEMS FOUND:")
        for p in problems:
            print("  -", p)
        return 1
    print("All checks consistent with Theorem C.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())