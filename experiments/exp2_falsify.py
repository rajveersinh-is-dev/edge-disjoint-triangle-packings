"""Experiment 2 (falsification): can a general position set do *worse* than
convex position?

Claim under test (Conjecture C1):
    for every n-point set P in general position,
        nu(P)  >=  floor((2n - 3) / 3),
    the value attained by every set in convex position.

Strategy: generate many families of configurations -- convex, nearly convex,
convex with one point pushed slightly inside, doubly-nested convex chains,
grids perturbed, "stacked"/Apollonian, random uniform, random clustered --
compute nu(P) exactly by integer programming, and report the minimum observed
per n.  Any configuration with nu(P) < floor((2n-3)/3) refutes C1.
"""

from __future__ import annotations

import csv
import math
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

import numpy as np

from planetri.geom import PointSet, is_in_convex_position
from planetri.solver import is_valid_packing, nu_exact, optimal_packing

RESULTS = os.path.join(os.path.dirname(__file__), "results")
N_MAX = 11


# --------------------------------------------------------------------------- #
# configuration families
# --------------------------------------------------------------------------- #


def convex_ring(n: int, r: float = 1000.0) -> np.ndarray:
    """Convex ring.
    
    Args:
        n:
        r (float):
    
    Returns:
        The computed result
    
    """
    ang = 2 * math.pi * np.arange(n) / n
    return np.stack([r * np.cos(ang), r * np.sin(ang)], axis=1)


def convex_ring_perturbed(n: int, eps: float, seed: int) -> np.ndarray:
    """Convex ring perturbed.
    
    Args:
        n:
        eps:
        seed:
    
    Returns:
        The computed result
    
    """
    rng = np.random.default_rng(seed)
    return convex_ring(n) + eps * rng.normal(size=(n, 2))


def convex_with_interior(n: int, seed: int, depth: float = 0.0) -> np.ndarray:
    """n-1 points on a ring plus one point at (1-depth)*centre."""
    rng = np.random.default_rng(seed)
    ring = convex_ring(n - 1)
    centre = np.array([[0.0, 0.0]]) if depth > 0 else np.array(
        [(1.0 - depth) * ring.mean(axis=0)]
    )
    return np.vstack([ring, centre])


def double_chain(n: int, seed: int, gap: float = 0.02) -> np.ndarray:
    """Two chains, one above the other: n/2 convex 'cups' back to back."""
    k = n // 2
    m = n - k
    xs = np.arange(k) / k
    up = np.stack([xs, gap + 1e-3 * np.sin(3 * np.pi * xs)], axis=1)
    xs2 = np.arange(m) / m
    down = np.stack([xs2, -gap - 1e-3 * np.cos(3 * np.pi * xs2)], axis=1)
    return np.vstack([up, down])


def clustered(n: int, seed: int, k: int = 3) -> np.ndarray:
    """Clustered.
    
    Args:
        n:
        seed:
        k (int):
    
    Returns:
        The computed result
    
    """
    rng = np.random.default_rng(seed)
    centres = rng.normal(size=(k, 2)) * 30
    idx = rng.integers(0, k, size=n)
    return centres[idx] + rng.normal(size=(n, 2))


def grid(n: int, seed: int) -> np.ndarray:
    """Grid.
    
    Args:
        n:
        seed:
    
    Returns:
        The computed result
    
    """
    r = int(math.ceil(math.sqrt(n)))
    pts = [(i, j) for i in range(r) for j in range(r)][:n]
    return np.array(pts, dtype=float) + 0.137


def stacked(n: int, seed: int) -> np.ndarray:
    """Points produced by repeatedly inserting into the largest existing face
    (an Apollonian-type configuration), which yields hull size 3."""
    rng = np.random.default_rng(seed)
    pts = [np.array([0.0, 0.0]), np.array([1000.0, 0.0]), np.array([0.0, 1000.0])]

    def faces(pts):
        """Faces.
        
        Args:
            pts (list):
        
        Returns:
            The computed result
        
        """
        from itertools import combinations

        out = []
        for i, j, k in combinations(range(len(pts)), 3):
            a, b, c = pts[i], pts[j], pts[k]
            d = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
            if d > 0:
                out.append((i, j, k, d))
        return out

    while len(pts) < n:
        fl = faces(pts)
        # pick the largest face with deterministic tie-breaking noise
        best = max(range(len(fl)), key=lambda t: fl[t][3] * (1 + 1e-3 * rng.random()))
        i, j, k, _ = fl[best]
        a, b, c = pts[i], pts[j], pts[k]
        g = (a + b + c) / 3.0 + rng.normal(size=2) * 3.0
        pts.append(g)
    return np.array(pts)


def random_uniform(n: int, seed: int) -> np.ndarray:
    """Random uniform.
    
    Args:
        n:
        seed:
    
    Returns:
        The computed result
    
    """
    rng = np.random.default_rng(seed)
    return rng.normal(size=(n, 2)) * 1000


def almost_collinear(n: int, seed: int, curve: float) -> np.ndarray:
    """Almost collinear.
    
    Args:
        n:
        seed:
        curve:
    
    Returns:
        The computed result
    
    """
    xs = np.linspace(0, 1000, n)
    return np.stack([xs, curve * (xs / 1000.0) ** 2], axis=1)


FAMILIES = {
    "convex": lambda n, s: convex_ring(n),
    "convex_perturbed": lambda n, s: convex_ring_perturbed(n, 1.0, s),
    "convex_interior_deep": lambda n, s: convex_with_interior(n, s, 0.0),
    "convex_interior_mid": lambda n, s: convex_with_interior(n, s, 0.5),
    "double_chain": double_chain,
    "clustered": clustered,
    "grid": grid,
    "stacked": stacked,
    "random_uniform": random_uniform,
    "parabola_flat": lambda n, s: almost_collinear(n, s, 1.0),
    "parabola_steep": lambda n, s: almost_collinear(n, s, 400.0),
}


def general_position_ok(ps: PointSet) -> bool:
    """No three points collinear."""
    from planetri.geom import orient

    n = ps.n
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                if orient(ps[i], ps[j], ps[k]) == 0:
                    return False
    return True


def main() -> int:
    """Entry point — parse arguments and run the main computation.
    
    Returns:
        int: Result of type int
    
    """
    os.makedirs(RESULTS, exist_ok=True)
    rows = []
    violations = []
    n0 = time.time()
    for n in range(3, N_MAX + 1):
        target = (2 * n - 3) // 3
        best_min = math.inf
        for fam, fn in FAMILIES.items():
            seeds = range(6) if fam not in ("convex",) else range(1)
            for s in seeds:
                coords = fn(n, s)
                # round to 4 decimals so that all arithmetic is exact
                coords = np.round(coords, 4)
                ps = PointSet([tuple(map(float, row)) for row in coords])
                if ps.n != n:
                    continue
                if not general_position_ok(ps):
                    continue
                val = nu_exact(ps)
                convex = is_in_convex_position(ps)
                rows.append(
                    {
                        "n": n,
                        "family": fam,
                        "seed": s,
                        "nu": val,
                        "target_convex": target,
                        "in_convex_position": int(convex),
                        "margin": val - target,
                    }
                )
                best_min = min(best_min, val)
                if val < target:
                    violations.append((n, fam, s, val, target))
        print(f"n={n:2d}  target(convex)={target:2d}  min nu observed={best_min:2d}"
              f"  {'OK' if best_min >= target else 'COUNTEREXAMPLE'}", flush=True)

    path = os.path.join(RESULTS, "falsification_general_position.csv")
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    print()
    print(f"configurations tested : {len(rows)}")
    print(f"elapsed               : {time.time() - n0:.1f}s")
    if violations:
        print("CONJECTURE C1 REFUTED. counterexamples:")
        for v in violations[:
            40]:
            print("   n=%d family=%s seed=%d  nu=%d < %d" % v)
        return 1
    print("Conjecture C1 survives all tested configurations "
          "(no set with nu(P) < floor((2n-3)/3)).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())