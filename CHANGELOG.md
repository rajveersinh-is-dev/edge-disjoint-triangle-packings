# Changelog

All notable changes to this project are documented here.

## [1.0.0] - 2026-10-04

Initial release.  This repository replaces an earlier project on edge-disjoint
triangulations of a convex polygon, which was discarded as a rediscovery of a
counting bound (see `docs/RESEARCH-NOTES.md`, section 1).

### Mathematics

- **Theorem 4.1.** For `n` points in convex position the largest number of pairwise
  non-crossing, edge-disjoint triangles is `floor((2n-3)/3)`.  Upper bound by edge count
  (Lemma 3.1a); lower bound by a two-parameter recurrence (Proposition 4.2) plus a
  floor-arithmetic induction (Lemma 4.3), with an explicit constructive version.
- **Theorem 5.1.** The convex extremal function is exactly the maximum independence number
  of a subcubic tree: `f(n) = a(n-2) = floor((2n-3)/3)`, where
  `a(N) = floor((2N+1)/3)` (Lemma 5.2: a vertex-cover count plus a caterpillar
  construction).
- **Lemma 3.4.** Every subcubic tree on `N` vertices is the weak dual of a triangulation of
  the convex `(N+2)`-gon.
- **Lemma 3.5.** Every `3`-cycle of a convex-polygon triangulation bounds a face, so the
  maximum packing inside a triangulation equals the independence number of its weak dual.
- **Proposition 6.2 / Lemma 6.3.** An Eulerian maximal planar graph packs `n-2` triangles
  and they partition its edges; Eulerian triangulations exist exactly for
  `n = 3` or `n >= 6` with `n != 7`.
- **Theorem 6.5.** `nu(P) <= n-2` for every `n`-point set in general position, with the
  equality cases classified.  The converse direction is verified exhaustively for `n <= 7`.
- **Conjecture 7.1 (unproved).** `nu(P) >= floor((2n-3)/3)` for every `n`-point set, with
  equality only in convex position.  Reported with a search for counterexamples and a
  discussion of why the obvious induction fails.

### Software

- `src/planetri/geom.py`: exact orientation and crossing predicates on exact points
  (decimal inputs are converted to exact integers, so `0.1 + 0.2 = 0.3` holds).
- `src/planetri/solver.py`: exact maximiser for arbitrary point sets, as maximum independent
  set in the triangle conflict graph, solved by 0/1 programming with HiGHS.
- `src/planetri/convex.py`: the recurrence, the closed form, and constructive packings
  (`convex_optimal_packing`, `convex_packing_avoids_edge`).
- `src/planetri/trees.py`: independence numbers of subcubic trees (`alpha_max_subcubic`,
  exact) and weak duals.
- 47 pytest tests aimed at catching wrong mathematics rather than covering code.

### Experiments

- `exp1_convex.py`: `f(n)` by dynamic program and by exhaustive triangulation enumeration
  with an exact ILP per triangulation; agreement for `n <= 10`.
- `exp2_falsify.py`: eleven structured configuration families, exact `nu(P)` for `n <= 11`,
  explicit search for counterexamples to Conjecture 7.1.
- `exp3_weakdual.py`: exhaustive verification of Lemmas 3.4 and 3.5.
- `exp4_general.py`: all labelled maximal planar graphs on `n <= 7`, Eulerian
  constructions for `n <= 24` with straight-line verified packings, and random-set checks
  of `nu(P) <= n-2`.
- `make_figures.py`, `run_all.py`.

### Documentation

- `paper/main.tex` (compiles with `tectonic`), `README.md`, `docs/RESEARCH-NOTES.md`
  (candidate survey, novelty audit, dead ends, and the bugs the falsification caught).