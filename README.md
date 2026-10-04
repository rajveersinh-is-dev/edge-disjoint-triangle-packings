# Edge-disjoint non-crossing triangle packings of planar point sets

How many triangles can you draw on a set of `n` points so that **no two share an edge**
and **no edge of one crosses an edge of another**?  The answer depends on the shape of
the point set by a factor of about `3/4`, and this repository determines exactly how.

For a finite set `P` of points in general position let `nu(P)` be the largest size of a
family of triangles spanned by `P` such that no two triangles share an edge and no edge of
one triangle crosses an edge of another.

```
  convex position:   nu(P) = floor((2n - 3) / 3)          [proved exactly]
  best configuration: nu(P) = n - 2                        [proved exactly, n = 3 or n >= 6, n != 7]
  convex hexagons therefore support 3 triangles,
  a straight-line octahedron supports 4.
```

## TL;DR

For `n` points in convex position the answer is `floor((2n-3)/3)`.  More interestingly,
that geometric extremal function is *not really geometric*: it equals the maximum
independence number of a subcubic tree on `n-2` vertices, which is `floor((2N+1)/3)`.  For
arbitrary point sets the sharp bound is `n-2`, and it is attained exactly when an Eulerian
triangulation on `n` vertices exists, i.e. for `n = 3` or `n >= 6` with `n != 7`.

## Research question

> Is convex position the *worst* configuration for packing non-crossing, edge-disjoint
> triangles?  Equivalently: is `floor((2n-3)/3)` a lower bound for `nu(P)` for every
> `n`-point set in general position?

We prove the convex value exactly, prove the sharp upper bound for arbitrary sets, and
classify the extremal configurations.  The lower bound for arbitrary configurations is
stated as **Conjecture 7.1 in the paper**: it survives an extensive search for
counterexamples but we do not prove it, and we say so plainly.

## Main results

| Statement | Result | Status |
|---|---|---|
| Theorem 4.1 | `f(n) = floor((2n-3)/3)` for convex `n`-point sets | proved |
| Theorem 5.1 | `f(n) = a(n-2) = floor((2n-3)/3)`, `a(N) = floor((2N+1)/3)` | proved |
| Lemma 3.4 | every subcubic tree on `N` vertices is the weak dual of a convex `(N+2)`-gon triangulation | proved |
| Lemma 3.5 | every `3`-cycle of a convex-polygon triangulation bounds a face | proved |
| Proposition 6.2 | an Eulerian maximal planar graph packs `n-2` triangles, and they partition its edges | proved |
| Lemma 6.3 | Eulerian triangulations exist iff `n = 3` or `n >= 6`, `n != 7` | proved |
| Theorem 6.5 | `nu(P) <= n-2`; equality cases classified | proved (`<=`), equality direction verified exhaustively for `n <= 7` |
| Conjecture 7.1 | `nu(P) >= floor((2n-3)/3)`, equality only in convex position | **conjecture, unproved** |

## Why this is interesting

The union of a packing is always a plane graph, so a counting argument bounds the answer by
`(edges)/(3)`.  In convex position that plane graph is forced to be *outerplanar*, so the
bound is `(2n-3)/3`; with interior points it may be a full planar graph, giving `(3n-6)/3 =
n-2`.  Convex position is therefore expected to be the *hardest* case, and the results
confirm the endpoints exactly.

The second surprise is structural: two elementary observations --- no vertex of a convex
point set lies inside a triangle of three others, and every subcubic tree is the weak dual
of a polygon triangulation --- convert the convex extremal problem into a two-line
combinatorial one.

## Key theorem

```python
from planetri.convex import convex_optimal_packing, convex_nu_closed_form
from planetri.trees import alpha_max_subcubic

assert convex_nu_closed_form(12) == 7
assert alpha_max_subcubic(10) == 7                    # a(10) = floor((2*10+1)/3)
assert convex_optimal_packing(12)[:3] == [(0, 1, 11), (1, 2, 3), (3, 4, 11)]
```

`convex_optimal_packing` is not a search: it backtracks through the recurrence in
Proposition 4.2 and returns an explicit packing, which the test suite validates with exact
geometric predicates.

## Computational verification

* `f(n)` is computed by two independent routes --- a two-parameter dynamic program and
  exhaustive enumeration of all `C_{n-2}` triangulations, each solved to optimality by 0/1
  programming --- which agree for `n <= 10`, and the dynamic program matches the closed form
  for `n <= 400`.
* Lemma 3.5 ("every `3`-cycle is a face") is verified exhaustively for `n <= 10`, and
  Lemma 3.4 ("every subcubic tree occurs") for `n <= 8`.
* `a(N) = floor((2N+1)/3)` is verified for `N <= 40` and cross-checked against 3000 random
  subcubic trees per `N`.
* Theorem 6.5: **all** `1 + 10 + 195 + 5712` labelled maximal planar graphs on `n = 3..7`
  vertices are enumerated and solved exactly.  For `n = 6` only Eulerian graphs reach
  `n-2`.
* Explicit Eulerian constructions for every `n <= 24`, each drawn with straight lines and
  its packing checked geometrically, not just combinatorially.
* Conjecture 7.1 is attacked with 11 structured configuration families plus the six families
  in the figures; **423 configurations** for `n <= 11` and more for `n <= 12`, none of them
  dipping below the convex value.

## Repository structure

```
src/planetri/          core library
  geom.py              exact orientation and crossing predicates on exact points
  solver.py            conflict graph + exact maximiser via 0/1 programming (HiGHS)
  convex.py            the recurrence, the closed form, constructive packings
  trees.py             independence numbers of subcubic trees, weak duals
tests/                 pytest suite (47 tests) aimed at catching bad mathematics
experiments/           five scripts; results/ holds their committed CSV outputs
  exp1_convex.py         f(n) by DP and by full enumeration
  exp2_falsify.py        Conjecture 7.1: search for counterexamples
  exp3_weakdual.py       structural lemmas, exhaustively
  exp4_general.py        Theorem 6.5, incl. all maximal planar graphs on n<=7
  make_figures.py        all four figures
  run_all.py             one-command reproduction
figures/               generated figures
paper/                 the manuscript (main.tex, compiles with tectonic)
docs/                  research notes: candidate survey, novelty audit, dead ends
```

## Reproducing

Requires Python >= 3.10, `numpy`, `scipy` and `networkx`.

```bash
pip install -e .
python -m pytest tests -q              # 47 tests, ~6 s
python experiments/run_all.py          # everything, ~5 min
python experiments/run_all.py --quick  # skips the n=7 planar scan
```

The paper compiles with `tectonic`:

```bash
cd paper && tectonic -X compile main.tex --keep-intermediates   # run twice
```

All experiments are deterministic; the few random searches use fixed seeds.

## Figures

| Figure | Question it answers |
|---|---|
| `fig_convex_packings.png` | what does an optimal packing look like? |
| `fig_three_functions.png` | do `f(n)`, `q(n)` and `a(n-2)` coincide? |
| `fig_configurations.png` | does any configuration dip below convex position? |
| `fig_extremal_gap.png` | how large is the convex/general gap for `n = 6`? |

## Limitations

These are stated in the paper too, and repeated here because they matter:

1. **Conjecture 7.1 is unproved.**  A finite search over configurations is not a proof.  We
   document in the paper exactly why the obvious induction fails (the arithmetic needs two
   triangles per peeled triple, and two edge-disjoint triangles need five vertices).
2. **One direction of Theorem 6.5 is only machine-verified.**  That `nu = n-2` forces an
   Eulerian triangulation is proved for `n <= 7` by exhaustion; the general implication is
   open here.
3. **The convex value is a counting bound.**  The upper bound is a two-line edge count, so
   the content of the convex part is the construction and the reduction to trees, not the
   value.
4. **No variants.**  Convex `k`-gons instead of triangles, weighted versions and
   overlap-permitting versions are all untouched.
5. **The straight-line embedding step is a standard fact whose bibliographic source we did
   not verify.**  It is used only to move from maximal planar graphs to point sets, and for
   `n <= 24` we avoid it entirely by exhibiting explicit drawings.

## Related work

* Aichholzer, Hackl and Korman (2017) pack edge-disjoint *plane spanning trees and paths*;
  those objects have `n-1` edges, so a counting argument settles their extremal question,
  unlike ours.
* Bernhart and Kainen (1979) give the book thickness of `K_n`; a different problem with a
  numerically related answer.
* Füredi, Jiang and Kostochka (2020) and O'Neill and Spiro (2021) study Turán-type and
  saturation problems for convex geometric hypergraphs, with no non-crossing constraint on
  the family.

We did not find a published result for the convex extremal value, for the identity
`f(n) = a(n-2)`, or for the classification in Part 3.  We searched under the terminology
"edge-disjoint triangles", "non-crossing triangles", "triangle packing", "convex position",
"outerplanar graph", "maximal outerplanar graph", "weak dual" and "independence number of
a subcubic tree".  No priority claim is made beyond that, and no "first" is claimed.

## Citation

```
Edge-disjoint non-crossing triangle packings of planar point sets:
a convex barrier and a reduction to subcubic trees
```

See `paper/main.tex`.  If you use this, please cite the repository and say so.

## License

MIT.  See `LICENSE`.