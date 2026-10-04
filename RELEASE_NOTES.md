## Research question

For a finite set `P` of points in general position in the plane, let `nu(P)` be the largest
number of triangles spanned by `P` such that **no two triangles share an edge** and **no
edge of one triangle crosses an edge of another**. How much does `nu(P)` depend on the shape
of `P`, and what are the extremal values?

The union of a packing is always a plane graph, so a counting argument bounds `nu(P)` by
`|E|/3`. In convex position that plane graph is forced to be *outerplanar* (`2n-3` edges);
with interior points it may be a full planar graph (`3n-6` edges). Convex position is
therefore expected to be the **hardest** configuration, and this release determines both
endpoints exactly.

## Main theorems (all proved in `paper/main.tex`)

| Statement | Result |
|---|---|
| **Theorem 4.1** | For `n` points in convex position, `f(n) = floor((2n-3)/3)`. |
| **Theorem 5.1** | `f(n) = a(n-2) = floor((2n-3)/3)`, where `a(N)` is the maximum independence number of a tree of maximum degree at most 3. |
| **Lemma 3.4** | Every subcubic tree on `N` vertices is the weak dual of a triangulation of the convex `(N+2)`-gon. |
| **Lemma 3.5** | Every `3`-cycle of a convex-polygon triangulation bounds a face, so the optimum inside `T` is `alpha(D(T))`. |
| **Lemma 5.2** | `a(N) = floor((2N+1)/3)` (a vertex-cover count, with a matching caterpillar construction). |
| **Proposition 6.2** | An Eulerian maximal planar graph packs `n-2` triangles, which partition its edges. |
| **Lemma 6.3** | Eulerian triangulations exist exactly for `n = 3` or `n >= 6` with `n != 7`. |
| **Theorem 6.5** | `nu(P) <= n-2`, with the equality cases classified. |

The headline reduction is Theorem 5.1: **in convex position the problem is not really
geometric.** Two elementary observations --- no vertex of a convex point set lies inside a
triangle of three others, and every subcubic tree is such a weak dual --- turn a geometric
extremal function into a two-line combinatorial one.

## The open part, stated honestly

**Conjecture 7.1 is not proved.** It claims that convex position is the worst configuration,
i.e. `nu(P) >= floor((2n-3)/3)` for every `n`-point set, with equality only in convex
position. It survives 423 tested configurations for `n <= 11` across eleven structured
families plus more for `n <= 12`, but a finite search is not a proof. The paper states
precisely why the obvious induction fails: peeling three points plus one triangle would
require `f(n) >= 1 + f(n-3)`, but `f(n) - f(n-3) = 2`, and two edge-disjoint triangles need
five vertices, so a three-point peel cannot supply two.

Two further caveats are marked in the paper:

* The implication `nu = n-2 => Eulerian triangulation` is verified exhaustively for
  `n <= 7` and is **not** proved in general.
* The straight-line embedding step that moves from maximal planar graphs to point sets is a
  standard fact whose bibliographic details we could not verify; it is flagged in the text
  and avoided entirely for `n <= 24` by exhibiting explicit drawings.

No peer review has taken place and nothing here has been submitted or published.

## Computational contribution

* `f(n)` computed by two independent routes --- a two-parameter dynamic program, and
  exhaustive enumeration of all `C_{n-2}` triangulations with an exact 0/1 program for each.
  They agree for `n <= 10`; the dynamic program matches the closed form for `n <= 400`.
* All `1 + 10 + 195 + 5712` labelled maximal planar graphs on `n = 3..7` enumerated and
  solved exactly. For `n = 6` **only** Eulerian graphs reach `n-2`; non-Eulerian ones reach
  only `3`.
* Explicit Eulerian triangulations for every `n <= 24`, each drawn with straight lines and
  its packing checked with exact geometric predicates.
* Both structural lemmas verified by exhaustion; `a(N)` cross-checked against 3000 random
  subcubic trees per `N`.
* 47 pytest tests aimed at catching wrong mathematics rather than at coverage.

## Reproducing

```bash
pip install -e .
python -m pytest tests -q            # 47 tests, ~6 s
python experiments/run_all.py        # every experiment, ~5 min
python experiments/run_all.py --quick   # skip the n=7 planar scan
cd paper && tectonic -X compile main.tex --keep-intermediates   # run twice
```

Everything is deterministic; the few random searches use fixed seeds. A clean-room check
(deleting `experiments/results/` and `figures/` and re-running) reproduces every number
quoted above.

## Novelty position

We did not find a published result for the convex extremal value, for the identity
`f(n) = a(n-2)`, or for the Eulerian classification, after searching OpenAlex, arXiv and
Crossref under the terminology "edge-disjoint triangles", "non-crossing triangles",
"triangle packing", "convex position", "outerplanar graph", "maximal outerplanar graph",
"weak dual" and "independence number of a subcubic tree". The nearest neighbours are plane
spanning-tree packing (Aichholzer--Hackl--Korman 2017), the book thickness of `K_n`
(Bernhart--Kainen 1979) and Turan-type problems for convex geometric hypergraphs
(Füredi--Jiang--Kostochka 2020); all are genuinely different problems. **No priority claim
is made beyond that, and no "first" is claimed.** Every reference was checked against
Crossref or OpenAlex: this corrected a wrong title and venue for Bernhart--Kainen, and
removed an unverifiable Fáry citation.

## Known limitations

1. Conjecture 7.1 is unproved.
2. One direction of the equality classification in Theorem 6.5 is machine-verified only, for
   `n <= 7`.
3. The convex value is a two-line counting bound; the content of that part is the
   construction and the tree reduction, not the value.
4. No variants: convex `k`-gons instead of triangles, weighted versions, and
   overlap-permitting versions are all untouched.

MIT licensed.