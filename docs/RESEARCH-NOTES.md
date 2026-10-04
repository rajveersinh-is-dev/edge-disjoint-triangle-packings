# Research notes: how this project was chosen, and what was abandoned

This file records the selection process, the novelty audit, and --- importantly --- the
candidates that were **rejected**.  It is kept in the repository so that a reader can judge
how much of the search space was actually explored, and so that the negative results are not
silently discarded.

## 1. Candidates considered and why they were rejected

Roughly twenty directions were screened across discrete geometry, extremal graph theory,
combinatorics on convex sets, and computational geometry.  The ones that were seriously
pursued and then dropped:

| Candidate | Why dropped |
|---|---|
| **Edge-disjoint triangulations of a convex polygon** (`tau(n)`, `kappa(n)`) | **Rediscovery.** The answer `floor(n/2)` follows from a two-line count (each triangulation uses `n-3` of the `n(n-3)/2` diagonals) matched by the classical double-fan. The same value is the book thickness of `K_n` (Bernhart–Kainen). No new mathematics. |
| Max # edge-disjoint non-crossing triangles, convex position | Initially the main candidate, and the DP was built and verified. It turned out to be exactly the edge-count bound `floor((2n-3)/3)`, i.e. the upper bound is free. **Kept, but demoted**: the content became the reduction to subcubic trees (`f(n) = a(n-2)`), not the value. |
| Same problem with convex *quadrilaterals* / `k`-gons | Edge count `floor((2n-3)/k)` looked like it might fail to be tight, but a short analysis suggested the counting bound is again the answer; the interesting variant was not isolated. Left as future work. |
| Vertex-disjoint non-crossing `k`-gons | Degenerates to `floor(n/k)`; no content. |
| "Convex geometric hypergraph" Turán problems | Heavily studied (Füredi–Jiang–Kostochka, O'Neill–Spiro, Frankl–Kupavskii); the specific questions we could formulate were either known or not precisely stateable. |
| Triangle packings in graphs with a crossing-number budget | The extremal question is open but far out of reach of the tools here. |
| Max # edge-disjoint plane spanning trees | Known (`floor(n/2)`), counting bound plus zig-zag. |
| **Covering all diagonals of a convex `n`-gon by triangulations** (an earlier direction in this workspace, files since deleted) | Genuinely rich (a covering number, plus an "overlap graph" on candidate faces whose chromatic number decides which families are realisable). Abandoned because the interesting structural hypothesis was never isolated and the tools needed were the same ones used here. |

## 2. Novelty audit

Searched (OpenAlex, arXiv API, Crossref) using: *"maximum number of edge-disjoint
non-crossing triangles"*, *"triangle packing geometric graph"*, *"convex geometric
hypergraph matching"*, *"plane graph triangle packing"*, *"partition edges maximal outerplanar
graph into triangles"*, *"independence number subcubic tree"*, *"weak dual triangulation"*.

Findings that shaped the paper:

* **Nothing** matched the convex-position extremal value or the identity `f(n) = a(n-2)`.
* The closest neighbours are plane spanning-tree/path packing and convex geometric
  hypergraph Turán problems; both are cited in the paper and both are genuinely different
  (different objects, different constraint).
* **A citation in the previous version of this repository was wrong.**  It cited
  "Bernhart and Kainen, *Book thickness of graphs*, Discrete Mathematics 27 (1979),
  295--301".  Crossref shows the real paper is *"The book thickness of a graph"*, *Journal of
  Combinatorial Theory, Series B* **27** (1979), 320--331.  Fixed.
* **Fáry (1949) could not be verified.**  Crossref and Numdam both resolve
  `BSMF_1949__77__128_0` to a *different* Fáry article.  Rather than invent a
  bibliographic entry, the paper states the straight-line embedding fact as standard, flags
  it, and avoids needing it for `n <= 24` by exhibiting explicit drawings.
* The Frankl–Kupavskii attribution for "Extremal problems for convex geometric hypergraphs
  and ordered hypergraphs" was wrong; Crossref gives **Füredi, Jiang and Kostochka**, *Can. J.
  Math.* **73** (2021), 1648--1666.  Fixed.

## 3. Method

```
  hypothesis  ->  small-instance enumeration  ->  pattern discovery  ->  proof attempt
      ^                                                                        |
      |                                                                        v
  new counterexample search  <-  computational verification  <-  proof or refutation
```

Concretely, the order was:

1. **Guess.** Non-crossing + edge-disjoint triangle packings on convex point sets.
2. **Enumerate.** Dynamic program for `f(n)`, plus exhaustive enumeration of all
   triangulations with an exact ILP for each.  Found `f(n) = floor((2n-3)/3)`, i.e. the
   edge-count bound --- so the value itself is not the contribution.
3. **Ask why.** Because every `3`-cycle of a convex-polygon triangulation is a face, the
   problem is an independence number in a weak dual tree.  Because every subcubic tree is
   such a weak dual, `f(n) = a(n-2)`, and `a(N) = floor((2N+1)/3)` has a two-line proof.
   *This* is the contribution.
4. **Falsify.** Eleven structured families of point sets; 423 configurations for `n <= 11`
   and more for `n <= 12`; no configuration dipped below the convex value.  The search for
   counterexamples to Conjecture 7.1 is a committed, re-runnable script, not an anecdote.
5. **Go the other way.** The union of a packing is planar, so `nu(P) <= n-2`; equality is a
   bipartite dual condition, i.e. an Eulerian triangulation; the existence question has the
   clean answer `n = 3` or `n >= 6, n != 7`.

## 4. Bugs the falsification effort caught

Recorded because they are the strongest evidence that the experiments were doing their job:

* `segments_cross` reported a *shared endpoint* as a crossing, which made every triangle
  conflict with every other and forced `nu(n) = 1` for all `n`.  Found by comparing the
  solver against the dynamic program on a convex pentagon.
* The triangulation enumerator emitted only fans while still producing the Catalan *count*,
  so the "exhaustive" check silently examined the wrong objects.  Found because the
  star triangulation `(1,3),(3,5),(1,5)` was missing from the enumeration.
* `triangulation_edges` stored the boundary edge `(n-1, 0)` unnormalised, so `3`-cycles were
  missed.
* The constructive packer's split arithmetic was off by one (`A+B = m+2`, not `m+1`), its
  right sub-polygon was shifted by `a` instead of `a-1`, and the root triangle used the
  wrong vertex.  Caught by validating every constructed packing geometrically.
* `alpha_max_subcubic` first used a Pareto frontier where a plain max-plus convolution
  suffices, giving `a(N) = N-1`; and an earlier version required exactly `s` non-empty
  child subtrees, giving `a(4) = 2` instead of `3`.  Caught by comparing against the star
  `K_{1,3}` and by randomised search.
* The Eulerian "connected sum" identified the wrong node set (12 vertices instead of 9 for
  `n = 9`), and used vertex triples that were not faces.  Caught by checking planarity and
  vertex count of every construction.

## 5. Dead ends worth recording

* **Proving Conjecture 7.1 by peeling.**  Peeling three hull vertices plus one triangle
  requires `f(n) >= 1 + f(n-3)`, but `f(n) - f(n-3) = 2`; and two edge-disjoint triangles
  need five vertices, so a three-point peel cannot supply two triangles.  The arithmetic is
  the obstruction.
* **Empty hull ears.**  The natural lemma --- "every point set has an *empty* hull ear, and an
  empty hull ear can be peeled" --- is false: a quadrilateral with one interior point in each
  of its four ear triangles has no empty hull ear at all.  (The convex case is fine, because
  there are no interior points.)
* **Segment split.**  Splitting a set by a line through two points makes the two halves share
  two vertices, so `|P_1| + |P_2| = n + 2` instead of `n + 1`, and the `$+1$` branch of the
  recurrence (which is worth exactly one triangle) cannot be used.  This is the gap noted in
  the paper.

## 6. Honest assessment of novelty

*Part 2* (the identity `f(n) = a(n-2)` and `a(N) = floor((2N+1)/3)`) is the part we believe
is genuinely new: it removes geometry from a geometric extremal problem in two elementary
steps.  *Part 3* (the Eulerian-triangulation characterisation of when `n-2` is attained) is
very likely a standard exercise in disguise, and we present it as a clean, completely proved
result rather than claiming novelty.  *Part 1* (the convex value) is a counting bound.  The
unproved conjecture is stated as such.