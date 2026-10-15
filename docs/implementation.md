# Implementation map

Where each formula of these notes is evaluated in the package, and what the
package does that the notes do not need to say.

## One model

The package evaluates only the four-cochain stacking model: the generalized
differential $\mathfrak d_k$ ([differential.md](differential.md)) and the
product ([product_lower.md](product_lower.md), [product_upper_k1_k3.md](product_upper_k1_k3.md),
[product_upper_k4_k6.md](product_upper_k4_k6.md)). The page differentials are
its components:

| arrow | component | file |
| --- | --- | --- |
| $d_2$ (`Dbar`, `D`), $d_3$ from row $-2$ (`Dtilde`) | the single-layer components $P(A)$, $\tau_k(0,B)=D(B)$, pure-C $g_k$ | `gap/cochain_space.gi` (`nativePrimary`), cup-$i$ words of `gap/hap_space.gi` |
| $d_3$ from row $0$ (`Tau`) | $\tau_k(A,B)$ with $\delta B=P(A)$ | `gap/page_operations.gi` → worker task `a_step` (`QDa`, `fsharp`) |
| $d_4$ (`Psi`) | $\mathcal K_k(0,b,c)$ with $\delta c=D(b)$ | `gap/page_operations.gi` → tasks `qd`, `atom_curvature` |
| $d_5$ (`T`) | $\mathcal K_k(A,B,C)=\delta_s^R\,\Pi\,\Omega_k$ on a full defining system | `gap/page_operations.gi` → task `a_step` (`pot`) |

The defining cochains are solved on the resolution $R$ and transported as
$\Lambda r+Hz$ with the raw comparison homotopy of the group bar
([natural_evaluation.md](natural_evaluation.md), Section 6; the page tasks
of `python/light_rows.py`): a page class needs only the chain homotopy
identity, not the side conditions or the verified retraction that the
extension rows need, and the worker checks the pairing identities per
degree instead of assuming them. For `Psi` the D component on $(0,a,b)$ is
evaluated in its closed $A=0$ form, (S1)–(S7) of
[secondary_operations.md](secondary_operations.md), which is $\delta_s$ of a
rational potential of the same class as $\Omega_k(0,a,b)$ and costs cup
words instead of a $\theta$ prism on every comparison simplex. The
secondary and tertiary pages describe these cochains as formulas; they are
not separate code paths.

Shortcuts taken before any evaluation: an arrow with a trivial source or
target is zero (`gap/pages.gi`); an even integral input of `Tau` gives
$\rho\,\mathrm{Dtilde}(\rho(A/2))$; `Psi`$_0$ vanishes on a class killed by $D$;
`T`$_0$ is the native Pontryagin-square Bockstein (`KOAHSS_NativeT0Page`).

## The extension problem

`gap/extension.gi` reads the E6 layers and assembles the presentation
([extension_assembly.md](extension_assembly.md)); `gap/relations.gi` decides
per relation: split at an odd prime, zero by the target layer, a primary
row, or a light row (`gap/light_rows.gi`, `python/light_rows.py`), which is
the residue of transported defining data. The A-over-D rows use the
squaring telescope and the one-fiber primitives of
[power_relations.md](power_relations.md) in degrees five and six
(`python/unary_tables.py`, tables in `data/`, generator in `tools/`); the
diagonal source on a one-fiber simplex $\tau$ is read as the pair source on
$\Delta\tau$ through the persistent pair-source tables, so the shipped and
cached pair values serve the one-fiber primitive.

## Transport and the worker

`gap/bar_transport.gi` is the comparison with the group bar of a finite
group (used when it retracts), `gap/cell_transport.gi` the comparison with
the cell complex ($fg=1$ by construction) for every other resolution, with
the normalized homotopy and the primary comparison tensors
(`gap/primary_transport.gi`, `gap/binary_tensors.gi`, `gap/native_coherence.gi`).
Two Python workers per backend (`gap/worker_model.gi`,
`python/stacking_worker.py`): one over the group-bar comparison in page
mode for the page differentials, one over the normalized comparison for
every degree of the extension problem. Each keeps its transport chains,
the lazily built bar cochains and the formula memo tables between requests
(`python/acceleration.py`); the universal-value bundle is read one table at
a time, when a table is first used.

## Two exact shortcuts in the transport

- The worker caches transport answers by labelled simplex; GAP caches the
  encoded answer by kind and anchored labels (`gap/cell_transport.gi`),
  since every labelled simplex of one anchored class has the same encoded
  chain and the encoding (anchoring and labelling each term) costs more
  than the lookup. Vertex characters are cached by label.
- The mod-two reduction of a lift of an even vector is the structural zero
  (`python/acceleration.py`): every value of such a lift is an integer
  combination of even entries. Without it a page defining system with an
  even A evaluates its B and C homotopies on thousands of simplices whose
  sources vanish identically (space group 147 at degree four spent two
  thirds of its page phase there).

## The normalized homotopy and the cell comparison

The extension rows need the comparison homotopy with the side conditions
of a strong deformation retraction. With $q = 1 - gf$ and $u = q h q$ the
package ships $H = u\,\partial\,u$ (`KOAHSS_ExtensionNormalizedTransport`),
evaluated as $H = q\,h\,(1 - h\partial)$. Since $f$ and $g$ are chain
maps with $fg = 1$ on every verified basis element, $q$ is an idempotent
chain map with $qg = 0$, and $\partial h + h\partial = 1 - gf$ gives
$\partial u = q(q - h\partial)q = q - u\partial$; hence
$H = u(q - u\partial) = u\,(1 - u\partial)$. The raw homotopy $h$ of the
cell comparison vanishes on the image of $g$: it is equivariant, and on
the anchored chain $g(e_j)$ it is the cone at the base vertex of
$g(e_j) - gfg(e_j) - h\partial g(e_j) = -hg(\partial e_j)$, which is zero
by induction on the degree ($hg = 0$ on the degree-zero generators). So
$hq = h$, $u = qh$, and $H = qh - q(hq)h\partial = q\,h\,(1 - h\partial)$:
the raw homotopy of the faces, the raw homotopy of
$\sigma - h(\partial\sigma)$ and one projection at the end, where
$u\,\partial\,u$ projected four times, each projection expanding every
term through $f$ and the $g$ chains, and acted on $\partial(u\sigma)$
with $(n+2)\,|u\sigma|$ terms. The answer is the same chain; on space
group 30 with the pin twist the degree-3 requests of the A row took 4176 of
the case's 5055 seconds in the form $u\,\partial\,u$. In page mode (group
bar without a retraction) $u\,\partial\,u$ is kept, since $fg = 1$ is not
available there.

The cell comparison (`gap/cell_transport.gi`) works on integers: a vertex
$[g, \text{type}]$ is the number $\text{index}(g)\cdot 2^{20} + \text{type}$
with the first position of $g$ in `R!.elts` as index (the index of the
resolution's boundary and homotopy words, of the contraction memo and of
the binary tensor engine), R terms carry element indices, and products,
inverses and the converted contraction words `R!.homotopy(n,[j,index])`
are memoized in dictionaries keyed by one integer (a simplex or a pair
is encoded in base $2^{40}$; GAP's sorted-list dictionaries then search in
the kernel, whereas its hash tables run interpreted code per lookup), so
translation, anchoring, sorting and lookups never touch a group element
after its first product or inverse. The verifier, the normalized transport and the primary
comparison read `indexedElements` on the comparison record; the group-bar
comparison keeps its element-valued contract. Two choices of the modified
contraction depend on an order, the right inverse of the boundary
coordinates (`SmithNormalFormIntegerMatTransforms`) and the gcd
representation of a private boundary; both keep the previous order, by
basis and then by the group element itself, since another right inverse
is another valid contraction and moves the presentation of an extension
within its gauge (space group 30 with the pin twist changed one B
coordinate of its A row when the keys were ordered by index). `bench/dev_hcheck.g`
(working copy only) writes the normalized homotopy of every simplex of the
comparison chains of a case, so that two builds can be compared term by
term.

The transport bounds (`KOAHSS_EXTENSION_TRANSPORT_LIMITS`) are the
budgets of the previous package's campaigns, `maxSupport` 131072 and
`maxTerms` 20000000 (the former 8192 and 2000000 refused space groups 115
and 194 with the pin twist), and the retained normalized answers are
bounded by `answerCacheEntries` 4096, `answerCacheTerms` 8000000 and
`answerCacheMaxTerms` 2000000, the encoded answers by
`encodedCacheTerms` 20000000 (formerly 128 answers, 200000 terms and
`maxTerms`, which computed large answers twice).

## Characters and contractions on space groups

The sign character of a vertex, $(-1)^{s(\text{path})}$ along the
contracting homotopy from the base vertex, is a homomorphism $G \to \pm1$
(two paths with the same ends differ by a boundary, on which the cocycle
$s$ vanishes). On a space-group resolution one contraction per element
costs up to a third of a second and grows with the distance of the element
from the origin, so `gap/bar_transport.gi` assembles the character of an
affine matrix from the characters of the generators and of the lattice
basis translations: an element is a lattice translation times a
representative of its linear part, a word in the generators found by
factorization in the finite point group (Cryst's `TranslationBasis` gives
the lattice). The first three assembled values are checked against the
contraction. The contraction memo (`gap/resolution_memo.gi`) keeps up to
16384 contractions; a degree-four space-group case needs a few thousand
distinct ones and the old bound of 2048 recomputed many of them. Space
group 230 with the pin twist went from 172 s to 106 s of page phase on the
login node (the old package: 85 s; its B-layer source needs lifts on fewer
faces, which is a difference of the model rather than of the glue).

## Cost of the page worker

On a case whose pages need no page differential beyond the native
shortcuts, no worker is started. Otherwise the page worker costs about half
a second to start (the Python interpreter and the formula modules) plus
the group-bar comparison chains of the degrees it pairs; on a large
space-group resolution the page differentials themselves then take tenths
of a second. This is the only fixed cost that the previous package did not
have on trivial cases; the extension worker is the same in both.

## Caches

Universal values (`python/universal_values.py`, `universal_sources.py`),
the scaffolding of the universal models (certificates, Smith reductions;
`acceleration.persist_scaffolding`) and computed degrees
(`gap/result_cache.gi`) persist in `cache/`.
