---
title: Single/double-trace flavor enumeration for the bounded odd vertex sector
status: derived
created: 2026-10-13
updated: 2026-10-13
bibtex_keys: [Alda:2026nrd, Ruiz-Femenia:2003jdx, Roig:2013baa]
provenance: analytical-derivation
---

# Abstract flavor tensors

## Inputs and scope

Stage A of the [selected plan](odd-breaking-low-energy-calibration-plan.md).
Use the [Ruiz operator basis](../../extracts/Ruiz-Femenia:2003jdx/vector-odd-parity.md),
[VVP basis](../../extracts/Ruiz-Femenia:2003jdx/vvp-three-pseudoscalar.md),
and [Roig trace/basis comparison](../../extracts/Roig:2013baa/anomalous-rcht-bases-and-double-trace.md).
Flavor matrices are abstract $3\times3$ Hermitian covariant blocks, with no
commutativity, tracelessness or equality assumed before the indicated projection.
The Lorentz components of two V fields remain separately labelled until Bose
symmetrization. This is the bounded one-odd-block vertex sector, not an
arbitrary nonlinear $U(3)$ generating functional or an EOM-minimal action.

The symbols in this enumeration are field slots, not the coupling combinations
called $A,B$ in the matching record:

- $A$: one tensor vector component, C-odd under transposition.
- $B$: $f_+$ for VJP, or a second tensor vector for VVP, also C-odd.
- $U$: one $u_\mu$ component or its indicated covariant derivative, C-even.
- $X$: the full dynamical $\chi_+$ block, C-even; it is not $X_{TC}$.
- $Y=i\chi_-$: Hermitian C-even odd source block, with exactly one spurion.
- In the separate $f_-$ sector, $B=f_-$ is C-even; $A$ remains C-odd.

The intrinsic-parity assignments are those of the cited covariant blocks.
All vertices have one Levi-Civita tensor and real action coefficients.

## Executable method

Run [the enumeration script](../../../scripts/enumerate_odd_breaking_vertices.py).
The [transcript](../../../data/vector_sector/odd_vertex_enumeration.md)
prints SymPy-rendered traces; the
[JSON catalogue](../../../data/vector_sector/odd_vertex_enumeration.json)
stores every trace word, projection relation and input hash.

Enumerate all permutations and partitions into one or two nonempty traces.
Canonicalize cyclic words and commuting scalar-trace order. SymPy `Trace`
does not canonicalize those finite label orbits, so the script implements
only that combinatorial enumeration; matrix products and traces use
`MatrixSymbol`, `MatMul`, `Trace` and ordinary SymPy expressions.

For the C-even three/four-block sectors, Hermitian/C projection gives
$(t+t^{\rm rev})/2$. For $Vf_-X$, C is odd on the product, giving
$i(t-t^{\rm rev})/2$. Equal reverse words vanish in the latter case.

## Three-block result

For $A,B,Z$ with $Z=U$ or $Y$, define

$$
F_0=\frac12\operatorname{tr}((AB+BA)Z),\quad
F_1=\operatorname{tr}(A)\operatorname{tr}(BZ),\quad
F_2=\operatorname{tr}(B)\operatorname{tr}(AZ),\quad
F_3=\operatorname{tr}(AB)\operatorname{tr}(Z).
$$

There is one single-trace and three double-trace tensors. Their VVP exchange
subspaces have dimensions three symmetric and one antisymmetric: the latter
is represented by $F_1-F_2$. It survives only when multiplied by an
antisymmetric Lorentz/momentum tensor so the full VVP vertex is Bose-even.
For $Vf_-X$ there is only one tensor in this trace ceiling,
$i\operatorname{tr}(A[B,X])/2$; double-trace charge projections vanish.
Normalization is declared: the published anticommutator trace is $2F_0$.

## Four-block one-spurion result

For $A,B,U,X$, the catalogue contains three C-even single traces and seven
C-even double traces. The double partitions are the four $3+1$ partitions
and three $2+2$ partitions. Its exchange decomposition is:

| Topology | Symmetric under $A\leftrightarrow B$ | Antisymmetric |
|---|---:|---:|
| Single trace | 2 | 1 |
| Double trace | 5 | 2 |
| Total | 7 | 3 |

Every explicit expression and its ordering is in the JSON; no fitted
coefficient is inferred. SymPy row reduction verifies both exchange subspaces.

The script constructs the polarized degree-four $3\times3$ Cayley-Hamilton
trace identity from permutation cycles. It contains nonzero three/four-trace
terms. It does not furnish a redundancy wholly within these one/two-trace
degree-four tensors when all block traces are independent. In particular,
dropping its higher-trace terms would create a spurious relation. The
catalogue does not apply an octet/traceless identity to the $U(3)$ singlet.

## Singlet and physical projections

Decompose $X=\widehat X+\operatorname{tr}(X)I/3$. The script sets $X=I$
in every four-block tensor and expresses the result exactly in the four
zero-spurion tensors. Thus for a constant singlet background the corresponding
corrections are rescalings of lower-block coefficients, not new separately
identifiable low-energy parameters. A derivative of the singlet spurion is
not removed by that argument. The trace topology before this projection
remains recorded separately.

Exact polynomial reductions with the repository's concrete diagonal flavor
matrices give:

| Flavor sector | Abstract count | Diagonal $A,T_P$, $B=Q$ | Two arbitrary diagonal vector slots |
|---|---:|---:|---:|
| $ABU$ | 4 | 3 | 4 |
| $ABY$ | 4 | 3 | 4 |
| $ABUX$ | 10 | 7 | 8 |

Here $Q=\operatorname{diag}(2/3,-1/3,-1/3)$, $X=2\chi$, with independent
odd quark-spurion entries; the $Y$ projection keeps the external density
independent. These are flavor-polynomial ranks, not actual-data ranks.
The concrete matrices are physical projection inputs, not replacements for
the abstract flavor reduction. `Tr(f_+)=0` for the electromagnetic source
eliminates singleton-$B$ traces. Neutral commutativity collapses three
single-trace placements to one, whereas charged matrices need not do so.

## Existing versus new double traces

The script inverts the cited Roig dictionary, obtaining

$$
d_1=-8n_f\kappa^{VV}_1,\qquad
d_2=\kappa^{VV}_2+n_f\kappa^{VV}_1,
$$

$$
c_2=c_6-4n_fM_V\kappa^V_{13},\qquad
c_3=M_V(\kappa^V_{14}+n_f\kappa^V_{13}).
$$

Thus a coefficient named “double trace” in a different EOM-reduced basis
can already be correlated with the original $d_1$ or $c_2,c_6$; simply
adding that literature coefficient as a second independent copy would
double count. These equations are a literature-basis audit, not imposed
off-shell source relations. Independent external density/axial terms are
retained here; applying the EOM dictionary to the dynamical Alda singlet
requires its source/anomaly contact terms and is deferred to normalized
source expansion in stage B.
