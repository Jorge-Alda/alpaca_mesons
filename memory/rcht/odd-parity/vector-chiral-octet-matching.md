---
title: Enlarged vector action in the chiral-octet VVP matching limit
status: derived
created: 2026-10-13
updated: 2026-10-15
bibtex_keys: [Ruiz-Femenia:2003jdx, Roig:2013baa, Knecht:2001qf]
provenance: analytical-derivation
---

# Chiral-octet correlator matching

## Inputs and assumptions

Use the declared [nonideal quadratic action](../foundations/vector-nonideal-mass-action.md),
[mass extension](../foundations/vector-mass-extension.md),
[source extension](../foundations/vector-leptonic-width.md), and the
[odd action](../../extracts/Ruiz-Femenia:2003jdx/vector-odd-parity.md).
The source-level correlator inputs are `vvp.tex`, equations `eq:VVPres`
and `eq:short1`, archived under Ruiz-Femenia:2003jdx. This is tree level,
leading electromagnetic order, the chiral octet and the vector-only
single-trace odd basis; no pseudoscalar resonance or independent local odd
counterterm is added. All masses are in GeV.
The [convention dictionary](levi-civita-convention.md) fixes the metric
$(+---)$ and $\epsilon^{0123}=+1$, hence $\epsilon_{0123}=-1$.
Displayed couplings are $c_i=-c_i^R$, $d_i=-d_i^R$ and signed TFF scalars
are the negatives of the legacy scalars at corresponding couplings.
The linked original script/transcript retain their source convention;
the results here are their action-preserving adaptation, not a fresh
correlator calculation.

The enlarged-action mass spurion is set to zero **before** taking the hard
limit for this calculation. This is the chiral QCD OPE's domain, distinct
from the finite-spurion physical TFF calculation in
[the companion record](vector-physical-tff-matching.md).
Define $M_0$ as the chiral octet pole mass and retain $M_V$ as the independent
normalization scale in $c_i/M_V$. Neither is an average fitted physical pole.

## Executable derivation

Run [the script](../../../scripts/match_vector_short_distance.py):

```text
.venv/bin/python scripts/match_vector_short_distance.py
```

The [executed transcript](../../../data/vector_sector/short_distance_derivation.md)
contains SymPy-rendered inputs, expansion coefficients and solved relations.

1. In the chiral limit, all $e_m^V$, $\kappa_V$, $z_V$, $y_V$ terms
   proportional to the quark spurion vanish. The singlet mass term remains.
   In the ideal isoscalar coordinates its inverse operator is
   $K=(M_0^2-s)I+\Delta_Vnn^T$, with $n=(\sqrt2,1)^T/\sqrt3$.
   The electromagnetic isoscalar source is
   $F_V(1/(3\sqrt2),-1/3)^T$, orthogonal to $n$.
   SymPy inverts $K$ and obtains the common octet pole response; the singlet
   mass drops out. The same decoupling applies to traceless octet currents.
2. Build the full published correlator bracket, including the $c_3$ and
   $d_2$ direct **pseudoscalar-density source** terms. A pseudoscalar pole
   residue alone would miss these terms. Its scalar convention is
   $\Pi_R=-\langle\bar\psi\psi\rangle_0\mathcal B_R/F^2$.
   Its canonical scalar is $\Pi=-\Pi_R$, since the tensor defining the
   correlator reverses. The following expansion instructions describe the
   archived source-convention calculation.
3. Substitute $p^2\to p^2/z$, $q^2\to q^2/z$, $r^2\to r^2/z$,
   with $z=\lambda^{-2}$. Expand through $z^2$. Cancel the $z$ coefficient
   and equate the $z^2$ bracket coefficient to
   $F^2(p^2+q^2+r^2)/(2p^2q^2r^2)$, the cited chiral OPE.
4. Clear denominators, extract polynomial coefficients in the three
   invariants and solve the resulting equations with SymPy.

## Result

The independently solved conditions, adapted to canonical couplings, are

$$
4c_3+c_1=0,\qquad c_1-c_2+c_5=0,\qquad
c_5-c_6=-\frac{N_CM_V}{64\sqrt2\pi^2F_V},
$$

$$
d_1+8d_2=\frac{N_CM_0^2}{64\pi^2F_V^2}
-\frac{F^2}{4F_V^2},\qquad
d_3=\frac{N_CM_0^2}{64\pi^2F_V^2}
-\frac{F^2}{8F_V^2}.
$$

For $M_0=M_V$ these reproduce the five adapted literature relations. Even-sector
spurion extensions do not generate replacement octet OPE relations at this
order, because they vanish in its chiral limit. Their fitted physical residues
must not be substituted for the chiral $F_V$ without a source expansion.

## Incompatibility with the current TFF prescription

The script independently takes the massless common-pole limit of the actual
normalized TFF kernel. Applying the correlator conditions gives

$$
\lim_{Q^2\to\infty}F_P(-Q^2,0)
=-\frac{\sqrt2q(Q,Q)F}{2M_V^2}.
$$

This is nonzero for the pion: $-F/(3M_V^2)$ in its current generator
normalization. It violates single-virtual $1/Q^2$ falloff. The mismatch is
also discussed in the [Roig source extract](../../extracts/Roig:2013baa/tff-correlator-compatibility.md).

Solving instead for simultaneous single- and symmetric-double-virtual
vanishing in the same massless common-pole kernel requires

$$
c_1-c_2+c_5=0,\qquad
c_5-c_6=-\frac{N_CM_V}{64\sqrt2\pi^2F_V},\qquad
d_3=\frac{N_CM_V^2}{64\pi^2F_V^2}.
$$

The $-F^2/(8F_V^2)$ term is lost. SymPy then obtains
$\lim Q^2F_P(-Q^2,-Q^2)=0$, whereas the pion symmetric OPE requires
$-2F/3$ in the fixed external-state phase used here; the original
positive-anomaly scalar phase is documented in the
[pion target extract](../../extracts/Knecht:2001qf/pion-tff-short-distance.md).
Vanishing alone is therefore insufficient matching. This is a result for
the declared common-pole/chiral kernel, not a theorem excluding every
single-vector theory with additional interactions.

## Verification and limits

Substitution into both correlator expansion coefficients gives exactly zero
residual. The common-limit TFF checks pass symbolically. No measured coupling
value or Gaussian theory covariance is inferred. The singlet-density OPE,
finite-quark-mass perturbative contributions and additional correlator hard
limits have not been rematched. This record establishes the octet limit and
an explicit failure of simultaneous matching in the present action, not a
complete finite-spurion $U(3)$ VVP prescription.
