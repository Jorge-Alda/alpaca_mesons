---
title: Flavor selection rules for axial and pseudoscalar ALP to four-pseudoscalar amplitudes
status: derived
created: 2026-10-23
updated: 2026-10-23
bibtex_keys: [Alda:2026nrd, Ecker:1988te, Ruiz-Femenia:2003jdx]
provenance: analytical-derivation
---

# Scope and inputs

This record gives flavor selection rules for the two axial/pseudoscalar-source
amplitudes in the [four-pseudoscalar calculation](alp-four-pseudoscalar-axial-pseudoscalar.md).
The source matrices are in the no-mixing basis, and charged external mesons
use outgoing conjugate generators. The physical neutral pseudoscalar
generators and nonideal $\omega$--$\phi$ basis are those of
Alda:2026nrd and the preliminary-plus-$G_V$ fit [@Alda:2026nrd;
@Ecker:1988te; @Ruiz-Femenia:2003jdx]. For the structural audit, their
diagonal entries are left independent and symbolic. The resulting counts are
generic flavor support before substituting the fitted neutral generators or
vector mixing angle; they do not define a different mixing prescription.
This is a tree-level flavor audit of the retained two-vector exchange
families. It does not include electromagnetic interactions, loops,
independent $VPPP$ operators, or other interactions outside that fit
truncation.

The reproducible symbolic audit and its tests are in
[`scripts/derive_alp_four_pseudoscalar.py`](../../../scripts/derive_alp_four_pseudoscalar.py).
Run

```text
.venv/bin/python scripts/derive_alp_four_pseudoscalar.py --selection-rules
.venv/bin/python scripts/derive_alp_four_pseudoscalar.py --test
```

# Trace conditions

Each $VPP$ vertex must have a nonzero flavor commutator trace,
$$
\operatorname{tr}\!\left(T_R[\lambda_i^\dagger,\lambda_j^\dagger]\right)\ne0.
$$
Equal mesons and pairs of diagonal neutral mesons therefore cannot form a
$VPP$ vertex. The commutator fixes the exchanged vector direction: for
example, $\pi^+\pi^-$ selects the neutral $\rho^0$ direction, whereas
$K^+K^-$ can select neutral $\rho^0,\omega,\phi$ directions. A $K\pi$ pair
selects a strange vector direction fixed by its charge and strangeness.

At the odd $aVV$ source, at least one of the derivative and pseudoscalar-mass
traces must be nonzero:
$$
\operatorname{tr}\!\left(\{T_R,T_S\}c_q^A\right),\qquad
\operatorname{tr}\!\left(\{T_R,T_S\}\{c_q^P,\chi\}\right).
$$
At the odd strong vertex in the axial chain, the corresponding traces are
$$
-\frac{1}{2f_\pi}\operatorname{tr}\!\left(\{T_R,T_S\}\lambda_j^\dagger\right),
\qquad
\frac{1}{2f_\pi}\operatorname{tr}\!\left(\{T_R,T_S\}
\{\lambda_j^\dagger,\chi\}\right).
$$
The axial $aVP$ vertex additionally requires
$$
\operatorname{tr}\!\left(T_R^\dagger[c_q^A,\lambda_i^\dagger]\right)\ne0.
$$
These trace tests determine allowed local vector routings. A zero of one
trace alone is not a selection rule if the other trace survives. Fitted
$d_i$ combinations, kinematics, or interference among external assignments
can still cancel a routing in the complete amplitude.

# Charge, strangeness, and admitted final states

The ALP is electrically neutral, so the four-pseudoscalar state must have
total charge zero. A flavor-diagonal source has total strangeness zero. With
the stated outgoing-generator convention, a single off-diagonal source
matrix unit selects
$$
c_q^{A,P}\propto E_{23}: S_{\rm final}=+1,
\qquad
c_q^{A,P}\propto E_{32}: S_{\rm final}=-1.
$$
A Hermitian $d$--$s$ source contains both orientations and therefore both
conjugate strangeness sectors. A diagonal source with $c_q^A\propto\mathbf1$
has no axial $aVP$ vertex because its commutator vanishes; its odd $aVV$
family need not vanish.

The symbolic trace enumeration finds 30 flavor-admitted four-pseudoscalar
multisets for a generic diagonal source, 29 for $E_{23}$, and 29 for $E_{32}$.
All have total charge zero and the required strangeness. No all-neutral final
state is admitted, because every routing needs two nonzero $VPP$
commutators.

The three diagonal-source channels
$\pi^+\pi^-\eta\eta$, $\pi^+\pi^-\eta\eta'$, and
$\pi^+\pi^-\eta'\eta'$ depend on the small unequal $u,d$ entries in the
physical $\eta,\eta'$ generators. For example,
$$
\operatorname{tr}\!\left(T_{\rho^+}
[\lambda_{\pi^+}^\dagger,\lambda_\eta]\right)
=\sqrt2(\eta_u-\eta_d),
$$
which is zero in the exact-isospin approximation and nonzero for the central
Alda:2026nrd generators. The audit keeps independent symbolic diagonal
entries for $\pi^0,\eta,\eta'$ to preserve this support; it does not silently
impose exact isospin.

# Exchanged vector directions

For diagonal sources, the odd $aVV$ vertex permits conjugate charged-vector
pairs $\rho^+\rho^-$, $K^{*+}K^{*-}$, and $K^{*0}\bar K^{*0}$, as well as
neutral-neutral pairs among $\rho^0,\omega,\phi$ when their source traces are
nonzero. For $E_{23}$, the allowed odd-source pairs are
$(K^{*+},\rho^-)$ and $(K^{*0},\rho^0)$, $(K^{*0},\omega)$,
$(K^{*0},\phi)$; $E_{32}$ gives the charge-conjugate set. The fitted
nonideal neutral-vector mixing determines the numerical $\omega$ and $\phi$
residues.

In the axial chain, all members of the vector nonet occur somewhere in the
generic support for each source sector. This does not mean that every
ordered vector pair is allowed: the axial commutator trace, the strong odd
trace, and the $VPP$ trace must all be nonzero for the particular routing.
The `--selection-rules` output gives the full multiset-to-vector assignments.

# Evidence boundary

The counts and routing list are generic symbolic trace support before
substituting the fitted neutral generators and before summing the full
momentum-dependent Bose permutations. They identify flavor-admitted states,
not guaranteed nonzero physical amplitudes. Fitted-coupling cancellations,
special kinematics, thresholds, and destructive interference can remove an
admitted channel. Identical-particle factorials enter widths, not these
amplitude-level selection rules.
