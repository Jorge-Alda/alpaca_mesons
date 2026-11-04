---
title: One-loop effects and renormalized constants in the quadratic U(3) Lagrangian
status: mixed
created: 2026-11-05
updated: 2026-11-05
bibtex_keys: [Herrera-Siklody:1996tqr, Kaiser:1998ds, Gasser:1983yg, Alda:2026nrd]
provenance: literature
---

# One-loop effects and renormalized constants in the quadratic U(3) Lagrangian

## Scope and inputs

The local NLO quadratic terms and conventions are those of [the NLO quadratic record](quadratic-lagrangian.md), built from [the retained NLO basis](lagrangian.md). The fields are the pseudoscalar nonet plus the ALP, and the source convention is the reduced $V_4=0$ basis. The loop statements below refer to the $U(3)$ nonet theory of [@Herrera-Siklody:1996tqr], with its stated dimensional-subtraction convention, and to the usual $SU(3)$ convention for the standard $L_i$ [@Gasser:1983yg]. They are not a new calculation of ALP-dependent loop graphs.

## What the loops add

At NLO, one-loop graphs formed from the leading chiral Lagrangian contribute alongside tree graphs from the local NLO Lagrangian. They have two distinct parts:

1. Ultraviolet-divergent local polynomials. These are absorbed into bare NLO couplings and leave finite, scale-dependent $L_i^r(\mu)$.
2. Finite nonlocal functions of external momenta and meson masses. These contain the physical threshold and unitarity-cut dependence and cannot be represented by a finite shift of the local constants. They must be added explicitly to two-point functions and amplitudes.

The one-loop effective action also gives meson self-energies, wavefunction corrections, and shifts of pole masses and mixing eigenvectors. For the $U(3)$ theory the anomaly term changes the divergence pattern: in general, counterterms can be induced at more than one chiral order. In the anomaly-free degenerate-mass example discussed in [@Herrera-Siklody:1996tqr], the octet mass and kinetic divergences have the same factor and can be removed by an octet field rescaling, while the singlet mass has a divergence without a matching singlet kinetic divergence. With the anomaly included, both octet and singlet masses receive divergent contributions. These points mean that the local quadratic Hessian in the [quadratic-lagrangian](quadratic-lagrangian.md) is not the complete NLO inverse propagator: one must add the loop self-energy matrix, which is generally momentum dependent and has meson thresholds. For ALP mixing, loops also correct the mixed two-point functions and hence the pole residues/LSZ projections; these terms are not obtained by merely replacing each $L_i$ by $L_i^r$.

The finite nonlocal terms are fixed once the leading theory, masses, sources, and subtraction prescription are specified. Their scale dependence cancels that of the renormalized local constants in physical observables. Their explicit expressions are process and kinematics dependent, so no single universal loop correction can be appended to the quadratic record without choosing a two-point function and evaluating its diagrams.

## Renormalization of the standard $L_i$

In the standard three-flavor $SU(3)$ chiral subtraction scheme, write

$$
L_i=L_i^r(\mu)+\Gamma_i\Lambda(\mu),\qquad
\Lambda(\mu)=\frac{\mu^{d-4}}{(4\pi)^2}\left[\frac{1}{d-4}-\frac12(\log4\pi+\Gamma'(1)+1)\right].
$$

Equivalently, at one loop,

$$
\mu\frac{d L_i^r}{d\mu}=-\frac{\Gamma_i}{16\pi^2},\qquad
L_i^r(\mu_2)=L_i^r(\mu_1)-\frac{\Gamma_i}{16\pi^2}\log\frac{\mu_2}{\mu_1}.
$$

For the five conventional constants occurring in the retained quadratic operators, the coefficients are

$$
\Gamma_4=\frac18,\qquad \Gamma_5=\frac38,\qquad
\Gamma_6=\frac{11}{144},\qquad \Gamma_7=0,\qquad
\Gamma_8=\frac{5}{48}.
$$

Thus explicitly,

$$
\begin{aligned}
L_4&=L_4^r(\mu)+\frac18\Lambda(\mu),&
L_5&=L_5^r(\mu)+\frac38\Lambda(\mu),\\
L_6&=L_6^r(\mu)+\frac{11}{144}\Lambda(\mu),&
L_7&=L_7^r(\mu),\\
L_8&=L_8^r(\mu)+\frac{5}{48}\Lambda(\mu).
\end{aligned}
$$

These coefficients are for the conventional $SU(3)$ renormalization scheme. They are not automatically the complete $U(3)$ counterterm prescription when the singlet and anomaly are dynamical. In the $U(3)$ treatment, the coefficient functions of the invariant singlet variable and the additional singlet-source operators also renormalize; see the full operator/counterterm construction in [@Herrera-Siklody:1996tqr]. In particular, the present records do not establish subtraction coefficients or running for $\Lambda_1,\Lambda_2,L_{18},L_{25}$ in the reduced large-$N_c$ basis. Those must be matched from the chosen $U(3)$ loop scheme, rather than inferred from the five $SU(3)$ values above.

## Constants entering the retained quadratic Lagrangian

The quadratic expansion already derived in the [quadratic-lagrangian](quadratic-lagrangian.md) contains the following local constants:

| Sector | Constants in the local quadratic terms | Role |
|---|---|---|
| Two-meson kinetic | $L_4^r,L_5^r,\Lambda_1^r,L_{18}^r$ | Octet and singlet-trace kinetic normalization and mixing |
| Two-meson mass | $L_6^r,L_7^r,L_8^r,\Lambda_2^r,L_{25}^r$ | Quark-mass and singlet/topological mass entries |
| Meson--ALP kinetic | $L_4^r,L_5^r,\Lambda_1^r,L_{18}^r$ | Derivative mixing entries |
| Meson--ALP mass | $L_6^r,L_7^r,L_8^r,\Lambda_2^r,L_{25}^r$ | Mass mixing entries |

Here $L_4^r,L_5^r,L_6^r,L_7^r,L_8^r$ mean the standard constants with the subtraction just displayed. Superscript $r$ on $\Lambda_{1,2}$ and $L_{18},L_{25}$ is a reminder that these are renormalized parameters in a specified $U(3)$ scheme; their explicit counterterm formulas are not supplied here. The ALP--meson terms are proportional to the ALP sources already shown in the[quadratic-lagrangian](quadratic-lagrangian.md). There is no claim that every listed constant has nonzero contribution for every source choice or flavor channel.

## Relation to the current quadratic record

Replacing the bare symbols by renormalized symbols in the existing local expressions is necessary but not sufficient for an NLO two-point function. The complete renormalized inverse propagator has the form

$$
\mathcal Q_{IJ}(s)=sK_{IJ}^{\rm LO}-M_{IJ,\rm LO}^2
+s\,\delta K_{IJ}^{\rm local}(\mu)-\delta M_{IJ,\rm local}^2(\mu)
+\Sigma_{IJ}^{\rm loop}(s;\mu),
$$

where $I,J$ range over the mesons and, when retained, the ALP. The sign convention for $\Sigma$ is defined by this displayed inverse-propagator equation. The renormalized local terms and loop self-energy are separately scale dependent; their sum is scale independent to the calculated order. The loop matrix includes both diagonal meson self-energies and meson--ALP mixed two-point functions if ALP-dependent vertices are included. Local polynomial pieces may be shifted between $\Sigma$ and the local constants by a finite scheme change, while the nonanalytic thresholds remain.

This note records the scope and renormalization structure, not explicit $\Sigma_{IJ}^{\rm loop}$ formulas. Such formulas require a specified ALP-source specialization and a diagrammatic calculation in the corresponding $U(3)$ scheme.
