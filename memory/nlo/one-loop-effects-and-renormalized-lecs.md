---
title: One-loop effects and renormalized constants in the quadratic U(3) Lagrangian
status: mixed
created: 2026-11-05
updated: 2026-11-05
bibtex_keys: [Herrera-Siklody:1996tqr, Kaiser:1998ds, Gasser:1983yg, Gasser:1984gg, Alda:2026nrd]
provenance: literature
---

# One-loop effects and renormalized constants in the quadratic U(3) Lagrangian

## Scope and inputs

The local NLO quadratic terms and conventions are those of [the NLO quadratic record](quadratic-lagrangian.md), built from [the retained NLO basis](lagrangian.md). The fields are the pseudoscalar nonet plus the ALP, and the source convention is the reduced $V_4=0$ basis. The loop statements below refer to the $U(3)$ nonet theory of [@Herrera-Siklody:1996tqr], with its stated dimensional-subtraction convention, and to the usual $SU(3)$ convention for the standard $L_i$, with the subtraction $\lambda(\mu)$ of [@Gasser:1983yg] and the $SU(3)$ coefficients $\Gamma_i$ of [@Gasser:1984gg]. They are not a new calculation of ALP-dependent loop graphs.

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

These coefficients are for the conventional $SU(3)$ scheme. They are not the $U(3)$ prescription when the singlet and anomaly are dynamical. In the $U(3)$ scheme used by [@Kaiser:1998ds],

$$
\Gamma_4=\frac18,\quad \Gamma_5=\frac38,\quad \Gamma_6=\frac1{16},\quad \Gamma_7=0,\quad \Gamma_8=\frac3{16}.
$$

Thus $L_6=L_6^r+\lambda/16$ and $L_8=L_8^r+3\lambda/16$ in this scheme, replacing the $SU(3)$ values above. Here

$$
L_i=L_i^r(\mu_\chi)+\Gamma_i\lambda(\mu_\chi),\qquad
\lambda(\mu_\chi)=\frac{\mu_\chi^{d-4}}{16\pi^2}\left[\frac1{d-4}-\frac12(\log4\pi+\Gamma'(1)+1)\right].
$$

Kaiser:1998ds states that $F,\tau,\Lambda_1,\Lambda_2$ do not acquire chiral-loop renormalization at this order. Hence $\Lambda_{1,2}$ have no pole in $\lambda$ at fixed QCD scale. They do depend on the QCD renormalization scale through the singlet axial-current anomalous dimension: $1+\Lambda_1$ transforms as $Z_A^2$ and $1+\Lambda_2$ as $Z_A$. The combinations $2L_5+3L_{18}$ and $2L_8-3L_{25}$ transform as $Z_A$.

The $U(3)$ chiral counterterms for the additional NLO constants are

$$
L_{18}=L_{18}^r(\mu_\chi)-\frac14\lambda(\mu_\chi),\qquad L_{25}=L_{25}^r(\mu_\chi),
$$

so $\Gamma_{18}=-1/4$ and $\Gamma_{25}=0$, consistent with the $n_l=3$ functions in [@Herrera-Siklody:1996tqr]. The vanishing chiral pole for $L_{25}$ does not remove its QCD-scale transformation through $2L_8-3L_{25}$.

For $L_{47}$, [@Herrera-Siklody:1996tqr] includes $O_{47}=D_\mu\hat\theta\langle C^\mu\rangle\langle M\rangle$ and gives a generally nonzero, singlet-function-dependent divergence. With that paper's definitions $\omega_i=W_i/W_1$ and $\omega_i^{(n)}=W_i^{(n)}/W_1$ (the superscript differentiates $W_i$, not $\omega_i$),

$$
L_{47}(X)=L_{47}^r(X,\mu_\chi)+B_{47}(X)\lambda(\mu_\chi),
$$

where

$$
\begin{aligned}
B_{47}={}&-\frac{n_l}{2}\omega_1^{(2)}\omega_2+\frac{i}{2}\omega_1^{(2)}\omega_3^{(1)}+\frac{n_l}{2}(\omega_1^{(1)})^2\omega_2\\
&+\frac14(n_l-n_l^3)\omega_1^{(1)}\omega_1^{(2)}\omega_2^{(1)}+i\left(-\frac12+\frac{n_l^2}{4}\right)\omega_1^{(1)}\omega_1^{(2)}\omega_3\\
&-\frac{i}{2}(\omega_1^{(1)})^2\omega_3^{(1)}-\frac14(n_l-n_l^3)(\omega_1^{(1)})^3\omega_2^{(1)}+\frac{i}{2}(1-n_l^2)(\omega_1^{(1)})^3\omega_3.
\end{aligned}
$$

The same paper then imposes its large-$N_c$ truncation to evaluate the counterterms. It expands

$$
W_1(X)=\frac{f^2}{4}(1+v_{12}X^2+\cdots),\quad
W_2(X)=\frac{f^2}{4}(1+v_{22}X^2+\cdots),\quad
W_3(X)=-i\frac{f^2}{4}(v_{31}X+\cdots),
$$

so at $X=0$, $\omega_1^{(1)}=0$, $\omega_1^{(2)}=2v_{12}$, $\omega_2=1$, and $\omega_3^{(1)}=-iv_{31}$. The general expression therefore reduces to

$$
B_{47}(0)=v_{12}(v_{31}-n_l)=v_{12}(v_{31}-3)\qquad(n_l=3).
$$

The large-$N_c$ rules in [@Herrera-Siklody:1996tqr] give $v_{12}=O(N_c^{-2})$ and $v_{31}=O(N_c^{-1})$. Hence $B_{47}(0)=O(N_c^{-2})$, beyond the paper's retained $O(1/N_c)$ correction to the $O(p^4)$ counterterms. Under that same truncation used for the other quoted $\Gamma_i$, the result is therefore

$$
\Gamma_{47}=0\quad\text{at the retained order of [@Herrera-Siklody:1996tqr].}
$$

This zero is a consequence of that specified truncation, not a general identity: at the next $1/N_c$ order, the general $B_{47}(X)$ need not vanish. The operator itself is present in [the singlet-source derivative basis](singlet-source-derivative-basis.md), so it must be kept with its matched counterterm if that higher order is retained.

## Constants entering the retained quadratic Lagrangian

The quadratic expansion already derived in the [quadratic-lagrangian](quadratic-lagrangian.md) contains the following local constants:

| Sector | Constants in the local quadratic terms | Role |
|---|---|---|
| Two-meson kinetic | $L_4^r,L_5^r,\Lambda_1^r,L_{18}^r$ | Octet and singlet-trace kinetic normalization and mixing |
| Two-meson mass | $L_6^r,L_7^r,L_8^r,\Lambda_2^r,L_{25}^r$ | Quark-mass and singlet/topological mass entries |
| Meson--ALP kinetic | $L_4^r,L_5^r,\Lambda_1^r,L_{18}^r,L_{47}^r$ | Derivative mixing entries |
| Meson--ALP mass | $L_6^r,L_7^r,L_8^r,\Lambda_2^r,L_{25}^r$ | Mass mixing entries |

Here $L_4^r,L_5^r,L_6^r,L_7^r,L_8^r,L_{18}^r,L_{25}^r$ use the $U(3)$ subtraction stated above; $\Lambda_{1,2}$ have no chiral-loop pole at this order, and $L_{47}$ has the function-valued counterterm $B_{47}(X)$. The ALP--meson terms are proportional to the ALP sources already shown in the [quadratic-lagrangian](quadratic-lagrangian.md) and [singlet-source derivative basis](singlet-source-derivative-basis.md). There is no claim that every listed constant has nonzero contribution for every source choice or flavor channel.

## Relation to the current quadratic record

Replacing the bare symbols by renormalized symbols in the existing local expressions is necessary but not sufficient for an NLO two-point function. The complete renormalized inverse propagator has the form

$$
\mathcal Q_{IJ}(s)=sK_{IJ}^{\rm LO}-M_{IJ,\rm LO}^2
+s\,\delta K_{IJ}^{\rm local}(\mu)-\delta M_{IJ,\rm local}^2(\mu)
+\Sigma_{IJ}^{\rm loop}(s;\mu),
$$

where $I,J$ range over the mesons and, when retained, the ALP. The sign convention for $\Sigma$ is defined by this displayed inverse-propagator equation. The renormalized local terms and loop self-energy are separately scale dependent; their sum is scale independent to the calculated order. The loop matrix includes both diagonal meson self-energies and meson--ALP mixed two-point functions if ALP-dependent vertices are included. Local polynomial pieces may be shifted between $\Sigma$ and the local constants by a finite scheme change, while the nonanalytic thresholds remain.

This note records the scope and renormalization structure, not explicit $\Sigma_{IJ}^{\rm loop}$ formulas. Such formulas require a specified ALP-source specialization and a diagrammatic calculation in the corresponding $U(3)$ scheme.

See also the declared power counting, loop classes, subtraction scheme and counterterm basis in [counting/INDEX.md](counting/INDEX.md).
