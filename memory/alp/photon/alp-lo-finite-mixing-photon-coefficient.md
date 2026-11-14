---
title: LO coefficient of a -> gamma gamma at finite ALP-meson mixing in dressed invariants
status: derived
created: 2026-11-14
updated: 2026-11-14
bibtex_keys: [Alda:2026nrd, Bauer:2020jbp]
provenance: analytical-derivation
---

# LO coefficient of $a\to\gamma\gamma$ at finite ALP-meson mixing in dressed invariants

## Inputs

- Photon coupling of the [dressed photon record](alp-dressed-photon-invariant.md): $\widehat\Theta_{\gamma,1}=c_\gamma+2N_c\operatorname{tr}(Q^2c_q^A)$, the first-order dressed photon coefficient, with $Q=\operatorname{diag}(2/3,-1/3,-1/3)$, $N_c=3$, and the normalization $\mathcal L\supset\frac{\alpha_{\rm em}}{4\pi}\big[c_\gamma\frac{a}{f_a}+\frac{N_c}{f_\pi}\operatorname{tr}(Q^2\Phi)\big]F\widetilde F$ of the [photon invariant record](alp-photon-invariant-from-mixing.md) ([@Bauer:2020jbp], [@Alda:2026nrd]).
- First-order dressed mass and gluonic coefficients of the [general redefinition record](../shared/alp-general-chiral-redefinition-invariants.md) and the [first-order record](../shared/alp-dressed-mass-first-order-invariants.md): $\widehat\Theta_1=c_{GA}$, $\widehat P_1=-\{m_q,c_q^{AP}\}$ with $c_q^{AP}=c_q^A-c_q^P$, and $m_q=M(0)$.
- Tadpole source of the [exact pole-equation record](../shared/alp-lo-exact-pole-equations.md), written in dressed form,
$$
\mathcal C=B_0\widehat P_1-\frac{m_0^2}{N_f}\widehat\Theta_1\mathbf 1
=-B_0\{m_q,c_q^{AP}\}-\frac{m_0^2}{N_f}c_{GA}\mathbf 1 .
$$
- Mass operator $\mathbb M_s(X)=\{B_0m_q,X\}+\frac{m_0^2}{N_f}\operatorname{tr}(X)\mathbf 1-sX$ on the neutral $U(3)$ space, $s$ the squared mass, with inverse $\mathbb M_s^{-1}$ defined by $\mathbb M_s(\mathbb M_s^{-1}(X))=X$.
- Physical ALP pole $m_a^2$: the exact root of the pole equation, not the bare $m_{a,0}^2$.

## Assumptions

- LO chiral Lagrangian, tree level, stable simple poles, no ALP–meson degeneracy at $m_a^2$.
- Redefinitions commute with $Q$ (the class of the [photon record's assumption 2](alp-dressed-photon-invariant.md)); see Validity.
- The photon vertex has no direct $c_q^A$ term in the original fields, consistent with the [photon invariant record](alp-photon-invariant-from-mixing.md) (the $c_q^A$ dependence enters through the kinetic mixing).

## Derivation

**Transformation that removes the ALP–meson mixing at the physical pole.** Define $\epsilon=f_\pi/f_a$. The exact alignment of the ALP pole is the field change
$$
P_i=P_i'+w_i\,a,\qquad
w_i=\epsilon\operatorname{tr}\Big(\lambda_i\Big[c_q^A+\mathbb M_{m_a^2}^{-1}\big(-B_0\{m_q,c_q^{AP}\}-\tfrac{m_0^2}{N_f}c_{GA}\mathbf 1\big)\Big]\Big),
$$
equivalently $\sum_iw_i\lambda_i=2\epsilon\big[c_q^A+\mathbb M_{m_a^2}^{-1}(\mathcal C)\big]$. The ALP pole has zero meson component in the primed frame. In the original basis the same vector is the pole eigenvector $(w_i,1)$, with the ALP normalization below. This relation was checked against the spectrum solver at the ALP pole to $10^{-15}$ at $\epsilon=1/2$ and $\epsilon=1/10$.

**Vertex in the primed frame.** Substituting $P_i=P_i'+w_ia$ into $\operatorname{tr}(Q^2\Phi)$ gives $\operatorname{tr}\big(Q^2\sum_iw_i\lambda_i\big)a=2\epsilon\operatorname{tr}\big(Q^2[c_q^A+\mathbb M_{m_a^2}^{-1}(\mathcal C)]\big)a$. The $a$ coefficient of the photon Lagrangian is therefore
$$
c_\gamma+2N_c\operatorname{tr}\big(Q^2c_q^A\big)+2N_c\operatorname{tr}\Big(Q^2\,\mathbb M_{m_a^2}^{-1}(\mathcal C)\Big)
=\widehat\Theta_{\gamma,1}+2N_c\operatorname{tr}\Big(Q^2\,\mathbb M_{m_a^2}^{-1}(\mathcal C)\Big).
$$

**Physical amplitude.** The external ALP leg carries the normalization of its pole, $1/\sqrt{1+2\epsilon^2\operatorname{tr}[\mathbb M_{m_a^2}^{-1}(\mathcal C)^2]}$, which equals $\sqrt{K'_{aa}}$ for the aligned field. The amplitude is
$$
\boxed{\mathcal M(a\to\gamma\gamma)=-\frac{e^2}{4\pi^2f_a}\,
\frac{\widehat\Theta_{\gamma,1}+2N_c\operatorname{tr}\!\big[Q^2\,\mathbb M_{m_a^2}^{-1}(\mathcal C)\big]}
{\sqrt{1+2\epsilon^2\operatorname{tr}\big[\big(\mathbb M_{m_a^2}^{-1}(\mathcal C)\big)^2\big]}}\;
\epsilon_{\mu\nu\rho\sigma}q_1^\mu\varepsilon_1^\nu q_2^\rho\varepsilon_2^\sigma}
$$
with $\epsilon=f_\pi/f_a$. The sign and factor follow the convention $i\mathcal M=-ie^2EC/(4\pi^2f_a)$ of the photon invariant record.

**Invariance.** Every ingredient is a dressed invariant or a physical mass: $\widehat\Theta_{\gamma,1}$, $\widehat P_1$, $\widehat\Theta_1$, $m_q$, $m_a$. So the amplitude is unchanged by the redefinitions of the Q-commuting class. The second-order coefficients of the dressed sources enter only the ALP mass entry, so for one ALP they enter only through the physical $m_a$. This is a consequence of the single-ALP vertex structure, not an assumption.

**Weak-mixing limit.** As $\epsilon\to0$, $m_a^2\to m_{a,0}^2$ and the normalization tends to 1. The numerator becomes $C_{a\gamma\gamma}(m_a)=\widehat\Theta_{\gamma,1}+2N_c\operatorname{tr}[Q^2\mathbb M_{m_a^2}^{-1}(\mathcal C)]$, the first-order coefficient of the [photon invariant record](alp-photon-invariant-from-mixing.md) and [dressed photon record](alp-dressed-photon-invariant.md). The first correction is of order $\epsilon^2$.

## Checks

From [`scripts/derive_alp_lo_photon_channels.py`](../../../scripts/derive_alp_lo_photon_channels.py), tested in [`scripts/test_derive_alp_lo_photon_channels.py`](../../../scripts/test_derive_alp_lo_photon_channels.py) (18 tests):

1. Canonical frame: the change $P_i=P_i'+\epsilon\operatorname{tr}(\lambda_ic_q^A)a$ gives kinetic matrix $\mathbf 1$ and mixed entry $-\epsilon\operatorname{tr}\big(\lambda_i\mathcal C\big)$, exact for symbolic $\epsilon$.
2. Normalization identity, exact for symbolic $s$ and $\epsilon$: $(w,1)^{\mathsf T}K(w,1)=1+2\epsilon^2\operatorname{tr}[(\mathbb M_s^{-1}(\mathcal C))^2]$, which equals $dD^{\rm eff}_{aa}/ds$ at the pole. Agreement with the solver to $10^{-15}$.
3. Weak-mixing expansion: the $\epsilon\to0$ limit equals the photon invariant record's coefficient, and the remainder is exactly of order $\epsilon^2$ (checked at $\epsilon=10^{-2},10^{-3},10^{-4}$).
4. Invariance: under constant axial $\kappa$, exact; under random nonlinear $g_{R,L}(x)=e^{i(xG_1+x^2G_2)}$ commuting with $Q$, the dressed invariants and the numerator are unchanged to $10^{-15}$ at the physical pole.
5. Ward identity: contraction with $\varepsilon_\gamma\to$ the photon momentum vanishes for the $\gamma\gamma$ structure.

## Reproducibility

```sh
.venv/bin/python scripts/derive_alp_lo_photon_channels.py
.venv/bin/python -m unittest discover -s scripts -p test_derive_alp_lo_photon_channels.py
```

Both ran on 2026-11-14 and passed.

## Validity and open points

- Redefinitions that do not commute with $Q$ turn the photon coupling into an $x$-dependent chiral charge. The records do not treat this case, so the invariance statement is limited to the Q-commuting class.
- The sign of the $\gamma\gamma$ WZW contribution is not tied to $\mathrm{orientation}=-1$ in any record. This is not checked here, and the $\gamma\gamma$ coefficient above is taken from the photon invariant record.
- The leg construction used in the $\gamma PP$ record is singular for a neutral pole with zero ALP component, such as the exact $K^0$ pole; see [the $\gamma PP$ record](../photon-two-pseudoscalar/alp-lo-finite-mixing-wzw-pp-amplitude.md).
- The LO model is uncontrolled where the ALP admixture is $O(1)$ or near-degenerate ([truncation record](../shared/alp-lo-truncation-limits.md)). Vector-meson conversion ([vector-mixing record](alp-photon-vector-mixing.md)) and resonance terms are not included.
