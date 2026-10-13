---
title: Normalized SM radiative amplitude with complete neutral-vector exchange
status: derived
created: 2026-10-13
updated: 2026-10-13
bibtex_keys: [Alda:2026nrd, Ruiz-Femenia:2003jdx, Ecker:1988te, Ecker:1989yg]
provenance: analytical-derivation
---

# SM radiative amplitude reconstructed from the action

## Inputs and assumptions

Use the VJP/VVP operators and normalizations from
[the operator digest](../../extracts/Ruiz-Femenia:2003jdx/vector-odd-parity.md)
and [the VVP digest](../../extracts/Ruiz-Femenia:2003jdx/vvp-three-pseudoscalar.md),
the even-parity mixing and canonical vector kinetic action in
[the foundations](../foundations/rcht-basics.md), and the already executed
[two-photon tensor derivation](radiative-tff-amplitude.md).
The same Fourier contractions produce the one-photon numerators here.

Tree level, canonical external vector residue, single-trace RChT action,
single vector nonet and common $F_V$ are assumptions. The internal pole masses
may be explicit inputs; a common $M_V$ is the baseline in the cited paper.
The source's octet action does not establish extra singlet/OZI matching or a
fully matched action for arbitrary split internal masses. No fit is performed.
Metric $(+---)$, source $\epsilon_{0123}=+1$, natural units, masses in GeV.

### Field conventions and Alda map

The canonical RChT vector octet is $V=\lambda^aV^a/\sqrt2$, with
$\operatorname{tr}(\lambda^a\lambda^b)=2\delta^{ab}$. Consequently its eighth
generator is $\operatorname{diag}(1,1,-2)/\sqrt6$. Ideal neutral vector matrices
are $T_\rho=\operatorname{diag}(1,-1,0)/\sqrt2$,
$T_\omega=\operatorname{diag}(1,1,0)/\sqrt2$,
$T_\phi=\operatorname{diag}(0,0,1)$. These already exist in
`odd_parity_rcht.ideal_vector_generators()`.

Pseudoscalar coefficient matrices remain explicit inputs and are not assigned
vector-state defaults.

The [indexed Alda appendix](../../extracts/Alda:2026nrd/appendix-a-expansion.md)
defines $\phi=\Phi_A/(2f_\pi)$; RChT uses
$u=\exp(i\Phi_R/(\sqrt2F))$. Therefore the explicit field map is
$\Phi_R=F\Phi_A/(\sqrt2f_\pi)$, and $u_\mu=-\partial_\mu\Phi_A/f_\pi$
at one-field order. The implementation requires explicit coefficients of
$\Phi_R$; Alda matrices must first use this map. No vector normalization is
inferred from Alda's pseudoscalar convention. Physical $\eta,\eta'$ residues,
decay constants, mixing and mass-spurion values remain explicit inputs.

## Derivation and flavor sums

Reproduce with

```text
.venv/bin/python scripts/derive_radiative_sm_amplitude.py
.venv/bin/python scripts/diagnose_radiative_rank.py
```

The executable uses the actual all-incoming Lorentz contractions of the
VJP/VVP action, including every Fourier derivative and both vector assignments.
For external vector matrix $T_V$, pseudoscalar coefficient $T_P$, and
$\chi=2B_0m_q$ define

$$
F_Q=\operatorname{tr}(\{T_V,Q\}T_P),\quad
F_H=\operatorname{tr}(\{T_V,Q\}\{T_P,\chi\}),\quad
F_K=\operatorname{tr}\big(T_V[[T_P,Q],\chi]\big).
$$

Matrices include the conjugation appropriate to each external leg; for
charged $V^+\to P^+\gamma$, use $T_V=E_{ji}$, $T_P=E_{ij}$. Crossing to
$P\to V\gamma$ retains the vector and pseudoscalar identities and masses.
The vector mass is $M_{\rm ext}$ in both decay orientations.

The linear electromagnetic building block is derived by expanding
$u=1+i\Phi_R/(\sqrt2F)$:

$$
f_-^{\mu\nu}=-\frac{\sqrt2 i e}{F}[\Phi_R,Q]F_\gamma^{\mu\nu}
+O(\Phi_R^2),\qquad \chi_+=2\chi+O(\Phi_R^2).
$$

Thus $\mathcal O^4_{VJP}$ contributes at one-pseudoscalar order. Its trace
vanishes for diagonal neutral $T_P$ and $K^0$ (equal $d,s$ charges), but
$F_K=-(\chi_s-\chi_u)$ for $K^{*+}\to K^+\gamma$. This mode accesses $c_4$.
For charged pion modes the term is proportional to isospin breaking.

For each canonical neutral internal vector $T_\alpha$, set

$$
w_\alpha=\operatorname{tr}(T_\alpha Q),\quad
F_{Q,\alpha}=\operatorname{tr}(\{T_V,T_\alpha\}T_P),\quad
F_{H,\alpha}=\operatorname{tr}(\{T_V,T_\alpha\}\{T_P,\chi\}).
$$

The photon-induced tensor response requires the **sum** over neutral vectors,
$\sum_\alpha w_\alpha T_\alpha/M_\alpha^2$ at $k_\gamma^2=0$. For common
mass, the script proves $\sum_\alpha w_\alpha T_\alpha=Q$, so the summed
exchange traces equal the direct charge/mass traces. A single rho does not
represent this response: $\rho^0\to\pi^0\gamma$ has zero rho charge trace but
a nonzero omega contribution $w_\omega F_{Q,\omega}=1/3$.

## Normalized result

For the source phase convention
$i\mathcal M=i\epsilon_{\alpha\beta\rho\sigma}
\epsilon_V^\alpha\epsilon_\gamma^\beta q^\rho p_P^\sigma G$, define

$$
\begin{aligned}
C_V&=F_Q\big[(c_1+c_2-c_5)m_P^2
 +(-c_1+c_2+c_5-2c_6)M_{\rm ext}^2\big]+4c_3F_H-4c_4F_K,\\
H_\alpha&=F_{Q,\alpha}\big[(d_1-d_3)m_P^2+d_3M_{\rm ext}^2\big]
 +4d_2F_{H,\alpha}.
\end{aligned}
$$

The action gives

$$
\boxed{G=\frac{2\sqrt2e}{FM_{\rm ext}M_V}C_V
 -\frac{4eF_V}{FM_{\rm ext}}\sum_\alpha
 \frac{w_\alpha H_\alpha}{M_\alpha^2}.}
$$

$G$ has dimension GeV$^{-1}$; $M_V$ normalizes the operators and
$M_{\rm ext}$ is always the external **vector** mass, even in $P\to V\gamma$.
$c_7$ and $d_4$ vanish for these tensor wavefunctions. The derivation verifies
the published omega-pion Eq. (27) exactly: canonical states give $F_Q=1$,
$F_H=2m_\pi^2$, $F_K=0$ at $\chi_u=\chi_d=m_\pi^2$, and the complete
exchange reduces to its rho contribution. This fixes the absolute normalization
without a calibration fit. Both VVP assignments place $d_1$ only in the
$m_P^2$ coefficient.

Spin sums and initial spin averages give

$$
\Gamma(V\to P\gamma)=\frac{|G|^2(M_{\rm ext}^2-m_P^2)^3}{96\pi M_{\rm ext}^3},
\qquad
\Gamma(P\to V\gamma)=\frac{|G|^2(m_P^2-M_{\rm ext}^2)^3}{32\pi m_P^3}.
$$

The numerical implementation is
[`radiative_sm.py`](../../../scripts/implementation/radiative_sm.py).
It takes explicit flavor matrices and real diagonal $\chi$ in GeV$^2$;
canonical ideal internal vectors and their three pole masses are explicit
model assumptions. It supplies both width orientations and rejects closed
channels. It does not replace the existing source-specific omega predictor.

## Rank experiment and limits

The ten-channel executable rank experiment uses canonical vector matrices,
explicit exploratory $\eta,\eta'$ coefficient matrices and reference masses,
and the stated LO isospin spurion
$\chi_u=\chi_d=m_\pi^2$, $\chi_s=2m_{K^+}^2-m_\pi^2$. These inputs do not
establish physical eta residues or the precision of LO quark-mass relations.
No data error or fit enters the experiment.

* With a common internal pole mass, OPE plus action-based radiative rows has
  rank **7 of 11**. The new directions are a correlated $c_3,d_2$ combination
  from mass-spurion mismatch and $c_4$ from the charged kaon.
* With explicitly split $(M_\rho,M_\omega,M_\phi)$ input poles, rank is
  **8 of 11**. The eighth unit-row singular value is $0.00357049$, while
  the three remaining null singular values are below $3\times10^{-16}$.
  It demonstrates algebraic accessibility; its small size supplies no
  evidence of a statistically useful constraint. Pole splitting here is a
  phenomenological extension whose accompanying symmetry-breaking operators
  and matching have not been established.
* In both runs, the exact unprobed directions are $c_7$, $d_4$ and
  $\delta c_2=\delta c_5=\delta c_6$, with all other entries unchanged.
  The common-pole case has one further null direction because its direct and
  exchange mass traces are proportional. Splittings can change that relative
  response; they are not an overall rescaling of each row.
* Adding the full common-pole neutral TFF rows at $Q^2=5$ GeV$^2$ with
  the **same independent** exploratory pseudoscalar matrices and $\chi$
  raises the common-pole radiative rank $7\to8$, and preserves the split-pole
  rank 8. This is sensitivity to a second $c_3,d_2$ combination when
  $F_H\ne2m_P^2F_Q$. The executed probe includes this comparison.
  If every TFF row instead assumes $F_H=2m_P^2F_Q$, it remains in the OPE
  span and adds zero rank, as established in the TFF derivation.

The next requirements for a physical coupled fit are source-grounded
pseudoscalar field/decay-constant inputs and a declared common-pole versus
consistently broken-mass action. The experimental covariance, PDG scale
factors and eta-prime resonant/nonresonant composition remain independent
limitations in the existing dataset plan. No saturation conclusion, fitted
coefficient, branching-fraction prediction or $\chi^2$ is reported.
