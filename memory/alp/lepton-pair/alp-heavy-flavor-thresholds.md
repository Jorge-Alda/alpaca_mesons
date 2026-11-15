---
title: Charm and bottom thresholds in the short-distance coefficient of the ALP form factor
status: mixed
created: 2026-11-15
updated: 2026-11-15
bibtex_keys: [Chetyrkin:1997un, Chetyrkin:1997sg, Larin:1993tq, Larin:1994va, Ludtke:2020moa, ParticleDataGroup:2024cfk]
provenance: analytical-derivation
---

# Scope, inputs, and assumptions

This record treats the charm and bottom quarks in the short-distance expansion of the doubly virtual ALP
transition form factor, for photon virtualities that are integrated up to infinity. It uses the
[leading hard-current expansion](alp-hard-current-ope.md), the
[axial-current matrix element](alp-axial-current-matrix-element.md) and the
[singlet evolution](alp-singlet-axial-current-evolution.md). Here $q^2<0$ is the photon virtuality,
$\mu=\sqrt{-q^2}$ is the scale of the Wilson coefficient, and $\langle\cdot\rangle$ is the flavor trace.

Inputs and provenance:

- Matching of the strong coupling at a quark threshold, in the $\overline{\rm MS}$ scheme with
  $\mu_h=m_h(m_h)$ and $\ell_h=\ln(\mu^2/\mu_h^2)$, with $\alpha_s\equiv\alpha_s^{(n_f)}$ in the $n_f$-flavor theory:
  $\alpha_s^{(n_f-1)}/\alpha_s^{(n_f)}=1-\frac{\alpha_s}{\pi}\frac{\ell_h}{6}+\bigl(\frac{\alpha_s}{\pi}\bigr)^2\bigl(\frac{\ell_h^2}{36}-\frac{19}{24}\ell_h+\frac{11}{72}\bigr)+\dots$
  [@Chetyrkin:1997sg, see the [extract](../../extracts/Chetyrkin:1997sg/alphas-threshold-matching.md);
  @Chetyrkin:1997un, see the [extract](../../extracts/Chetyrkin:1997un/alphas-decoupling-constant.md)].
- Beta-function coefficients $\beta_0^{(n_f)}=\frac14[11-\frac23n_f]$, $\beta_1^{(n_f)}=\frac1{16}[102-\frac{38}3n_f]$ for
  $\mu^2\,d(\alpha_s/\pi)/d\mu^2=-\beta_0(\alpha_s/\pi)^2-\beta_1(\alpha_s/\pi)^3$ [@Chetyrkin:1997un].
- Singlet axial-current anomalous dimension and evolution [@Larin:1993tq], as in the
  [singlet evolution record](alp-singlet-axial-current-evolution.md). Whether a heavy quark decouples from the singlet
  axial part is not established for this process; the only statement available, for a weak-current correlator, is that
  it need not [@Larin:1994va].
- NLO coefficient $1-\alpha_s(\mu)/\pi$ of the leading local operator [@Ludtke:2020moa].
- Masses [@ParticleDataGroup:2024cfk, see the [extract](../../extracts/ParticleDataGroup:2024cfk/quark-and-tau-masses.md)]:
  $m_c(m_c)=1.2730\pm0.0046$ GeV, $m_b(m_b)=4.183\pm0.007$ GeV, $m_\tau=1776.93\pm0.09$ MeV,
  $m_e=0.51099895$ MeV, and the top-quark mass $172.57\pm0.29$ GeV from direct measurements (a kinematic reconstruction,
  usually interpreted as the pole mass).

Assumptions (not results): the heavy quark is treated as a free quark in the hard propagator at leading order in
$\alpha_s$; thresholds are placed at $\mu_c=m_c(m_c)$ and $\mu_b=m_b(m_b)$; one-loop couplings are continuous at the thresholds
(the constant $\tfrac{11}{72}$ is neglected in the evolution); the singlet matrix element is continuous across each threshold.

# Heavy-quark mass in the local coefficient

The [script](../../../scripts/derive_alp_heavy_flavor_threshold.py) repeats the Wick and Clifford derivation of the hard-current
record with the massive propagator $S_F(k)=i(\slashed k+m)/(k^2-m^2+i0)$ and projects the result by traces onto the 16 Clifford
basis matrices. For one flavor of charge $Q_h$ and mass $m_h$,
$$
\boxed{i\!\int\! d^4x\,e^{iq\cdot x}T\{J^\mu_{\rm em}(x/2)J^\nu_{\rm em}(-x/2)\}\supset
\frac{2i\,Q_h^2}{q^2-m_h^2+i0}\,\epsilon^{\mu\nu\rho\sigma}q_\rho\,\bar h\gamma_\sigma\gamma_5h
-\frac{2m_h Q_h^2}{q^2-m_h^2+i0}\,g^{\mu\nu}\,\bar hh .}
$$
The terms linear in $m_h$ combine to $m_h\{\gamma^\mu,\gamma^\nu\}=2m_hg^{\mu\nu}$; their antisymmetric part vanishes,
and there is no pseudoscalar, vector, or $\sigma^{\rho\sigma}$ term. The second term has no $\epsilon$ tensor and therefore does not
contribute to the scalar $F_{a\gamma^*\gamma^*}$. The axial term reduces to the massless result at $m_h=0$. Relative to the
massless coefficient its factor is
$$
\frac{q^2}{q^2-m_h^2}=\frac{|q^2|}{|q^2|+m_h^2}\;\longrightarrow\;
\begin{cases}1-\dfrac{m_h^2}{|q^2|}+\dfrac{m_h^4}{|q^2|^2}-\dots,&|q^2|\gg m_h^2,\\[2mm]
\dfrac{|q^2|}{m_h^2}-\dfrac{|q^2|^2}{m_h^4}+\dots,&|q^2|\ll m_h^2.\end{cases}
$$
The limit $|q^2|\ll m_h^2$ is only the free-quark expression: there the heavy quark is not a light field of the
low-energy theory, and the local expansion in terms of $\bar h\gamma_\sigma\gamma_5h$ is not a controlled description.

**Heavy-flavor axial matrix element is a named input.** The contribution of $h=c,b$ to the operator term is
$$
F_{a\gamma^*\gamma^*}\supset\frac{2i\,Q_h^2}{(q^2-m_h^2)\,m_a^2}\,p_\sigma\langle0|\bar h\gamma^\sigma\gamma_5h|a(p)\rangle ,
$$
and $\langle0|\bar h\gamma^\sigma\gamma_5h|a(p)\rangle$ is not supplied by the three-flavor chiral theory of the
[matrix-element record](alp-axial-current-matrix-element.md), where only $u,d,s$ appear. No source read gives it. The
ALP couplings to heavy quarks, and the term independent of $q^2$ that a derivative coupling to a heavy quark produces in the hard
loop, belong to task 3.1.

# Strong-coupling matching at the thresholds

The script verifies that the logarithmic coefficients of the matching relation follow from the renormalization-group equations of the
$n_f$ and $n_f-1$ flavor theories with $\mu_h$ independent of $\mu$: the coefficient at order $\alpha_s$ is
$\beta_0^{(n_f)}-\beta_0^{(n_f-1)}=-\tfrac16$, the coefficient of $\ell_h^2$ at order $\alpha_s^2$ is $\tfrac1{36}$, and the coefficient of $\ell_h$ at order $\alpha_s^2$ is $-\tfrac{19}{24}$ when the constant at order $\alpha_s$ is zero. All agree exactly with the cited relation. The constant
$\tfrac{11}{72}$ is not fixed by this condition. At $\mu=\mu_h$ there is no constant at order $\alpha_s$, so the coupling is continuous at
one loop, and
$$
\alpha_s^{(n_f-1)}-\alpha_s^{(n_f)}=\tfrac{11}{72}\,\alpha_s\Bigl(\frac{\alpha_s}{\pi}\Bigr)^2,\qquad
\Bigl[1-\frac{\alpha_s}{\pi}\Bigr]_{\text{below}}-\Bigl[1-\frac{\alpha_s}{\pi}\Bigr]_{\text{above}}=-\frac{11}{72}\Bigl(\frac{\alpha_s}{\pi}\Bigr)^3 .
$$
The jump of the NLO coefficient is therefore two orders beyond its accuracy; the coefficient is continuous at the order kept.
The script also checks that the form with the running mass $m_h(\mu)$ in [@Chetyrkin:1997un] (coefficient $-\tfrac{11}{24}$) agrees with
the form above only for $\mu^2\,dm_h/d\mu^2=-(\alpha_s/\pi)\,m_h$ at one loop; this is the standard value but it is not
recorded from a source in this repository.

# Singlet evolution across the thresholds

With $\mu_0<\mu_c$ and $n=3,4,5$ active flavors in the three intervals, the logarithm of the ratio of the singlet axial-current
matrix element at $\mu>\mu_b$ to that at $\mu_0$ is
$$
\boxed{\ln\frac{\langle0|\sum_f\bar q_f\gamma^\mu\gamma_5q_f|a(p)\rangle_\mu}{\langle0|\sum_f\bar q_f\gamma^\mu\gamma_5q_f|a(p)\rangle_{\mu_0}}
=\frac2\pi\Bigl[\tfrac13\bigl(\alpha_s(\mu_c)-\alpha_s(\mu_0)\bigr)+\tfrac{12}{25}\bigl(\alpha_s(\mu_b)-\alpha_s(\mu_c)\bigr)+\tfrac{15}{23}\bigl(\alpha_s(\mu)-\alpha_s(\mu_b)\bigr)\Bigr],}
$$
where in each interval $\alpha_s$ is the one-loop coupling of the theory with the corresponding number of flavors, and for lower
$\mu$ the sum stops at that scale. The coefficient in front of each bracket is $2n/(\pi\beta_0^{(n)})$ with $\beta_0^{(n)}=11-\frac23n$ in the form used for the singlet evolution. The script checks exactly that the evolution equation
$\mu^2\,d\ln(\text{ratio})/d\mu^2=-6C_Fn\,(\alpha_s/4\pi)^2$ holds in each interval, that the coupling and the ratio are
continuous at both thresholds, and that the interval with $n=3$ reproduces the fixed-flavor result of the singlet record. The derivative of the logarithm
jumps at each threshold by $-6C_F(\alpha_s/4\pi)^2=-\alpha_s^2/(2\pi^2)$ per added flavor, because the anomalous dimension
depends on $n$.

**Continuity of the matrix element across a threshold is an assumption.** No matching coefficient for the singlet axial current is
available in the sources read, and [@Larin:1994va] states non-decoupling in a different process. A finite matching constant would multiply the ratio at
each threshold. The four-flavor and five-flavor regimes are flagged accordingly.

# Active flavors for $2m_e<m_a<2m_\tau$

With the PDG values above, $2m_e=1.022$ MeV and $2m_\tau=3.554$ GeV, the ALP-mass range lies below $\mu_b=4.183$ GeV and spans
$\mu_c=1.273$ GeV. The number of active flavors is a function of $\mu=\sqrt{-q^2}$ (and therefore of the
matching scale):

| Scale $\mu=\sqrt{-q^2}$ | Active flavors | Comment |
|---|---|---|
| $\mu<m_c(m_c)=1.273$ GeV | 3 ($u,d,s$) | scale of the low-energy matrix element, $\mu_0\sim1$ GeV |
| $1.273\ \text{GeV}<\mu<m_b(m_b)=4.183$ GeV | 4 | includes the whole range of $m_a$ above $1.273$ GeV |
| $4.183\ \text{GeV}<\mu<m_t$ | 5 | reached only by the virtuality integral, not by $m_a$ |
| $\mu\gtrsim m_t$ ($\approx172.6$ GeV, direct measurements) | outside the supported domain | electroweak effects are not treated |

For the matching scale, which must not lie below $m_a$: if $m_a\gtrsim1.273$ GeV (up to $2m_\tau=3.554$ GeV) the matching scale lies in the
four-flavor region, so charm must be included as an active flavor there and the supported ALP-mass range cannot use three massless
flavors at the matching point. For $m_a<1.273$ GeV the matching scale can lie in either region. The marker $m_t$ above is
nominal and has not been derived; where the top threshold and electroweak corrections limit the supported domain is a decision that remains open.

# Remaining named inputs

1. The heavy-quark matrix elements $\langle0|\bar h\gamma^\sigma\gamma_5h|a(p)\rangle$ for $h=c,b$, and the ALP couplings to $c$ and $b$.
2. The matching coefficient of the singlet axial current at the thresholds.
3. The identification of the free-quark mass $m_h$ in the hard propagator with $m_h(m_h)$ at the order beyond leading.
4. The upper limit of the supported virtuality domain (top threshold and electroweak effects).

# Verification

Run
```text
.venv/bin/python scripts/derive_alp_heavy_flavor_threshold.py
.venv/bin/python scripts/derive_alp_heavy_flavor_threshold.py --test
```
Fourteen unit tests pass. They cover the massive Dirac-matrix projection (including the comparison with the massless script), the
heavy-flavor factor and its expansions, the renormalization-group consistency of the matching logarithms, the size of the
jump, the piecewise singlet evolution and its continuity, the jumps of the derivatives, and the printed table of active flavors with symbolic thresholds.
Numerical values of the masses enter only the table above, from the cited extract.
