---
title: Charm and bottom contribution to the ALP two-photon form factor from their ALP couplings
status: derived
created: 2026-11-15
updated: 2026-11-15
bibtex_keys: [Alda:2026nrd, Messerli:2025rnv, Bauer:2021mvw, Alda:2024cxn]
provenance: analytical-derivation
---

# Scope, inputs, and assumptions

This record gives the contribution of a heavy quark $h=c,b$ to $F_{a\gamma^*\gamma^*}(q_1^2,q_2^2)$ at one loop and leading order in
$\alpha_s$, in terms of its ALP couplings. It replaces the vacuum-to-ALP matrix element $\langle0|\bar h\gamma^\sigma\gamma_5h|a(p)\rangle$,
which the three-flavor chiral theory does not supply. $F_h$ is the contribution of the flavor $h$ to the form factor; its place in the total is set out in the [assembly record](alp-high-energy-tff-assembly.md).
Conventions: [shared conventions](alp-high-energy-tff-conventions.md); $i\mathcal M=ie^2RF_{a\gamma^*\gamma^*}$
[@Messerli:2025rnv, Eq. (9)]; photon momenta $q_1,q_2$ with $p=q_1+q_2$, $p^2=m_a^2$; $N_c=3$; $Q_h$ the charge of $h$.

ALP couplings, extending the partonic Lagrangian of [@Alda:2026nrd, Eq. (lagr\_partons)] to $h$:
$$
\mathcal L\supset\frac{\partial_\mu a}{f_a}\,\bar h\gamma^\mu\bigl(c_{hh}^V+c_{hh}^A\gamma_5\bigr)h
+i\frac{a}{f_a}\,2m_h\,c_{hh}^P\,\bar h\gamma_5h .
$$
The sign of the pseudoscalar term is that of [@Alda:2026nrd] (the diagonal case of $+i\frac af_a\bar q\{c_q^P,m_q\}\gamma_5q$).
The vector coupling does not contribute (executed in the script, see the checks below).
Assumptions: free quark loop (no gluon exchange, corrections of relative order $\alpha_s$); the mass $m_h$ is a single parameter whose
scheme is not specified at this order. The closed forms below are real for $m_a<2m_h$; the continuation above threshold is given in its own section.

# Result

The [script](../../../scripts/derive_alp_heavy_quark_triangle.py) computes the pseudoscalar triangle from Dirac traces and Feynman
parameters: both photon orderings have the numerator trace $4i\,m_h\,\epsilon^{\mu\nu\rho\sigma}q_{1\rho}q_{2\sigma}$, independent of the loop momentum.
For the derivative coupling it applies the axial Ward identity with a Pauli-Villars regulator of mass $M$: each regulated piece obeys the naive identity, so the
amplitude is $[2m_h\Gamma_P(m_h)-2M\Gamma_P(M)]$, and the $M\to\infty$ limit of the second piece is the anomaly term. The sum, with the overall sign fixed by the direct trace calculation, is expressed with the [loop function $B_4$](alp-loop-function-b4.md),
$\tau_a=4m_h^2/m_a^2$, $\tau_i=4m_h^2/q_i^2$:
$$
\boxed{F_h=\frac{2N_cQ_h^2}{4\pi^2f_a}\Bigl[c_{hh}^P+\bigl(c_{hh}^A-c_{hh}^P\bigr)B_4(\tau_a,\tau_1,\tau_2)\Bigr],}
$$
$$
B_4(\tau_a,\tau_1,\tau_2)=1-\int_0^1\!dx\int_0^{1-x}\!dy\;\frac{2m_h^2}{m_h^2-xy\,m_a^2-x(1-x-y)\,q_2^2-y(1-x-y)\,q_1^2-i0}.
$$
Here $x$ labels the quark line between the photon-2 vertex and the ALP vertex (paired with $q_2^2$), $y$ the line between the ALP vertex and
the photon-1 vertex (paired with $q_1^2$), and $1-x-y$ the line between the two photons. Only $c_{hh}^A-c_{hh}^P$ and $c_{hh}^P$ appear, as the invariance below requires.
In the form used for the derivation, $F_h=\frac{2N_cQ_h^2}{4\pi^2f_a}[c_{hh}^A-(c_{hh}^A-c_{hh}^P)(1-B_4)]$.

Checks (all executed in the script):

- **Anchor.** For $m_h\to0$, $B_4\to1$ and $F_h\to2N_cQ_h^2c_{hh}^A/(4\pi^2f_a)$, the coefficient of $2N_c\langle Q^2c_q^A\rangle$ in
  $C_{a\gamma\gamma}=c_\gamma+2N_c\langle Q^2(c_q^A+\mathbb M_a^{-1}(\bm{\mathcal C}))\rangle$ of the [LSZ record](../../rcht/odd-parity/odd-alp-source-lsz-assembly.md), with the same sign.
- **Decoupling.** For $m_h\to\infty$, $B_4\to0$ and $F_h\to\frac{2N_cQ_h^2}{4\pi^2f_a}c_{hh}^P$: the derivative coupling decouples and the
  pseudoscalar coupling remains.
- **Rotation invariance.** $c_\gamma/(4\pi^2f_a)+F_h$ is unchanged under $c_{hh}^A\to c_{hh}^A+\kappa$, $c_{hh}^P\to c_{hh}^P+\kappa$,
  $c_\gamma\to c_\gamma-2N_cQ_h^2\kappa$ for arbitrary $m_h,q_1^2,q_2^2,m_a^2$. The shift of $c_\gamma$ is the one that makes the sum invariant; it is not derived from the path-integral measure.
- **Real photons.** For $q_1^2=q_2^2=0$, $B_4=B_1(\tau_a)=1-\tau_af^2(\tau_a)$, with $f=\arcsin(1/\sqrt{\tau_a})$ for $m_a<2m_h$. The bracket of $F_h$ is then
  $c_{hh}^P+(c_{hh}^A-c_{hh}^P)B_1(\tau_a)$, the structure of the $c_{gg}$ formula $c_{GP}+\sum c_{qq}^{AP}B_1$ of [@Alda:2026nrd] with the photon replaced by a gluon, and the structure
  $c_{\gamma\gamma}+\sum N_cQ_f^2c_{ff}B_3$ of [@Bauer:2021mvw; @Alda:2024cxn] for one real photon, to which $B_4(\tau_a,\tau_1,\infty)=B_3(\tau_a,\tau_1)$ reduces.
- **Vector coupling.** With the vertex $(c_{hh}^V/f_a)\slashed p$ and both photon orderings, the numerator traces have no Levi-Civita structure (their projection on
  $\epsilon^{\mu\nu\rho\sigma}q_{1\rho}q_{2\sigma}$ vanishes identically, whereas the same projection of the pseudoscalar numerator gives $4im_h$), and in addition the two orderings
  cancel, $N_1^{\mu\nu}(k)+N_2^{\mu\nu}(-k-p)=0$, by charge conjugation (Furry's theorem). Hence $F_h$ does not depend on $c_{hh}^V$. The cancellation of the orderings uses a shift
  of the loop momentum in a linearly divergent integral, which is allowed in dimensional regularization or in the Pauli-Villars difference; the vanishing of the coefficient of $R$ does not need it.

# Limits and continuation

The limits at large virtuality, at low energy, for one real photon, and the imaginary part above the pair threshold are properties of
$B_4$ and are collected, with their checks, in the [function record](alp-loop-function-b4.md). In terms of them: at large virtuality $F_h\to\frac{2N_cQ_h^2}{4\pi^2f_a}c_{hh}^A$
with corrections of order $\frac{m_h^2}{|q^2|}\ln\frac{|q^2|}{m_h^2}$ (equal virtualities) or $\frac{m_h^2}{|q^2|}\ln^2\frac{|q^2|}{m_h^2}$ (one real photon); at low energy
$B_4\to-\frac{m_a^2+q_1^2+q_2^2}{12m_h^2}$, so $F_h\to\frac{2N_cQ_h^2}{4\pi^2f_a}c_{hh}^P$ plus corrections suppressed by $1/m_h^2$; and above the pair threshold
$\operatorname{Im}F_h=\frac{2N_cQ_h^2}{4\pi^2f_a}(c_{hh}^A-c_{hh}^P)\operatorname{Im}B_4$ with $\operatorname{Im}B_4\le0$ (for $m_h^2\to m_h^2-i0$).

# Relation to the local operator expansion

For $|q^2|\gg m_h^2$ the triangle splits into a constant, $\frac{2N_cQ_h^2}{4\pi^2f_a}c_{hh}^A$ (the term without a matrix element of a light operator), and a remainder of
relative order $m_h^2/|q^2|$. The one-loop vacuum-to-ALP axial matrix element of $h$ from the same couplings, inserted in the local expansion with the free-quark coefficient of the
[heavy-flavor record](alp-heavy-flavor-thresholds.md), reproduces at leading power in $1/|q^2|$ and first order in $m_a^2/m_h^2$ the sign, the factor $2m_h^2/|q^2|$, the whole dependence on $\ln m_h^2$, and
the term $m_a^2/(3|q^2|)$ of the expansion of $F_h$. It does not reproduce a remainder that is polynomial in $m_h^2$, independent of $m_a$, and comes from loop momenta of order $|q^2|$. Hence for the heavy
flavor the local expansion with the derivative-free axial operator and the free-quark coefficient alone does not give $F_h$: other operators of relative order $m_h^2/|q^2|$ (for example $m_h\bar hi\gamma_5h$
and derivative operators) contribute. The one-loop matrix element has a logarithmic ultraviolet divergence whose finite part depends on the regularization, and its explicit form is not needed for $F_h$ and
is not recorded here. For $h=c,b$ the contribution to the form factor is the triangle above; the axial matrix-element term is not added for these flavors. For $u,d,s$ the hadronic description of the
[axial-current record](alp-axial-current-matrix-element.md) is used, as set out in the [assembly record](alp-high-energy-tff-assembly.md).

# Open points

1. **Charmonium region.** The free-charm triangle has a $c\bar c$ cut for $m_a>2m_c$ and does not contain resonances; the supported ALP-mass range reaches
   $2m_\tau=3.554$ GeV, which spans $\eta_c$ and $J/\psi$. The treatment of charm in this region is a decision still open.
2. **Mass scheme.** $m_h$ is taken as a single parameter; whether to use $m_h(m_h)$ or a pole mass is not fixed at this order.
3. **Gluonic and QCD corrections.** Terms of relative order $\alpha_s$ and the ALP–gluon coupling $c_G$ in the heavy-quark sector are not included.
4. **Convention checks.** The photon-dual orientation ($\epsilon^{0123}=+1$, canonical in the repository following Bauer:2020jbp) and the shift $c_\gamma\to c_\gamma-2N_cQ_h^2\kappa$
   are fixed by the [photon phase record](../photon/alp-photon-odd-action-phase-matching.md); the shift is now also obtained from the loop, because rotation invariance of $c_\gamma/(4\pi^2f_a)+F_h$
   requires exactly it, given the Pauli-Villars anomaly term. The identification $T^{\mu\nu}\leftrightarrow i\mathcal M=ie^2RF$ follows from the Wick contraction bookkeeping,
   $i\mathcal M=ie^2\varepsilon^*_\mu\varepsilon^*_\nu T^{\mu\nu}$ for outgoing photons with $T^{\mu\nu}=i\int d^4x\,e^{iq\cdot x}\langle0|T\{J^\mu(x/2)J^\nu(-x/2)\}|a\rangle$, and the
   triangle and the local expansion agree in sign in the soft-region comparison above.
5. **Soft-region comparison.** The comparison above was done with a cutoff on the loop momentum and with the odd integral set to zero (as in dimensional regularization); the scheme-dependent constant affects only the part that is not reproduced, not $F_h$.
6. **Notation.** $B_4$ is the approved notation; see the [function record](alp-loop-function-b4.md).

# Verification

Run
```text
.venv/bin/python scripts/derive_alp_heavy_quark_triangle.py
.venv/bin/python scripts/derive_alp_heavy_quark_triangle.py --test
```
Thirty-eight unit tests pass: the trace numerator and the symmetry of the integral, the heavy-quark and massless limits (massless shown for the two spacelike closed forms, not for general
kinematics), the real-photon closed form by differentiation in $m_a^2/m_h^2$, the Pauli-Villars derivative-coupling result and the anchor including its sign, the rotation invariance, the
expansions listed above, and the soft-region comparison.
