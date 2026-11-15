---
title: Vacuum-to-ALP quark axial-current matrix elements at LO ChPT
status: derived
created: 2026-11-15
updated: 2026-11-15
bibtex_keys: [Alda:2026nrd, Hoferichter:2020lap, Bauer:2020jbp]
provenance: analytical-derivation
---

# Scope, inputs, and assumptions

This record gives the vacuum-to-ALP matrix element of the quark axial current weighted
by a diagonal flavor matrix $X$,
$\langle0|\bar q\gamma^\mu\gamma_5Xq|a(p)\rangle$, at leading order (LO) in chiral
perturbation theory (ChPT), in the no-mixing source basis of Alda:2026nrd. The case
needed for the short-distance expansion is $X=Q^2$, where $Q=\operatorname{diag}(\tfrac23,-\tfrac13,-\tfrac13)$
is the quark-charge matrix. It supplies the matrix element left open in the
[leading hard-current expansion](alp-hard-current-ope.md). Conventions are those of the
[shared-conventions record](alp-high-energy-tff-conventions.md). In this record $q^2$
is the photon virtuality (spacelike for $q^2<0$) and $Q$ always denotes the charge matrix.
All flavor matrices are diagonal in the $(u,d,s)$ basis, as in the
[tree TFF](alp-doubly-virtual-tff.md), and are written in matrix notation.

Inputs and their provenance:

- ALP coupling $\frac{\partial_\mu a}{f_a}\bar q\gamma^\mu(c_q^V+c_q^A\gamma_5)q$ and the
  vielbein $u_\mu=2\frac{\partial_\mu a}{f_a}c_q^A-\frac1{f_\pi}\partial_\mu\Phi+\dots$
  [@Alda:2026nrd, `main_arXiv.tex` lines 228--240 and 452], with $f_\pi=93$ MeV
  and $\Phi=\sum\lambda_PP$ ([neutral generators](../../extracts/Alda:2026nrd/neutral-meson-generators.md)).
- Leading chiral Lagrangian $\frac{f_\pi^2}4\langle u_\mu u^\mu+\chi_+\rangle$ plus the
  topological term [@Alda:2026nrd, lines 550--557; [extract](../../extracts/Alda:2026nrd/lo-topological-action.md)].
- The no-mixing parametrization of [@Alda:2026nrd, Eqs. (III.30), (III.33), Sec. III E], as
  assembled in the [LSZ record](../../rcht/odd-parity/odd-alp-source-lsz-assembly.md): the physical ALP
  has no mesonic component at order $1/f_a$, and the coefficient of $\partial_\mu a$ in $u_\mu$ is
  $$
  u_\mu\big|_a=-\frac{2}{f_a}\,\mathbb M_a^{-1}(\bm{\mathcal C})\,\partial_\mu a .
  $$
  Here $\langle\cdot\rangle$ is the flavor trace,
  $$
  \mathbb M_a(\bm A)=B_0\{m_q,\bm A\}+\frac{m_0^2}{N_f}\langle\bm A\rangle\mathbf 1-m_a^2\bm A,\qquad
  \bm{\mathcal C}=-B_0\{m_q,c_q^{AP}\}-\frac{m_0^2}{N_f}c_{GA}\mathbf 1,
  $$
  with $c_q^{AP}=c_q^A-c_q^P$ and $c_{GA}=c_G+\langle c_q^A\rangle$ [@Alda:2026nrd, Sec. III E].
  Both $\bm{\mathcal C}$ and $\mathbb M_a$ are invariant under ALP-dependent quark-field redefinitions, so the
  coefficient of $\partial_\mu a$ above is too.

Assumptions (not results): one ALP at tree level and order $1/f_a$; diagonal sources;
the no-mixing parametrization, in which the meson overlap is already contained in
$\mathbb M_a^{-1}(\bm{\mathcal C})$, so that $\langle0|\Phi|a\rangle=0$ at this order and no second meson
pole is added. The constant $f_\pi$ is the LO chiral coupling of the effective Lagrangian, not a meson
decay constant used as an ALP input; no $\eta$ or $\eta'$ quantity enters.

# Derivation

The paired [script](../../../scripts/derive_alp_axial_current_matrix_element.py)
adds a diagonal external axial source $A_\mu$ with
$\mathcal L_{\rm QCD}\supset\bar q\gamma^\mu\gamma_5A_\mu q$, so that
$\bar q\gamma^\mu\gamma_5Xq=\partial\mathcal L/\partial A_\mu$ contracted with $X$ in flavor
space. Because the ALP coupling has the same form with $A_\mu=c_q^A\partial_\mu a/f_a$, the
source enters the vielbein as $u_\mu\supset2A_\mu$. The terms $\chi_+$ and the topological
term do not depend on $A_\mu$ at this order, so only $\frac{f_\pi^2}4\langle u_\mu u^\mu\rangle$
contributes, and
$$
\bar q\gamma^\mu\gamma_5Xq\;\to\;f_\pi^2\bigl\langle X\,u^\mu\bigr\rangle
=-\frac{2f_\pi^2}{f_a}\bigl\langle X\,\mathbb M_a^{-1}(\bm{\mathcal C})\bigr\rangle\,\partial^\mu a
-f_\pi\bigl\langle X\,\partial^\mu\Phi\bigr\rangle.
$$
With $\langle0|a(x)|a(p)\rangle=e^{-ip\cdot x}$ the term with $\Phi$ drops out
(no-mixing parametrization) and
$$
\boxed{\langle0|\bar q\gamma^\mu\gamma_5Xq|a(p)\rangle=\frac{2if_\pi^2}{f_a}\,p^\mu\bigl\langle X\,\mathbb M_a^{-1}(\bm{\mathcal C})\bigr\rangle,}
\qquad
p_\mu\langle0|\bar q\gamma^\mu\gamma_5Xq|a(p)\rangle=\frac{2if_\pi^2}{f_a}\,m_a^2\bigl\langle X\,\mathbb M_a^{-1}(\bm{\mathcal C})\bigr\rangle.
$$

## Isovector, octet, and singlet parts

For $X=Q^2$, write $Q^2=c_3\lambda_3+c_8\lambda_8+c_0\lambda_0$ with $\lambda_0=\sqrt{2/3}\,\mathbf 1$,
$\langle\lambda_i\lambda_j\rangle=2\delta_{ij}$, so that
$c_3=\tfrac16$, $c_8=\tfrac1{6\sqrt3}$, $c_0=\tfrac2{3\sqrt6}$, and the diagonal entries of $Q^2$ are
$\tfrac49,\tfrac19,\tfrac19$. Then
$$
\bigl\langle Q^2\,\mathbb M_a^{-1}(\bm{\mathcal C})\bigr\rangle=
\underbrace{\tfrac16\bigl\langle\lambda_3\,\mathbb M_a^{-1}(\bm{\mathcal C})\bigr\rangle}_{\text{isovector}}
+\underbrace{\tfrac1{6\sqrt3}\bigl\langle\lambda_8\,\mathbb M_a^{-1}(\bm{\mathcal C})\bigr\rangle}_{\text{octet}}
+\underbrace{\tfrac29\bigl\langle\mathbb M_a^{-1}(\bm{\mathcal C})\bigr\rangle}_{\text{singlet}},
$$
and each part multiplies $\frac{2if_\pi^2}{f_a}p^\mu$ in the boxed matrix element. The
charge weights are those of the hard-current record.

## Form with the mass and topological terms displayed

The definition of $\mathbb M_a$ applied to $\mathbb M_a^{-1}(\bm{\mathcal C})$ gives
$m_a^2\,\mathbb M_a^{-1}(\bm{\mathcal C})=B_0\{m_q,\mathbb M_a^{-1}(\bm{\mathcal C})\}
+\frac{m_0^2}{N_f}\langle\mathbb M_a^{-1}(\bm{\mathcal C})\rangle\mathbf 1-\bm{\mathcal C}$. Inserting $\bm{\mathcal C}$,
$$
p_\mu\langle0|\bar q\gamma^\mu\gamma_5Xq|a(p)\rangle=\frac{2if_\pi^2}{f_a}
\Bigl\langle X\Bigl(B_0\bigl\{m_q,\;c_q^{AP}+\mathbb M_a^{-1}(\bm{\mathcal C})\bigr\}
+\frac{m_0^2}{N_f}\bigl(c_{GA}+\langle\mathbb M_a^{-1}(\bm{\mathcal C})\rangle\bigr)\mathbf 1\Bigr)\Bigr\rangle .
$$
The first term comes from the quark masses, the second from the topological term. Dividing
by $m_a^2$ gives the matrix element itself and makes the factor $m_a^{-2}$ explicit; the form
with $\mathbb M_a^{-1}(\bm{\mathcal C})$ above is finite for $m_a\to0$, and the explicit $m_a^{-2}$
cancels once $\mathbb M_a^{-1}(\bm{\mathcal C})$ is evaluated.

# Contribution of the axial-current operator to the hard coefficient

Inserting the boxed matrix element with $X=Q^2$ in the local expansion of the
[hard-current record](alp-hard-current-ope.md), for nearly symmetric hard photons,
$$
\boxed{F_{a\gamma^*\gamma^*}\;\supset\;\frac{2i}{q^2m_a^2}\,p_\mu\langle0|\bar q\gamma^\mu\gamma_5Q^2q|a(p)\rangle
=-\frac{4f_\pi^2}{f_a\,q^2}\bigl\langle Q^2\,\mathbb M_a^{-1}(\bm{\mathcal C})\bigr\rangle,}
$$
with no $m_a^{-2}$ singularity. The trace depends on the sources but not on $q^2$, so this term
falls as $1/q^2$ for $|q^2|\to\infty$.

**This is only one contribution to the large-virtuality TFF, not the complete result.**
It is the matrix element of the operator that multiplies the leading Wilson coefficient.
It depends on the sources only through $c_q^{AP}$ and $c_{GA}$ (via $\bm{\mathcal C}$), so it is
invariant under the ALP-dependent quark-field redefinitions by itself. The explicit ALP couplings
also enter the hard region directly, and these terms are not contained in the matrix element above:

- **Derivative axial coupling $c_q^A$ and the photon coefficient $c_\gamma$.** A hard quark
  loop with the derivative coupling $\partial_\mu a\,\bar q\gamma^\mu\gamma_5q$ is the
  anomalous triangle. For massless quarks and a soft ALP its divergence is the anomaly, so
  it gives a term independent of $q^2$, which combines with $c_\gamma$ into
  $c_\gamma+2N_c\langle Q^2c_q^A\rangle$, invariant under the redefinitions. This is the combination in
  $C_{a\gamma\gamma}=c_\gamma+2N_c\langle Q^2(c_q^A+\mathbb M_a^{-1}(\bm{\mathcal C}))\rangle$ of the
  [LSZ record](../../rcht/odd-parity/odd-alp-source-lsz-assembly.md), with the $\mathbb M_a^{-1}(\bm{\mathcal C})$ part
  carried by the matrix element above. Such a constant would dominate over the $1/q^2$ term at
  large virtuality, and it is consistent with the light-fermion loop structure in
  Bauer:2020jbp. Heavy quarks decouple from it. **Not derived here**; it is part of
  the assembly in task 3.1 and must be derived there, including its dependence on the quark
  masses and on the ALP mass.
- **Pseudoscalar coupling $c_q^P$.** The hard-loop insertion of
  $\frac{a}{f_a}\bar q\,i\gamma_5\{c_q^P,m_q\}q$ is suppressed by powers of $m_q/\sqrt{|q^2|}$ for light quarks,
  and contributes through the $m_a^2/q^2$ and heavy-quark terms (tasks 2.3 and 3.3). Its low-energy
  effect is already in $\bm{\mathcal C}$ through $c_q^{AP}$.
- **Gluonic coupling $c_G$.** The operator $\alpha_sG\widetilde G$ appears in the expansion
  with a coefficient of relative order $\alpha_s$, and through the anomaly its matrix element is
  related to the singlet axial current. Its low-energy effect is in $\bm{\mathcal C}$ through
  $c_{GA}$ (the combination $c_{GP}=c_G+\langle c_q^P\rangle$ enters the gluon coupling of the ALP in
  [@Alda:2026nrd]). The treatment of this operator at order $\alpha_s$ is task 2.2.

The statements in this list are expectations from the anomalous Ward identity and from
Bauer:2020jbp, not results of this record. The matrix element is the unrenormalized LO
ChPT value; scale, scheme, singlet-current treatment, $\alpha_s$ corrections, and heavy flavors
are tasks 2.2 and 2.3.

# Verification

Run
```text
.venv/bin/python scripts/derive_alp_axial_current_matrix_element.py
.venv/bin/python scripts/derive_alp_axial_current_matrix_element.py --test
```
Six unit tests pass: the plane-wave derivative sign, the matrix element, the
flavor decomposition and charge weights, the definition of $\mathbb M_a$ and the mass-plus-topological form
on the solved three-flavor system, the hard coefficient with a check that no
$\eta$-type symbol appears, and the neutral-pion method check. In that last check
the same step with $\Phi=\pi^0\lambda_3$ gives
$\langle0|\bar q\gamma^\mu\gamma_5Xq|\pi^0(p)\rangle=if_\pi p^\mu\langle X\lambda_3\rangle$ and hence, for $X=Q^2$,
$F=-2f_\pi/(3q^2)$, which agrees with the pion limit in Hoferichter:2020lap
(see the [Hoferichter:2020lap extracts](../../extracts/Hoferichter:2020lap/INDEX.md)).
It checks the method only and is not an input for the ALP.

The script prints its results in flavor components, with its own symbol names; these are
the diagonal entries of the matrix expressions above. The relative sign between
$\mathbb M_a^{-1}(\bm{\mathcal C})$ and $\partial_\mu\Phi/f_\pi$ is fixed by the source Lagrangian
and vielbein cited above; the script does not rederive $u_\mu|_a=-2\mathbb M_a^{-1}(\bm{\mathcal C})\,\partial_\mu a/f_a$
or the statement $\langle0|\Phi|a\rangle=0$.
