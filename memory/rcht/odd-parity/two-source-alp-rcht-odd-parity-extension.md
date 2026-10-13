---
title: Two-source ALP anomaly seed and odd-parity RChT extension
status: derived
created: 2026-09-28
updated: 2026-09-30
bibtex_keys: [Alda:2026nrd, Bai:2025fvl, Ruiz-Femenia:2003jdx, Ecker:1988te, Ecker:1989yg]
provenance: analytical-derivation
---

# Two-Source ALP Anomaly Seed and Odd-Parity RChT Extension

## Scope, evidence, and obstruction

This record extends the [no-vector photon invariant](../../alp/photon/alp-photon-invariant-from-mixing.md) to two vector-like sources and interfaces it with the [repository RChT vectors](../foundations/rcht-photon-vector-diagonalization.md). The calculation is performed by [`scripts/derive_two_source_rcht_anomaly.py`](../../../scripts/derive_two_source_rcht_anomaly.py). Inputs are Alda's Appendix-A CCWZ and source formulas in [`appendix-a-expansion.md`](../../extracts/Alda:2026nrd/appendix-a-expansion.md) [@Alda:2026nrd], the vector-only background-field limitations in [`no-axial-vector-and-rcht-matching.md`](../../extracts/Bai:2025fvl/no-axial-vector-and-rcht-matching.md) [@Bai:2025fvl], and the independent odd-parity operator basis in [`vector-odd-parity.md`](../../extracts/Ruiz-Femenia:2003jdx/vector-odd-parity.md) [@Ruiz-Femenia:2003jdx].

> **WZW anomaly matching does not uniquely predict the complete $aVV$ vertex in antisymmetric-tensor RChT.** Independent VVP and VJP/homogeneous coefficients $d_i,c_i$, together with singlet/double-trace and $SU(3)$-breaking terms, require short-distance constraints and/or data. This obstructs a unique numerical result, but not the general extension below.

Literature evidence consists of the cited CCWZ/source definitions and RChT operator bases. Every boxed extension and change-of-basis formula below is a derived result. No value for an independent odd-parity resonance coefficient is inferred from the anomaly.

## Conventions, units, and assumptions

Natural units are used. The source couplings $g_A,g_B$, generators, and all $c$ coefficients are dimensionless. $a,f_a,f_\pi,B_0,m_q,M_i$ and derivatives have mass dimension one. The Hermitian source generators are $T_A,T_B$. Flavor traces use $\operatorname{tr}(\lambda_i\lambda_j^\dagger)=2\delta_{ij}$.

Use one unordered-pair, symmetric coefficient convention,

$$
\mathcal L_{a\mathcal F\mathcal F}
=\frac{1}{16\pi^2}\frac{a}{f_a}
\sum_{A\leq B}g_Ag_B C_{AB}\,
F^A_{\mu\nu}\widetilde F^{B\mu\nu},
\qquad C_{AB}=C_{BA}.
$$

Thus every unordered pair, diagonal or unlike, occurs once. For the photon, $g_\gamma=e$, so the diagonal normalization is $e^2/(16\pi^2)=\alpha_{\rm em}/(4\pi)$. This convention is used consistently and is not silently changed by a symmetry factor.

The derivation retains one ALP at $O(1/f_a)$, tree-level fields, and a complete $U(3)$ pseudoscalar basis. It assumes $\mathbb M_a^{-1}$ exists, so it excludes exact ALP--meson degeneracies with nonzero overlap. The quadratic photon-vector transformation is exact only in the already specified truncated tree-level RChT quadratic model. No ideal vector mixing, on-shell reduction, short-distance relation, flavor commutativity, or numerical approximation is imposed.

Alda Eq. (III.37) states, for every physical pseudoscalar $i$,

$$
|m_a^2-m_i^2|\gg\frac{f_\pi}{f_a}
\left|\operatorname{tr}(C\lambda_i^\dagger)\right|
\approx\frac{f_\pi}{f_a}m_i^2.
$$

Thus the excluded real-mass neighborhoods scale as $f_\pi/f_a$ and are generally negligible in measure or scan resolution when $f_a\gg f_\pi$, but the inequality is a local, pointwise validity condition. A finite caller-supplied $\mathbb M_a^{-1}(C)$ is necessary but not sufficient: each scan point must be tested with its exact overlap, and a violation must be masked/flagged or treated with coupled mixing. Matrix magnitude alone is not the test. For stable or sufficiently narrow pseudoscalars, a retained violating point requires coupled diagonalization. For unstable states, real-mass proximity and complex-pole proximity are distinct, and the relevant replacement is a coupled complex-propagator/pole analysis [@Alda:2026nrd].

## 1. Symmetric two-source anomaly seed

Because the Lagrangian contains each unordered pair only once, define the
source multiplicity

$$
s_{AB}\equiv\frac{2}{1+\delta_{AB}}.
$$

It is $s_{AA}=1$ for a diagonal pair and $s_{AB}=2$ for $A\neq B$, where it
combines the two ordered cross terms in the background-field anomaly.
The symmetric two-source extension of the no-vector WZW seed is

$$
\boxed{
C_{AB}^{\rm seed}=c_{AB}
+s_{AB}N_c\operatorname{tr}\!\left(\{T_A,T_B\}[c_q^A+\mathbb M_a^{-1}(C)]\right).}
$$

Under Alda's axial reparametrization,

$$
c_q^A\to c_q^A+\kappa,
\qquad
c_{AB}\to c_{AB}-s_{AB}N_c
\operatorname{tr}(\{T_A,T_B\}\kappa),
$$

the script expands the transformed seed and obtains

$$
C_{AB}^{\rm seed\prime}-C_{AB}^{\rm seed}=0.
$$

For $A=B=\gamma$, $T_A=T_B=Q$, $g_A=g_B=e$, and $c_{AA}=c_\gamma$, hence

$$
\boxed{C_{\gamma\gamma}^{\rm seed}
=c_\gamma+2N_c\operatorname{tr}\!\left(Q^2[c_q^A+\mathbb M_a^{-1}(C)]\right),}
$$

exactly reproducing the [repository photon result](../../alp/photon/alp-photon-invariant-from-mixing.md), including its sign and normalization.

For distinct sources, $s_{AB}=2$, so the corresponding unordered coefficient is

$$
C_{AB}^{\rm seed}=c_{AB}
+2N_c\operatorname{tr}\!\left(\{T_A,T_B\}[c_q^A+\mathbb M_a^{-1}(C)]\right),
\qquad A\neq B.
$$

For hadronic RChT vectors, $c_{AB}$ is **not** a UV Wilson coefficient fixed by Alda. It parametrizes the local background-field/anomaly convention used when vector-like sources are matched onto resonances. Independent homogeneous resonance operators remain in addition to this seed.

## 2. Physical-ALP insertions by completeness

Alda's physical-ALP meson components are

$$
w_i^{(a)}=\frac{f_\pi}{f_a}
\operatorname{tr}[(c_q^A+\mathbb M_a^{-1}(C))\lambda_i^\dagger].
$$

Completeness, rather than an assumed matrix replacement, gives

$$
\sum_iw_i^{(a)}\lambda_i
=\frac{2f_\pi}{f_a}(c_q^A+\mathbb M_a^{-1}(C)).
$$

At the Goldstone vacuum, Alda Appendix A gives a direct term $2(\partial_\mu a/f_a)c_q^A$ and a meson term $-(\partial_\mu\Phi)/f_\pi$. Their LSZ sum is therefore

$$
\boxed{
u_\mu\big|_{a_{\rm physical}}
=-2\frac{\partial_\mu a}{f_a}\mathbb M_a^{-1}(C).}
$$

For the pseudoscalar source, Alda's convention is $c_q^{AP}=c_q^A-c_q^P$. From $\chi=\chi_S+i\chi_P$ and $\chi_-=u^\dagger\chi u^\dagger-u\chi^\dagger u$, the direct vacuum term is

$$
\chi_-\big|_{a,\rm direct}
=4iB_0\frac{a}{f_a}\{c_q^P,m_q\}.
$$

Expanding $u=\exp[i\Phi/(2f_\pi)]$ fixes the meson sign,

$$
\chi_-\big|_{\Phi}
=-\frac{2iB_0}{f_\pi}\{\Phi,m_q\}.
$$

The same completeness sum then yields

$$
\boxed{
\chi_-\big|_{a_{\rm physical}}
=-4iB_0\frac{a}{f_a}
\{c_q^{AP}+\mathbb M_a^{-1}(C),m_q\}.}
$$

This agrees with the no-mixing variable $\widehat c_q^P=-c_q^{AP}-\mathbb M_a^{-1}(C)$ and fixes the sign independently through the direct-plus-meson expansion.

### Insertions in the four VVP operators

Ruiz--Femenia's operators are retained as distinct Lorentz structures [@Ruiz-Femenia:2003jdx]:

$$
\begin{aligned}
\mathcal O^1_{VVP}&=\epsilon_{\mu\nu\rho\sigma}
\langle\{V^{\mu\nu},V^{\rho\alpha}\}\nabla_\alpha u^\sigma\rangle,\\
\mathcal O^2_{VVP}&=i\epsilon_{\mu\nu\rho\sigma}
\langle\{V^{\mu\nu},V^{\rho\sigma}\}\chi_-\rangle,\\
\mathcal O^3_{VVP}&=\epsilon_{\mu\nu\rho\sigma}
\langle\{\nabla_\alpha V^{\mu\nu},V^{\rho\alpha}\}u^\sigma\rangle,\\
\mathcal O^4_{VVP}&=\epsilon_{\mu\nu\rho\sigma}
\langle\{\nabla^\sigma V^{\mu\nu},V^{\rho\alpha}\}u_\alpha\rangle.
\end{aligned}
$$

For $V=T_AV_A+T_BV_B$ and no additional background field in $\nabla$, define the normalized labeled symmetrization

$$
\mathcal S_{AB}[X_AV_B]
\equiv\frac{X_AV_B+X_BV_A}{1+\delta_{AB}}.
$$

It retains both assignments for distinct vectors and avoids double counting for
identical vectors. The physical-ALP pieces are

$$
\begin{aligned}
\mathcal O^1_{ABa}&=-\frac{2}{f_a}\epsilon_{\mu\nu\rho\sigma}
(\partial_\alpha\partial^\sigma a)
\left\langle\{T_A,T_B\}\mathbb M_a^{-1}(C)\right\rangle
\mathcal S_{AB}\!\left[V_A^{\mu\nu}V_B^{\rho\alpha}\right],\\
\mathcal O^2_{ABa}&=+\frac{4B_0}{f_a}\epsilon_{\mu\nu\rho\sigma}a
\left\langle\{T_A,T_B\}\{c_q^{AP}+\mathbb M_a^{-1}(C),m_q\}\right\rangle
\mathcal S_{AB}\!\left[V_A^{\mu\nu}V_B^{\rho\sigma}\right],\\
\mathcal O^3_{ABa}&=-\frac{2}{f_a}\epsilon_{\mu\nu\rho\sigma}
(\partial^\sigma a)
\left\langle\{T_A,T_B\}\mathbb M_a^{-1}(C)\right\rangle
\mathcal S_{AB}\!\left[(\partial_\alpha V_A^{\mu\nu})V_B^{\rho\alpha}\right],\\
\mathcal O^4_{ABa}&=-\frac{2}{f_a}\epsilon_{\mu\nu\rho\sigma}
(\partial_\alpha a)
\left\langle\{T_A,T_B\}\mathbb M_a^{-1}(C)\right\rangle
\mathcal S_{AB}\!\left[(\partial^\sigma V_A^{\mu\nu})V_B^{\rho\alpha}\right].
\end{aligned}
$$

The normalized $\mathcal S_{AB}$ exchanges the complete labeled vector factors,
preserving which factor carries each Lorentz index and derivative. Consequently
**all four $d_i$ operators can contribute** when the full Alda pseudoscalar
source is retained; in particular, $d_2$ must not be omitted. Their distinct
derivative structures have not been combined, integrated by parts, or reduced
with on-shell equations.

## 3. Exact quadratic change of basis for anomalous bilinears

Use field-strength variables

$$
X_{\rm old}=(F_{\rm old},\mathbf U)^T,
\qquad U_i=M_iV_i,
\qquad
X_{\rm phys}=(F_\gamma,\mathbf W)^T.
$$

With the definitions in the [quadratic diagonalization](../foundations/rcht-photon-vector-diagonalization.md) and [interaction expansion](../foundations/rcht-vpp-vpa-photon-pa-expansion.md),

$$
H=K_V^{1/2}O,
\qquad
\boldsymbol\rho=\frac{O^T\boldsymbol\epsilon}{\sqrt Z},
$$

the existing exact map is

$$
\boxed{
X_{\rm old}=S X_{\rm phys},\qquad
S=\begin{pmatrix}
Z^{-1/2}&-\boldsymbol\rho^T\\
Z^{-1/2}\boldsymbol\epsilon&H
\end{pmatrix}.}
$$

For a general symmetric anomalous bilinear matrix

$$
\Lambda_{\rm old}=
\begin{pmatrix}\lambda_{\gamma\gamma}&\boldsymbol\ell^T\\
\boldsymbol\ell&\Lambda_V\end{pmatrix},
$$

the matrix convention is defined independently of the unordered coefficients
in Sec. 1 by

$$
\mathcal L_\Lambda
=\mathcal N\frac{a}{f_a}\frac12
X_{\rm old}^T\Lambda_{\rm old}\widetilde X_{\rm old},
$$

where the dual acts componentwise and $\mathcal N$ is a common dimensionless
normalization. Since every component of $X_{\rm old}$ has mass dimension two,
$\Lambda_{\rm old}$ is dimensionless. Symmetry of $\Lambda_{\rm old}$ makes an
off-diagonal matrix entry occur twice inside the one-half bilinear, whereas a
diagonal entry occurs once with the explicit factor $1/2$. Thus, if this matrix
is matched to an unordered coefficient matrix in the same field normalization,
$\Lambda_{AB}=C_{AB}$ for $A\neq B$ but
$\Lambda_{AA}=2C_{AA}$; couplings and any tensor-to-$U_i=M_iV_i$ matching
factors must also be included in $\Lambda$. The symbols $\Lambda_{AB}$ and
$C_{AB}$ must therefore not be identified entry by entry without this conversion.

direct symbolic congruence gives

$$
\boxed{\Lambda_{\rm phys}=S^T\Lambda_{\rm old}S}
$$

with blocks

$$
\boxed{
\lambda_{\gamma\gamma}^{\rm phys}
=\frac{\lambda_{\gamma\gamma}+2\boldsymbol\ell^T\boldsymbol\epsilon
+\boldsymbol\epsilon^T\Lambda_V\boldsymbol\epsilon}{Z},}
$$

$$
\boxed{
\boldsymbol\lambda_{\gamma W}^{T\,\rm phys}
=\frac{1}{\sqrt Z}\left[
(\boldsymbol\ell^T+\boldsymbol\epsilon^T\Lambda_V)H
-(\lambda_{\gamma\gamma}+\boldsymbol\epsilon^T\boldsymbol\ell)
\boldsymbol\rho^T\right],}
$$

$$
\boxed{
\Lambda_{WW}^{\rm phys}
=H^T\Lambda_VH
-\boldsymbol\rho\boldsymbol\ell^TH
-H^T\boldsymbol\ell\boldsymbol\rho^T
+\lambda_{\gamma\gamma}\boldsymbol\rho\boldsymbol\rho^T.}
$$

This constant-matrix congruence applies to field-strength bilinears and the anomaly seed. The derivative $d_1,d_3,d_4$ structures must each be transformed at operator level; because derivatives occupy inequivalent slots, they cannot all be represented by one constant $\Lambda$.

## Verification

Run

```text
.venv/bin/python scripts/derive_two_source_rcht_anomaly.py
```

The script constructs a complete normalized $U(3)$ basis, derives the LSZ sums, checks both CCWZ insertions including the $\chi_-$ sign, proves seed invariance and exact photon reduction, and derives the three blocks of $S^T\Lambda S$. All assertions must print as passed.
