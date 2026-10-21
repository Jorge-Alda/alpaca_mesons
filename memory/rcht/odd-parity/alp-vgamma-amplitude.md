---
title: No-mixing ALP amplitudes to omega gamma and phi gamma
status: derived
created: 2026-10-21
updated: 2026-10-21
bibtex_keys: [Alda:2026nrd, Ruiz-Femenia:2003jdx, Ecker:1988te, Ecker:1989yg]
provenance: analytical-derivation
---

# Inputs and scope

This record derives the retained mesonic odd-parity contribution to
$a\to\omega\gamma$ and $a\to\phi\gamma$. The ALP couplings are supplied
already in Alda:2026nrd's no-mixing basis. The odd-parity Lagrangian contains
the original seven
VJP and four VVP operators, with the canonical tensor convention
$\epsilon^{0123}=+1$; the coefficient and tensor signs include the adapter
from the opposite orientation used in the indexed original-operator
contractions [see the
[Levi-Civita convention record](levi-civita-convention.md)]. The even-parity
Lagrangian contributes the fitted electromagnetic source with its stated
$\chi_+$ insertions. Vector pole masses and the mixing angle are physical
inputs; this calculation adds no mass term to determine them. The result is
tree level and leading order in $e$, with narrow on-shell external vectors.
It does not add an ALP--meson pole term or an independent local
anomalous $aV\gamma$ term.

Use the independent diagonal no-mixing couplings $c_q^A$ and $c_q^P$ of
Alda:2026nrd, with $\chi=\operatorname{diag}(2B_0m_u,2B_0m_d,2B_0m_s)$.
The odd vertices contain the derivative coupling as $2c_q^A/f_a$ and the
mass insertion as $-2\{c_q^P,\chi\}/f_a$. No relation between these two
couplings is imposed [@Alda:2026nrd; see the
[two-vector ALP derivation](alp-vv-amplitude.md)].

# Even-sector photon response

Use $Q=\operatorname{diag}(2/3,-1/3,-1/3)$ and the fitted even-parity
electromagnetic source
$$
S_\gamma=F_VQ+z_V\{\chi_+,Q\}
+y_V\operatorname{tr}(\chi_+Q)\mathbf1,
\qquad R_i=\operatorname{tr}(T_iS_\gamma).
$$
The single-trace source fit has $y_V=0$; the double-trace source option keeps
$y_V$. Take $M_\rho,M_\omega,M_\phi$ and $\theta_V$ from the selected physical
vector inputs. With $T_{\omega_I}=\operatorname{diag}(1,1,0)/\sqrt2$ and
$T_{\phi_I}=\operatorname{diag}(0,0,1)$, use
$$
T_\rho=\operatorname{diag}(1,-1,0)/\sqrt2,\quad
T_\omega=\cos\theta_VT_{\omega_I}-\sin\theta_VT_{\phi_I},\quad
T_\phi=\sin\theta_VT_{\omega_I}+\cos\theta_VT_{\phi_I}.
$$
This amplitude calculation adds no vector-mass operator or mass parameter;
the physical pole masses and mixing angle are inputs. Masses are in GeV,
$F_V$ in GeV, $z_V,y_V$ in GeV$^{-1}$, and the odd-operator normalization
is $M_V=0.770$ GeV. The source convention is the one recorded in the
[quadratic photon--vector derivation](../foundations/rcht-photon-vector-diagonalization.md)
[@Ecker:1988te; @Ecker:1989yg].

At leading electromagnetic order the diagonal photon--vector equations give
the zero-virtuality response
$$
\mathcal W_\gamma(0)=\sum_{i=\rho,\omega,\phi}
  \frac{R_i}{M_i^2}T_i,
\qquad M_i^2r_i=\sqrt2eR_i,
\qquad V_{\mu\nu}=\sqrt2e\,\mathcal W_\gamma(0)F_{\mu\nu}.
$$
This retains the three physical pole masses and their fitted source residues;
it is not a common-pole approximation.

# Odd-parity amplitudes

For either channel, let $T_V$ be $T_\omega$ or $T_\phi$ and let
$M_{\rm ext}$ be the mass of that outgoing vector:
$$
M_{\rm ext}=M_\omega\quad\text{for }a\to\omega\gamma,
\qquad M_{\rm ext}=M_\phi\quad\text{for }a\to\phi\gamma.
$$
$M_V=0.770$ GeV remains the odd-operator normalization scale; it is not the
outgoing-vector mass. The amplitude coefficient in
$$
\mathcal M(a\to V\gamma)=G_{V\gamma}\,
\epsilon_{\mu\nu\rho\sigma}\varepsilon_V^\mu
\varepsilon_\gamma^\nu q_V^\rho q_\gamma^\sigma
$$
is
$$
\begin{aligned}
G_{V\gamma}=\frac{e}{M_{\rm ext}}\Bigg\{&
-\frac{8}{M_Vf_a}\Bigg[
\operatorname{tr}(T_VQc_q^A)\Big(
c_1(m_a^2-M_{\rm ext}^2)+c_2(m_a^2+M_{\rm ext}^2)
-c_5(m_a^2-M_{\rm ext}^2)-2c_6M_{\rm ext}^2\Big)\\
&\hspace{27mm}+4c_3\operatorname{tr}
  \left(T_VQ\{c_q^P,\chi\}\right)\Bigg]\\
&+\frac{8\sqrt2}{f_a}\Bigg[
\operatorname{tr}(T_V\mathcal W_\gamma(0)c_q^A)
  \Big((d_1-d_3)m_a^2+d_3M_{\rm ext}^2\Big)\\
&\hspace{27mm}+4d_2\operatorname{tr}
  \left(T_V\mathcal W_\gamma(0)\{c_q^P,\chi\}\right)\Bigg]\Bigg\}.
\end{aligned}
$$
The first bracket is the direct VJP contribution. The second is VVP with one
vector leg converted to the real photon by the diagonalized response; its
formula includes both VVP vector assignments. The factor $1/M_{\rm ext}$ is
the external-vector tensor normalization. The formula applies for
$m_a>M_{\rm ext}$; no width or branching fraction is included.

# Fit-null projection

Use coefficient order $(c_1,\ldots,c_7,d_1,\ldots,d_4)$. The stored primary
fit nullspace has four directions, represented exactly by shifts in $c_7$,
$d_4$, and the two VJP combinations
$$
\delta c_5=\delta c_2+\delta c_1,\qquad
\delta c_6=\delta c_2,
$$
with all other shifts zero. A basis for the latter plane is a shift with only
$\delta c_1=\delta c_5$ and a shift with
$\delta c_2=\delta c_5=\delta c_6$. Differentiating both physical amplitude
rows with respect to the eleven coefficients and multiplying by these four
columns gives
$$
\begin{pmatrix}
\partial G_{\omega\gamma}/\partial(c_1,\ldots,d_4)\\
\partial G_{\phi\gamma}/\partial(c_1,\ldots,d_4)
\end{pmatrix}N=0_{2\times4}.
$$
The executable also projects the floating null columns stored for every
primary source hypothesis onto this exact subspace; the largest Frobenius
distance is $9.84\times10^{-15}$. Thus both amplitudes are insensitive to the
fit null directions at fixed even-sector inputs. In addition, $c_4$ is absent
because the electromagnetic-compatible neutral source has $f_-=0$; it is not
one of these four fit-null directions. The $c_7$ on-shell wedge contraction
and neutral $d_4$ contraction vanish separately. This result does not assign
errors to identifiable combinations or cover nuisance refits, additional odd
operators, loops, or unstable-vector line shapes.

# Reproduction

Run the executable SymPy derivation:

```text
.venv/bin/python scripts/derive_physical_alp_vgamma.py
```

It derives the Lorentz contractions, ALP source matrices, even source,
diagonal pole response, physical amplitude gradients, and fit-null
projections. Its raw Ward-identity, response-stationarity and
universal/common-mass identities are retained as assertions in the script.
