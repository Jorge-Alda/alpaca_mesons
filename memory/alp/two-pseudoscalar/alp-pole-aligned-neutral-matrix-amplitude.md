---
title: Pole-aligned finite-mixing neutral two-pseudoscalar amplitude
status: derived
created: 2026-11-12
updated: 2026-11-12
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# Physical two-pseudoscalar amplitude at finite mixing

## Inputs and assumptions

Use the LO CCWZ Lagrangian of the [joint expansion](../shared/alp-lo-joint-vertex-expansion.md), the [pole-alignment map](../shared/alp-pole-alignment-chiral-generator.md), and the [meson operator inverse](../shared/alp-lo-finite-mixing-operator-inverse.md), with the conventions of `Alda:2026nrd` [@Alda:2026nrd]. Natural units, tree-level LO, a stationary vacuum with $M(0)=m_q$ Hermitian and diagonal and $P_0=0$ are assumed. The sources commute with the quark charge matrix; the off-diagonal flavour-changing vector source is restricted to the Hermitian $d$–$s$ pair. General Taylor derivatives of the sources are retained. Assume $\widehat\Theta(0)=0$, a purely quadratic independent ALP potential (so its third derivative vanishes), and a cancelled ALP tadpole. No weak mixing, widths, loops or resonances are included.

The three squared masses $s_a,s_i,s_j$ are physical pole positions of the full quadratic system, including the second derivative of the mass source. They are not identified with $B_0(m_q+m_{q'})$. A flavour $K^0$ is not generally a mass eigenstate; the two daughter labels identify physical pole branches. This formula covers simple poles with a nonzero ALP component and invertible meson operator. A zero-ALP pole requires a separately normalized meson eigenvector, as in the [earlier derivation](alp-explicit-finite-mixing-neutral-amplitude.md); degenerate poles require the full residue subspace. Overall pole-state signs are conventional. The physical amplitude below excludes the overall Feynman-rule factor $i$.

## Invariant source and inverse operator

The [dressed mass source](../shared/alp-general-chiral-redefinition-invariants.md) and [dressed topological source](../shared/alp-dressed-topological-invariant.md) are $\widehat M(x)=W_R^\dagger(x)M(x)W_L(x)=\widehat S(x)+i\widehat P(x)$ and $\widehat\Theta(x)=\Theta(x)+\int_0^x\langle C_A(t)\rangle dt$, where $x=a/f_a$, $\partial_x W_{R,L}=iC_{R,L}W_{R,L}$, and $W_{R,L}(0)=\mathbf1$. Subscripts on the hatted sources denote derivatives at $x=0$, not coefficients of powers of $x$. This is a source representative in which both dressed connections vanish. The input conventions are $C_{R,n}=C_{V,n}+C_{A,n}$, $C_{L,n}=C_{V,n}-C_{A,n}$ and $\chi_n=2B_0(S_n+iP_n)$.

Differentiate $W_R^\dagger M W_L$ at $x=0$ and insert $M(0)=m_q$:
$$
\widehat M_1=S_1+iP_1-i(C_{V,0}+C_{A,0})m_q+i m_q(C_{V,0}-C_{A,0}).
$$
Separating its Hermitian scalar and pseudoscalar parts gives
$$
\widehat S_1=S_1+i[m_q,C_{V,0}],\qquad
\widehat P_1=P_1-\{m_q,C_{A,0}\},\qquad
\widehat\Theta_1=\Theta_1+\langle C_{A,0}\rangle.
$$

Using the requested notation for the established operator,
$$
\mathbb M_s(X)=B_0\{m_q,X\}+\frac{m_0^2}{N_f}\langle X\rangle\mathbf1-sX,
$$
its invariant argument is
$$
\mathcal C=B_0\widehat P_1-\frac{m_0^2}{N_f}\widehat\Theta_1\mathbf1
=B_0\big(P_1-\{m_q,C_{A,0}\}\big)-\frac{m_0^2}{N_f}\big(\Theta_1+\langle C_{A,0}\rangle\big)\mathbf1.
$$
For the linear source parametrization of `Alda:2026nrd`, $P_1=\{c_q^P,m_q\}$ and $\Theta_1=c_G$, this reduces exactly to its invariant
$$
\mathcal C=-B_0\{m_q,c_q^{AP}\}-\frac{m_0^2}{N_f}c_{GA}\mathbf1.
$$
The existing analytical inverse applies at each $s$; its regularity conditions are given in the [operator-inverse record](../shared/alp-lo-finite-mixing-operator-inverse.md). The general source changes its argument as shown and changes the scalar pole equation through higher source derivatives. It does not change the meson operator above.

The meson part of the quadratic equation at a pole with ALP coordinate $a_s$ is
$$
\mathbb M_s(\Phi_s)=\frac{2f_\pi}{f_a}\mathcal C\,a_s.
$$
Where $\mathbb M_s$ is invertible, setting $a_s=1$ and applying the inverse gives
$$
\Phi_s=\frac{2f_\pi}{f_a}\mathbb M_s^{-1}(\mathcal C).
$$
The remaining quadratic equation determines the pole position. With $\epsilon=f_\pi/f_a$ and $m_{a,0}^2$ the coefficient from the independent quadratic ALP potential, eliminating $\Phi_s$ gives
$$
0=s-m_{a,0}^2-\epsilon^2\left[-B_0\langle\widehat S_2\rangle+\frac{2m_0^2}{N_f}\widehat\Theta_1^2\right]
+2\epsilon^2\langle\mathcal C\,\mathbb M_s^{-1}(\mathcal C)\rangle.
$$
The dressed kinetic term is canonical. With the meson normalization $\langle\lambda_i\lambda_j\rangle=2\delta_{ij}$, the squared norm of this unnormalized direction is
$$
1+\frac12\langle\Phi_s^2\rangle
=1+2\frac{f_\pi^2}{f_a^2}\left\langle\big(\mathbb M_s^{-1}(\mathcal C)\big)^2\right\rangle.
$$
Thus the ALP coordinate of the normalized pole direction is
$$
\left[1+2\frac{f_\pi^2}{f_a^2}\left\langle\big(\mathbb M_s^{-1}(\mathcal C)\big)^2\right\rangle\right]^{-1/2}.
$$
The inverse may be large near a meson pole. Consequently the explicit powers of $1/f_a$ in the result below do not define a weak-mixing expansion.

## Derivation of the physical amplitude

1. In the dressed representative the derivative connections vanish. Use the CCWZ building blocks $u_\mu$ and $\chi_+$, with $u=\exp(i\Phi/2f_\pi)$ and Hermitian $\Phi$:
$$
u_\mu=i(u^\dagger\partial_\mu u-u\partial_\mu u^\dagger),\qquad
\chi_+=2B_0\left[u^\dagger(\widehat S+i\widehat P)u^\dagger+u(\widehat S-i\widehat P)u\right].
$$
The terms of the LO Lagrangian needed here are
$$
\mathcal L_{\rm LO}=\frac{f_\pi^2}{4}\langle u_\mu u^\mu+\chi_+\rangle
-\frac{f_\pi^2m_0^2}{N_f}\left(\widehat\Theta-\frac{i}{2}\langle\log u-\log u^\dagger\rangle\right)^2+\mathcal L_a.
$$
The independent ALP potential $\mathcal L_a$ has no cubic term by assumption. Expanding the CCWZ blocks and using cyclicity of the flavour trace gives
$$
\begin{aligned}
u_\mu&=-\frac{\partial_\mu\Phi}{f_\pi}+O(\Phi^3),\\
\langle\chi_+\rangle&=4B_0\left\langle\widehat S+\frac{\widehat P\Phi}{f_\pi}-\frac{\widehat S\Phi^2}{2f_\pi^2}-\frac{\widehat P\Phi^3}{6f_\pi^3}\right\rangle+O(\Phi^4),\\
\widehat\Theta-\frac{i}{2}\langle\log u-\log u^\dagger\rangle&=\widehat\Theta+\frac{\langle\Phi\rangle}{2f_\pi}.
\end{aligned}
$$
The kinetic term has no cubic contribution, and $\widehat P_0=0$ removes the source-independent $\Phi^3$ term. Insert $\widehat S(x)=\sum_n x^n\widehat S_n/n!$ and the analogous series for $\widehat P$ and $\widehat\Theta$, with $x=a/f_a$. The terms of total degree three in $a,\Phi$ are
$$
\begin{aligned}
\mathcal L_{\rm LO}^{(3)}={}&-\frac{B_0}{2f_a}a\langle\widehat S_1\Phi^2\rangle
+\frac{f_\pi B_0}{2f_a^2}a^2\langle\Phi\widehat P_2\rangle
-\frac{f_\pi m_0^2\widehat\Theta_2}{2N_f f_a^2}a^2\langle\Phi\rangle\\
&+\frac{f_\pi^2B_0}{6f_a^3}a^3\langle\widehat S_3\rangle
-\frac{f_\pi^2m_0^2}{N_f f_a^3}a^3\widehat\Theta_1\widehat\Theta_2.
\end{aligned}
$$
These are derivatives of the full dressed functions, so in general $\widehat S_3$, $\widehat P_2$ and $\widehat\Theta_2$ cannot be inferred from the first-order invariant $\mathcal C$ alone. The [second-order source record](../shared/alp-dressed-mass-second-order-invariants.md) and [general source record](../shared/alp-general-chiral-redefinition-invariants.md) specify their construction. No expansion in the mixing has been made.

2. For the pole-aligned calculation, apply the exact symmetric axial map
$$
g_R=\exp[-ix\mathbb M_{s_a}^{-1}(\mathcal C)],\qquad g_L=g_R^\dagger,\qquad u'=g_R u h^\dagger=h u g_L^\dagger,
$$
where $h$ is the CCWZ compensator. Its linear field relation is
$$
\Phi'=\Phi-\frac{2f_\pi}{f_a}\mathbb M_{s_a}^{-1}(\mathcal C)a.
$$
The parent direction found above therefore has zero meson component in the primed frame. The map has no quadratic field term, so at a stationary vacuum its nonlinear terms do not modify the cubic contact. The individual terms of this contact are not separate decay amplitudes.

The relevant transformed derivatives in that frame are
$$
\begin{aligned}
C'_{V,0}&=0,& C'_{A,0}&=-\mathbb M_{s_a}^{-1}(\mathcal C),& C'_{A,1}&=0,& S'_1&=\widehat S_1,\\
P'_2&=\widehat P_2-2\{\mathbb M_{s_a}^{-1}(\mathcal C),\widehat S_1\},\\
S'_3&=\widehat S_3+3\{\mathbb M_{s_a}^{-1}(\mathcal C),\widehat P_2\}
-3\{\mathbb M_{s_a}^{-1}(\mathcal C),\{\mathbb M_{s_a}^{-1}(\mathcal C),\widehat S_1\}\},\\
\Theta'_1&=\widehat\Theta_1+\langle\mathbb M_{s_a}^{-1}(\mathcal C)\rangle,&
\Theta'_2&=\widehat\Theta_2.
\end{aligned}
$$
Each daughter's transformed meson matrix is $(2f_\pi/f_a)[\mathbb M_{s_n}^{-1}(\mathcal C)-\mathbb M_{s_a}^{-1}(\mathcal C)]$. Substitution of these expressions into the transformed cubic contact cancels the parent-specific differences and yields the symmetric expression below.

Indeed, using the normalizations shown above, the primed pole directions satisfy
$$
 a=\sum_r\frac{a_r}{\sqrt{1+2\epsilon^2\langle[\mathbb M_{s_r}^{-1}(\mathcal C)]^2\rangle}},\qquad
 \Phi'=2\epsilon\sum_r\frac{[\mathbb M_{s_r}^{-1}(\mathcal C)-\mathbb M_{s_a}^{-1}(\mathcal C)]a_r}{\sqrt{1+2\epsilon^2\langle[\mathbb M_{s_r}^{-1}(\mathcal C)]^2\rangle}}.
$$
Inserting these into the inverse linear relation gives
$$
 \Phi=\Phi'+2\epsilon\mathbb M_{s_a}^{-1}(\mathcal C)a
 =2\epsilon\sum_r\frac{\mathbb M_{s_r}^{-1}(\mathcal C)a_r}{\sqrt{1+2\epsilon^2\langle[\mathbb M_{s_r}^{-1}(\mathcal C)]^2\rangle}}.
$$
Thus the cubic contact in the primed variables equals the displayed dressed cubic polynomial evaluated on these restored directions. The parent-dependent differences cancel before the pole-field coefficient is extracted.

3. To see the substitution explicitly, insert the sum of the three pole fields into the dressed cubic expression using the normalized directions derived above:
$$
a=\sum_{r\in\{a,i,j\}}\left[1+2\frac{f_\pi^2}{f_a^2}\left\langle\big(\mathbb M_{s_r}^{-1}(\mathcal C)\big)^2\right\rangle\right]^{-1/2}a_r,
$$
$$
\Phi=\frac{2f_\pi}{f_a}\sum_{r\in\{a,i,j\}}\left[1+2\frac{f_\pi^2}{f_a^2}\left\langle\big(\mathbb M_{s_r}^{-1}(\mathcal C)\big)^2\right\rangle\right]^{-1/2}\mathbb M_{s_r}^{-1}(\mathcal C)\,a_r.
$$
Differentiate the substituted fields with respect to three independently labelled pole fields. The product rule gives the following three identities, with $r,t\in\{a,i,j\}$:
$$
\begin{aligned}
\partial_{a_a}\partial_{a_i}\partial_{a_j}\langle\widehat S_1 a\Phi^2\rangle\big|_0
&=\frac{4f_\pi^2}{f_a^2}\prod_r\left[1+2\frac{f_\pi^2}{f_a^2}\langle(\mathbb M_{s_r}^{-1}\mathcal C)^2\rangle\right]^{-1/2}
\sum_{r<t}\left\langle\widehat S_1\{\mathbb M_{s_r}^{-1}\mathcal C,\mathbb M_{s_t}^{-1}\mathcal C\}\right\rangle,\\
\partial_{a_a}\partial_{a_i}\partial_{a_j}(a^2\Phi)\big|_0
&=\frac{4f_\pi}{f_a}\prod_r\left[1+2\frac{f_\pi^2}{f_a^2}\langle(\mathbb M_{s_r}^{-1}\mathcal C)^2\rangle\right]^{-1/2}\sum_r\mathbb M_{s_r}^{-1}(\mathcal C),\\
\partial_{a_a}\partial_{a_i}\partial_{a_j}a^3\big|_0
&=6\prod_r\left[1+2\frac{f_\pi^2}{f_a^2}\langle(\mathbb M_{s_r}^{-1}\mathcal C)^2\rangle\right]^{-1/2}.
\end{aligned}
$$
Inserting these derivatives into the five cubic terms and adding them constructs the physical amplitude:
$$
\mathcal M(a_{\rm phys}\to P_iP_j)=\left.\partial_{a_a}\partial_{a_i}\partial_{a_j}\mathcal L_{\rm LO}^{(3)}\right|_{a_r=0}.
$$
The derivatives of $a^2\langle\Phi\widehat P_2\rangle$ and $a^2\langle\Phi\rangle$ use the second identity; the last two terms use the third. The result, with the common normalization factored out, is
$$
\begin{aligned}
\mathcal M(a_{\rm phys}\to P_iP_j)={}&\frac{f_\pi^2}{f_a^3}
\prod_{r\in\{a,i,j\}}\left[1+2\frac{f_\pi^2}{f_a^2}\left\langle\big(\mathbb M_{s_r}^{-1}(\mathcal C)\big)^2\right\rangle\right]^{-1/2}\\
&\times\Bigg\{-2B_0\sum_{r<t}\left\langle\widehat S_1\left\{\mathbb M_{s_r}^{-1}(\mathcal C),\mathbb M_{s_t}^{-1}(\mathcal C)\right\}\right\rangle\\
&\qquad+2B_0\left\langle\widehat P_2\sum_r\mathbb M_{s_r}^{-1}(\mathcal C)\right\rangle
-\frac{2m_0^2}{N_f}\widehat\Theta_2\left\langle\sum_r\mathbb M_{s_r}^{-1}(\mathcal C)\right\rangle\\
&\qquad+B_0\langle\widehat S_3\rangle-\frac{6m_0^2}{N_f}\widehat\Theta_1\widehat\Theta_2\Bigg\}.
\end{aligned}
$$
The pair sum runs over the three distinct unordered pairs from $\{a,i,j\}$, and the other sums run over all three poles. This is symmetric under exchange of the daughters. For identical daughters, no extra factor enters the amplitude; the phase space carries $1/2!$. The formula includes every ALP admixture of the daughters, and the full source dependence needed by the LO cubic contact. Its phase is fixed by choosing the ALP component of each pole positive.

## Verification

The [CCWZ building-block derivation](../../../scripts/derive_alp_ccwz_arbitrary_order.py) expands the defining products for $u_\mu$ and $\chi_+$. The [pole-amplitude derivation](../../../scripts/derive_alp_invariant_pole_aligned_amplitude.py) differentiates the cubic polynomial and prints the five contributions to the physical amplitude. Its [unit test](../../../scripts/test_derive_alp_invariant_pole_aligned_amplitude.py) compares the compact expression with the CCWZ cubic generator and the pole-aligned contact for noncommuting Hermitian $U(3)$ sources. The [source-consistent mixed-pole tests](../../../scripts/test_derive_alp_explicit_neutral_amplitude.py) compare the two routes using the same quadratic and cubic sources; their on-shell scalar products are algebraic and need not lie in real decay phase space. General nonconstant source transport is not separately compared.

```sh
.venv/bin/python scripts/derive_alp_invariant_pole_aligned_amplitude.py
.venv/bin/python scripts/derive_alp_ccwz_arbitrary_order.py --order 1 --mesons 3
.venv/bin/python -m unittest discover -s scripts -p test_derive_alp_invariant_pole_aligned_amplitude.py
.venv/bin/python -m unittest discover -s scripts -p test_derive_alp_explicit_neutral_amplitude.py
```

## Leading weak-mixing limit

Take $f_a\gg f_\pi$ at fixed $f_\pi$, quark masses, bare ALP mass and source derivatives, and write $\epsilon=f_\pi/f_a$. Assume the unmixed ALP pole is separated from the relevant meson eigenvalues by fixed nonzero gaps. Choose normalized Hermitian eigenmatrices of the unmixed meson operator:
$$
\mathbb M_0(\lambda_n)=m_{n,0}^2\lambda_n,\qquad
\langle\lambda_n\lambda_m\rangle=2\delta_{nm},\qquad
\mathbb M_0(X)=B_0\{m_q,X\}+\frac{m_0^2}{N_f}\langle X\rangle\mathbf1.
$$
The eigenvalues $m_{n,0}^2$ describe the unmixed meson operator; they are not the physical finite-mixing pole masses. Expanding $\mathcal C$ in this orthonormal basis and applying $\mathbb M_s^{-1}$ gives
$$
\mathcal C=\frac12\sum_n\lambda_n\langle\lambda_n\mathcal C\rangle,\qquad
\mathbb M_s^{-1}(\mathcal C)=\frac12\sum_n
\frac{\lambda_n\langle\lambda_n\mathcal C\rangle}{m_{n,0}^2-s}.
$$
For a simple meson eigenvalue with $\langle\lambda_i\mathcal C\rangle\ne0$, insert the singular $n=i$ term into the scalar pole equation above and write $s_i=m_{i,0}^2+\delta s_i$. The leading equation and its solution are
$$
0=m_{i,0}^2-m_{a,0}^2-\frac{\epsilon^2\langle\lambda_i\mathcal C\rangle^2}{\delta s_i}+O(\epsilon^2),\qquad
\delta s_i=\frac{\epsilon^2\langle\lambda_i\mathcal C\rangle^2}{m_{i,0}^2-m_{a,0}^2}+O(\epsilon^4).
$$
Substitution into the pole direction and its normalization yields
$$
\frac{2\epsilon\mathbb M_{s_i}^{-1}(\mathcal C)}
{\sqrt{1+2\epsilon^2\langle[\mathbb M_{s_i}^{-1}(\mathcal C)]^2\rangle}}
=\operatorname{sgn}\!\left[(m_{a,0}^2-m_{i,0}^2)\langle\lambda_i\mathcal C\rangle\right]\lambda_i+O(\epsilon^2).
$$
The sign follows from choosing a positive ALP component for each mixed pole. It can be absorbed into the phase of the daughter state. The parent ALP component is $1+O(\epsilon^2)$, its meson direction is $O(\epsilon)$, and a meson-like daughter's ALP component is $O(\epsilon)$. Inserting these orders and the normalized daughter eigenmatrices into the cubic contact gives
$$
\mathcal M(a_{\rm phys}\to P_iP_j)
=-\frac{B_0}{2f_a}\left\langle\widehat S_1\{\lambda_i^\dagger,\lambda_j^\dagger\}\right\rangle+O(\epsilon^3).
$$
The dagger denotes the outgoing flavour direction and the phases have been chosen as in the earlier [first-order two-pseudoscalar result](alp-flavour-violating-two-pseudoscalar-amplitudes.md). A zero overlap cannot be normalized through the ALP component; its normalized eigenmatrix is used directly. In a degenerate meson eigenspace, choose normalized eigenmatrices adapted to the mixing source and its orthogonal complement, then use the full pole equations for the mixed direction. Near an ALP–meson degeneracy this expansion is not uniform. The remainder denotes order in $f_\pi/f_a$ at fixed dimensionful inputs.

The general source satisfies $\widehat S_1=S_1+i[m_q,C_{V,0}]$. In the linear-source convention of `Alda:2026nrd`, $S_1=i[c_q^S,m_q]$ and $C_{V,0}=c_q^V$, so
$$
\widehat S_1=-i[c_q^{VS},m_q].
$$
Substituting this into the preceding amplitude gives the earlier result, in its vertex-factor convention,
$$
i\mathcal M=-\frac{B_0}{2f_a}\left\langle
\{\lambda_i^\dagger,\lambda_j^\dagger\}[c_q^{VS},m_q]\right\rangle.
$$
To compare with the alternative expression in that record, first use trace cyclicity and the definition of $\mathbb M_0$:
$$
\begin{aligned}
B_0\left\langle\{\lambda_i^\dagger,\lambda_j^\dagger\}[c_q^{VS},m_q]\right\rangle
&=B_0\left\langle c_q^{VS}[m_q,\{\lambda_i^\dagger,\lambda_j^\dagger\}]\right\rangle\\
&=\left\langle c_q^{VS}\left([\mathbb M_0(\lambda_i^\dagger),\lambda_j^\dagger]-[\lambda_i^\dagger,\mathbb M_0(\lambda_j^\dagger)]\right)\right\rangle.
\end{aligned}
$$
The singlet contributions vanish separately because they commute with either generator. Applying $\mathbb M_0(\lambda_n^\dagger)=m_{n,0}^2\lambda_n^\dagger$ now gives
$$
B_0\left\langle\{\lambda_i^\dagger,\lambda_j^\dagger\}[c_q^{VS},m_q]\right\rangle
=(m_{i,0}^2-m_{j,0}^2)\left\langle c_q^{VS}[\lambda_i^\dagger,\lambda_j^\dagger]\right\rangle.
$$
The singlet term cancels in this identity. It reproduces the earlier result with *unmixed* eigenvalues, without identifying them with the physical finite-mixing masses. For $a\to K^0\pi^0$, in the isospin-symmetric limit, use $\lambda_{K^0}^\dagger=\sqrt2 E_{sd}$, $\lambda_{\pi^0}=\operatorname{diag}(1,-1,0)$ and the $d$–$s$ source entry. Direct evaluation gives
$$
i\mathcal M(a\to K^0\pi^0)=\frac{B_0(m_s-m_d)(c_q^{VS})_{ds}}{\sqrt2 f_a}+O(\epsilon^3).
$$
This is the earlier flavour-state limit, not a general finite-mixing amplitude. The neutral kaon flavour directions are degenerate in the unmixed strong mass operator, so this check combines the relevant real directions rather than identifying either one with $K^0$ alone.

The [weak-mixing spectrum derivation](../shared/alp-lo-weak-mixing-expansion.md) gives the pole asymptotics. The [SymPy derivation](../../../scripts/derive_alp_two_pseudoscalar_weak_limit.py) derives the eigenbasis inverse and pole shift in a two-state example, differentiates the cubic expression, and evaluates the mass-operator identity and the $K^0\pi^0$ result. Its [unit test](../../../scripts/test_derive_alp_two_pseudoscalar_weak_limit.py) checks these identities and the leading power counting:

```sh
.venv/bin/python scripts/derive_alp_two_pseudoscalar_weak_limit.py
.venv/bin/python -m unittest discover -s scripts -p test_derive_alp_two_pseudoscalar_weak_limit.py
```
