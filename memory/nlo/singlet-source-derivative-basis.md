---
title: Singlet-source derivative operators in Alda notation
status: derived
created: 2026-11-02
updated: 2026-11-05
bibtex_keys: [Kaiser:1998ds,Herrera-Siklody:1996tqr,  Alda:2026nrd]
provenance: analytical-derivation
---

# Singlet-source derivative operators in Alda notation

## Assumptions and inputs

Use the Kaiser singlet derivatives and their field convention from [the source extract](../extracts/Kaiser:1998ds/singlet-source-derivative-terms.md), and Alda:2026nrd's definitions of $X_{TC}$ and $u_\mu$ from [the topological-source record](../extracts/Alda:2026nrd/lo-topological-action.md) and [Appendix A](../extracts/Alda:2026nrd/appendix-a-expansion.md). Write $Q_\mu=\langle u_\mu\rangle$ and $X_\mu=\partial_\mu X_{TC}$. Lorentz indices in products below are contracted. The derivative basis is

$$
\mathcal L_{\rm singlet,deriv}=V_3D_\mu\psi D^\mu\psi+V_4D_\mu\psi D^\mu\theta+V_5D_\mu\theta D^\mu\theta.
$$

The covariant-derivative identities, in the project's conventions, are

$$
D_\mu\psi=-Q_\mu,\qquad D_\mu\theta=Q_\mu-X_\mu.
$$

The second follows from Alda:2026nrd's $\partial_\mu X_{TC}=Q_\mu-2c_{GA}\partial_\mu a/f_a$ and the source relation $D_\mu\theta=2c_{GA}\partial_\mu a/f_a$ [@Kaiser:1998ds; @Alda:2026nrd].

## Conversion and field convention

Substituting gives

$$
\begin{aligned}
V_3D_\mu\psi D^\mu\psi&=V_3 Q_\mu Q^\mu,\\
V_4D_\mu\psi D^\mu\theta&=V_4\left(X_\mu Q^\mu-Q_\mu Q^\mu\right),\\
V_5D_\mu\theta D^\mu\theta&=V_5\left(X_\mu X^\mu-2X_\mu Q^\mu+Q_\mu Q^\mu\right).
\end{aligned}
$$

Thus the complete expression in the Alda operator basis is

$$
\mathcal L_{\rm singlet,deriv}
=V_3Q_\mu Q^\mu
+V_4(X_\mu Q^\mu-Q_\mu Q^\mu)
+V_5(X_\mu-Q_\mu)(X^\mu-Q^\mu).
$$

Kaiser:1998ds removes the $V_4$ structure by the allowed meson-field redefinition $U\to U\exp[i f(\psi+\theta)]$. In Alda notation, the removable structure is therefore the whole combination $X_\mu Q^\mu-Q_\mu Q^\mu$, not $X_\mu Q^\mu$ alone. Removing it shifts the coefficients assigned to the other singlet derivative terms. The coefficient of the standalone cross monomial is not itself invariant under this change of field coordinates.

After choosing $V_4=0$, the $V_5$ term remains as the fixed square

$$
V_5(X_\mu-Q_\mu)(X^\mu-Q^\mu)=V_5D_\mu\theta D^\mu\theta.
$$

It is not a standalone $X_\mu X^\mu$ term: when expanded in Alda monomials it includes the cross term and the $Q_\mu Q^\mu$ term with correlated coefficients. Alda:2026nrd's Eq. (III.16) gives

$$
Q_\mu-X_\mu=2c_{GA}\frac{\partial_\mu a}{f_a},
$$

so for a dynamical ALP this remaining square is

$$
4V_5c_{GA}^2\frac{\partial_\mu a\partial^\mu a}{f_a^2}.
$$

It has no meson–ALP kinetic mixing at $O(1/f_a)$; it first changes the ALP kinetic normalization at $O(1/f_a^2)$. With an external $\theta$ source, the corresponding source-derivative contact term remains and may matter for source correlators. Thus eliminating $V_4$ does not eliminate the source-derivative-squared contribution, but it must not be treated as an isolated $(\partial X_{TC})^2$ monomial: inside the retained $V_5$ square its coefficient is tied to the cross and $Q_\mu Q^\mu$ terms.

## Verification

The script output explicitly gives the $V_4$ operator as $X_\mu Q^\mu-Q_\mu Q^\mu$, the remaining $V_5$ operator as $(X_\mu-Q_\mu)^2$, and its ALP-source form as $4V_5c_{GA}^2(\partial a)^2/f_a^2$. The script runs with the project environment.

## Scope beyond the retained quadratic basis

At the two-derivative level, the relevant singlet derivative structures are $D_\mu\psi D^\mu\psi$, $D_\mu\psi D^\mu\theta$, and $D_\mu\theta D^\mu\theta$. At order $p^4$, there are further source-derivative operators in the unrestricted $U(3)$ basis. Kaiser:1998ds describes 57 potentials before the reduction to the subset used for masses and decay constants; its $O_{46}$ and $O_{52}$ source-derivative terms are removable by a higher-order meson-field redefinition [the source extract](../extracts/Kaiser:1998ds/singlet-source-derivative-terms.md) [@Kaiser:1998ds]. Other terms in the full list involve additional meson derivatives or external field strengths, and their relevance depends on the observables and sources being retained. Therefore, within the project's retained quadratic NLO basis the relevant source-derivative contribution is the $V_5$ square, while a claim about the complete $p^4$ generating functional requires specifying the source and observable content.

## The $O_{46}$, $O_{47}$, $O_{50}$, and $O_{52}$ source terms

For the operators in the full source-dependent $p^4$ basis, let $C_\mu=U^\dagger D_\mu U$, $M=U^\dagger\chi+\chi^\dagger U$, and $\hat\theta=i\theta$. Their structures are listed in Appendix C of Herrera-Siklody:1996tqr [@Herrera-Siklody:1996tqr]:

$$
\begin{aligned}
O_{46}&=D_\mu\hat\theta\,\langle C^\mu M\rangle,&
O_{47}&=D_\mu\hat\theta\,\langle C^\mu\rangle\langle M\rangle,\\
O_{50}&=D_\mu\hat\theta D^\mu\hat\theta\,\langle M\rangle,&
O_{52}&=iD_\mu D^\mu\hat\theta\,\langle M\rangle.
\end{aligned}
$$

In Alda notation, $D_\mu\theta=Q_\mu-X_\mu=2c_{GA}\partial_\mu a/f_a$, with $Q_\mu=\langle u_\mu\rangle$ and $X_\mu=\partial_\mu X_{TC}$ [@Alda:2026nrd]. At the real-mass vacuum, $M$ has a constant term proportional to $\chi_0$, and the linear meson part of $\langle M\rangle$ vanishes. The executable expansion in [the quadratic-operator script](../../scripts/derive_theta_operator_quadratic_terms.py) gives

$$
\begin{aligned}
O_{46}^{(a\Phi)}&=-\frac{4c_{GA}}{f_af_\pi}\partial_\mu a\,
\langle\partial^\mu\Phi\,\chi_0\rangle,\\
O_{47}^{(a\Phi)}&=-\frac{4c_{GA}}{f_af_\pi}\partial_\mu a\,
\langle\partial^\mu\Phi\rangle\langle\chi_0\rangle,\\
O_{50}^{(a\Phi)}&=0.
\end{aligned}
$$

The single-trace $O_{46}$ and the $O_{52}$ direction are the two structures removed in Kaiser:1998ds by

$$
U\longrightarrow U\exp\!\left[
f_1D_\mu\theta\,U^\dagger D^\mu U
+i f_2D_\mu D^\mu\theta
\cdots\right].
$$

The first term in the exponent changes the meson coordinates along the covariant vector $D_\mu\hat\theta\,C^\mu$ and removes the $O_{46}$ coefficient. The singlet phase shift proportional to $D^2\hat\theta$ removes $O_{52}$. These are changes of operator basis: the variation of the lower-order Lagrangian generates the terms that cancel $O_{46}$ and $O_{52}$, while shifting other local coefficients. They do not make the corresponding quadratic bilinears vanish in the original coordinates.

$O_{47}$ is a distinct double-trace structure and is not removed by these two redefinition directions. It gives a singlet kinetic mixing proportional to $\langle\chi_0\rangle$; it is absent from the retained basis only if the additional source-dependent operator is omitted by the chosen large-$N_c$/operator truncation. $O_{50}$ is not removed by the cited redefinitions either, but it contains two factors of $D_\mu\theta$ and hence starts at $O(1/f_a^2)$ for a dynamical ALP. It cannot contribute to the one-ALP kinetic mixing at $O(1/f_a)$; as an external-source operator it remains relevant to source contact terms.

## Explicit meson-field shift for the $O_{47}$ bilinear

The expressions below first show the quadratic expansion of $O_{47}$ with unit coefficient. In the convention $\mathcal L=\sum_i L_iO_i$, the contribution to the Lagrangian is multiplied by $L_{47}$:

$$
\mathcal L_{47}^{(a\Phi)}
=-\frac{4c_{GA}}{f_af_\pi}\,
\partial_\mu a\,\langle\partial^\mu\Phi\rangle\langle\chi_0\rangle.
$$

With $\Phi=\lambda_3\pi_3+\lambda_8\eta_8+\lambda_0\eta_0$ and $\lambda_0=\sqrt{2/3}\,\mathbf1$, this is

$$
\mathcal L_{47}^{(a\Phi)}
=-\frac{4\sqrt6\,c_{GA}\langle\chi_0\rangle}{f_af_\pi}
\partial_\mu\eta_0\,\partial^\mu a.
$$

Accordingly, the mixed kinetic column $\Delta K_i^{(47)}$ includes the factor $L_{47}$. Under

$$
P_i=P_i'-\rho_i a,
$$

the blocks transform as $K'_{i\hat a}=K_{i\hat a}-K_{ij}\rho_j$ and $M_{i\hat a}^{2\prime}=M^2_{i\hat a}-M^2_{ij}\rho_j$. Thus choosing

$$
\rho_i=(K^{-1})_{ij}\Delta K_j^{(47)}
$$

removes only the $O_{47}$ addition to the kinetic block. In the leading canonical meson basis $K_{ij}=\delta_{ij}$, this is a pure singlet shift with $\rho_0=-4\sqrt6\,L_{47}c_{GA}\langle\chi_0\rangle/(f_af_\pi)$ and $\rho_3=\rho_8=0$. With a nontrivial meson kinetic matrix, $K^{-1}$ can spread the shift over the neutral basis.

The contribution is transferred to mass mixing,

$$
\Delta M_{i\hat a}^{2\prime}=-M^2_{ij}\rho_j,
$$

and to any direct linear source $\mathcal J=V_{\hat a}a+V_iP_i$, whose coefficient becomes $V'_{\hat a}=V_{\hat a}-V_i\rho_i$. The invariant source is unchanged:

$$
M^2K^{-1}(K_{\cdot\hat a}-K\rho)
-(M^2_{\cdot\hat a}-M^2\rho)
=M^2K^{-1}K_{\cdot\hat a}-M^2_{\cdot\hat a}.
$$

This shift removes the $O_{47}$ kinetic mixing from the chosen coordinates; it does not remove the operator's contribution from the full theory. The executable calculation is in [the O47 field-shift script](../../scripts/derive_o47_field_shift.py); run it with .venv/bin/python scripts/derive_o47_field_shift.py.

There is a coefficient-function subtlety for $O_{52}$. Its parity-even coefficient is odd in $X_{TC}$, so near the real-mass vacuum it starts as $F_{52}(X_{TC})=F'_{52}(0)X_{TC}+\cdots$. After integration by parts,

$$
F'_{52}(0)X_{TC}D^2\theta\,\langle M\rangle
\doteq
F'_{52}(0)D_\mu\theta\left[
(\partial^\mu X_{TC})\langle M\rangle
+X_{TC}\partial^\mu\langle M\rangle
\right],
$$

where $\doteq$ denotes equality up to a total derivative, with signs understood in the $O_{52}$ convention above. Using $X_\mu=Q_\mu-D_\mu\theta$, the first term shifts the $O_{47}$ and $O_{50}$ coefficients together; the second has no $a\Phi$ bilinear about the real-mass vacuum. Thus setting $O_{52}$ to zero does not justify setting $O_{47}$ to zero. The retained quadratic-basis decision about $O_{47}$ must be made separately.

The SymPy script expands the four operators about the real-mass vacuum and checks their field and $1/f_a$ orders. It also displays the integration-by-parts reshuffling of the parity-even $O_{52}$ coefficient. Run .venv/bin/python scripts/derive_theta_operator_quadratic_terms.py.

### Notation

The building blocks used in Herrera-Siklody:1996tqr adapted to the notation of Alda:2026nrd are

$$
\begin{aligned}
C_\mu &= U^\dagger D_\mu U = -u u^\dagger u_\mu u\\
M &= U^\dagger \chi + \chi^\dagger U = u^\dagger \chi_+ u\\
N &= U^\dagger \chi - \chi^\dagger U = u^\dagger \chi_- u\\
D_\mu\hat\theta &= iD_\mu\theta = i(\langle u_\mu\rangle -\partial_\mu X_{TC}) = 2i c_{GA} \frac{\partial_\mu a}{f_a}
\end{aligned}
$$

### Operators

The operators that modify the kinetic terms must include exactly two derivatives and no field without a derivative. The operators from Herrera-Siklody:1996tqr fulfilling these conditions, and where at least one of the derivatives is the singlet $D_\mu\theta$, are:

$$
\begin{aligned}
O_{46} &= (\langle u_\mu\rangle - \partial_\mu X_{TC})\langle u^\mu \chi_+\rangle,\\
O_{47} &= (\langle u_\mu\rangle - \partial_\mu X_{TC})\langle u^\mu\rangle \langle\chi_+\rangle,\\
O_{50} &= -(\langle u_\mu\rangle - \partial_\mu X_{TC})(\langle u^\mu\rangle - \partial^\mu X_{TC})\langle\chi_+\rangle,\\
O_{52} &= -\nabla^\mu(\langle u_\mu\rangle - \partial_\mu X_{TC})\langle\chi_+\rangle.
\end{aligned}
$$

## Reproducibility

The SymPy script substitutes the covariant derivatives, expands each operator in the $(Q_\mu,X_\mu)$ basis, and applies the Alda source identity to the remaining $V_5$ square. Run .venv/bin/python scripts/derive_singlet_source_derivative_basis.py.
