---
title: Singlet-source derivative operators in Alda notation
status: derived
created: 2026-11-02
updated: 2026-11-02
bibtex_keys: [Kaiser:1998ds, Alda:2026nrd]
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

## Reproducibility

The SymPy script substitutes the covariant derivatives, expands each operator in the $(Q_\mu,X_\mu)$ basis, and applies the Alda source identity to the remaining $V_5$ square. Run .venv/bin/python scripts/derive_singlet_source_derivative_basis.py.
