---
title: Vector three-pion contact and its IBP/EOM reduction
status: derived
created: 2026-10-26
updated: 2026-10-26
bibtex_keys: [GomezDumm:2010ra, Alda:2026nrd, Ecker:1988te, Ecker:1989yg]
provenance: analytical-derivation
---

# Assumptions and inputs

This record tests whether the local odd-parity structure
$$
\mathcal C_{V3\pi}=
i\epsilon_{\mu\nu\alpha\beta}
\left\langle(\nabla_\gamma V^{\mu\nu})
\{u^\gamma,[u^\alpha,u^\beta]\}\right\rangle
$$
is independent modulo total derivatives and the leading Goldstone equation
of motion (EOM). External sources are set to zero after forming the
covariant expressions. The basis comparison uses the five local $VPPP$
operators of [@GomezDumm:2010ra], whose construction uses the leading
Goldstone EOM. The CCWZ and trace conventions follow the repository's
[@Alda:2026nrd] conventions.

The candidate is a symmetry-allowed local odd-parity resonance structure:
it is C-even, Hermitian after the displayed factor of $i$, and has the
required Levi-Civita tensor. At cubic order, use
$u_\mu=-\partial_\mu\Phi/f_\pi+O(\Phi^3)$ and the fixed off-shell tensor
projection
$$
V_i^{\mu\nu}=\frac{q^\mu\varepsilon_i^\nu-q^\nu\varepsilon_i^\mu}{M_i},
\qquad q\cdot\varepsilon_i=0.
$$
For $T_i=\operatorname{diag}(t_n,t_n,t_s)$ and $s=q^2$, the symbolic
contraction gives $-32i s t_n/(M_i f_\pi^3)$ before Fourier factors,
multiplying $\epsilon(p_+,p_-,p_0,\varepsilon_i)$. The vertex vanishes for
ideal $\phi$ and the exact-isospin $\rho^0$ projection, but is nonzero for
an isoscalar with a nonzero nonstrange generator component.

For definiteness, the leading $SU(N_f)$ EOM used here is
$$
\nabla_\mu u^\mu=\frac{i}{2}\left(
\chi_- -\frac{\langle\chi_-\rangle}{N_f}\right).
$$
The subtraction is important for a traceless $SU(N_f)$ Goldstone field.
This is not the $U(3)$ nonet EOM: a singlet field and its topological
potential add terms and must be treated separately.

# Operator basis and reduction

With $h^{\alpha\gamma}=\nabla^\alpha u^\gamma+
\nabla^\gamma u^\alpha$, the relevant members of the published basis are
$$
\begin{aligned}
\mathcal O_2={}&i\epsilon_{\mu\nu\alpha\beta}
\left\langle V^{\mu\nu}\left(
h^{\alpha\gamma}u^\beta u_\gamma-
u_\gamma u^\beta h^{\alpha\gamma}\right)\right\rangle,\\
\mathcal O_3={}&i\epsilon_{\mu\nu\alpha\beta}
\left\langle V^{\mu\nu}\left(
u_\gamma h^{\alpha\gamma}u^\beta-
u^\beta h^{\alpha\gamma}u_\gamma\right)\right\rangle,\\
\mathcal O_4={}&\epsilon_{\mu\nu\alpha\beta}
\left\langle\{V^{\mu\nu},u^\alpha u^\beta\}\chi_-\right\rangle.
\end{aligned}
$$
Applying covariant integration by parts to move $\nabla_\gamma$ off the
vector tensor, using the zero-source Maurer--Cartan relation to exchange
antisymmetrized derivatives of $u_\mu$, and then applying the displayed
EOM gives, modulo a total derivative,
$$
\boxed{\mathcal C_{V3\pi}\simeq
-\mathcal O_2-\mathcal O_3+\mathcal O_4
-\frac{2}{N_f}\epsilon_{\mu\nu\alpha\beta}
\langle\chi_-\rangle
\langle V^{\mu\nu}u^\alpha u^\beta\rangle.}
$$
The last term is a double-trace EOM remainder. Thus the candidate is not a
new independent single-trace $VPPP$ structure: it lies in the span of
$\mathcal O_2,\mathcal O_3,\mathcal O_4$ plus the displayed double trace.
It is not zero: its nonzero on-shell vertex is represented by the
corresponding combination of basis couplings.

In the isoscalar, pion-only projection, $\chi_-$ has no trace component at
the relevant three-pion order, so the double-trace remainder vanishes. In
that restricted channel,
$$
\mathcal C_{V3\pi}\simeq-\mathcal O_2-\mathcal O_3+\mathcal O_4.
$$
The shorter relation obtained from a nonet EOM that omits the trace
subtraction must not be assumed for a $U(3)$ theory with a singlet
potential.

# Reproducibility

The executable
[derivation script](../../../scripts/derive_wzw_gamma_three_pion_and_vector_contact.py)
performs the covariant IBP/Maurer--Cartan reduction, projects the leading
EOM, and solves for the coefficients in the five-operator basis plus the
double-trace remainder. Run:

    .venv/bin/python scripts/derive_wzw_gamma_three_pion_and_vector_contact.py

Six focused SymPy checks pass. They verify the stated reduction under the
listed assumptions; they do not fix the physical coefficients of the
independent $VPPP$ basis operators.
