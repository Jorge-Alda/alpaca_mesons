---
title: Conditional three-pion discontinuity and subtracted ALP TFF representation
status: derived
created: 2026-10-25
updated: 2026-10-25
bibtex_keys: [Alda:2026nrd, Messerli:2025rnv, Niecknig:2012sj]
provenance: analytical-derivation
---

# Scope and inputs

This derives the normalization and dispersive bookkeeping of the charged
three-pion cut. It does not supply the complete ALP production matrix element
or its numerical spectral weight. The inputs are the unnormalized scalar
and epsilon convention in the [tree TFF](alp-doubly-virtual-tff.md), the
strict [spin-one pole residue](alp-three-pion-vector-pole-factorization.md),
the [decay partial-wave equation](alp-three-pion-vector-partial-wave.md), and
the [numerical-input audit](alp-three-pion-numerical-input-status.md).
The source convention is Alda:2026nrd; the dispersive strategy is
Messerli:2025rnv and the decay equation is Niecknig:2012sj
[@Alda:2026nrd; @Messerli:2025rnv; @Niecknig:2012sj]. The bibliography
entries were inspected. No additional literature normalization is imported.

Use natural units, momenta in GeV, metric $(+---)$,
$\epsilon^{0123}=+1$, relativistic states
$\langle p'|p\rangle=2E_p(2\pi)^3\delta^3(\mathbf p'-\mathbf p)$,
and $S=1+iT$. For the decay crossing write $p=q_1+q_2$ and
$s_i=q_i^2$, with $p^2=m_a^2$. Reversing both photon momenta from the
all-incoming tree convention leaves its epsilon tensor unchanged. The user
specifies $2m_e<m_a<2m_\tau$ and the complete two-photon lepton loop;
the formula below therefore cannot be limited to a finite spacelike region.

# Unitarity and scalar extraction

After removing the momentum-conserving delta function, $S^\dagger S=1$
gives $T-T^\dagger=iT^\dagger T$. Strip one common electromagnetic
coupling from each one-current amplitude, as in
$\mathcal M(\gamma^*\to3\pi)=e\varepsilon_\mu
\langle3\pi,\mathrm{out}|J_{\rm em}^\mu|0\rangle$.
For boundary values connected by this physical cut, the charged three-pion
contribution is
$$
\epsilon^{\mu\nu\rho\sigma}q_{1\rho}q_{2\sigma}
\operatorname{Disc}_{s_1}F_{a\gamma^*\gamma^*}
=i\int d\Phi_3(q_1;p_+,p_-,p_0)
\langle\pi^+\pi^-\pi^0,\mathrm{out}|J_{\rm em}^\mu|0\rangle^*
\langle\pi^+\pi^-\pi^0,\mathrm{out}|J_{\rm em}^\nu|a\rangle.
$$
The transition matrix element on the right is crossed with the same
external-state convention as the two-current amplitude. In particular,
converting the stored $a\gamma^*\to3\pi$ amplitude into this matrix element
requires its crossing/LSZ phase and removal of its factor $e$. This equation
does not set the unresolved phase $\zeta_i$ to one. The complex conjugate
belongs to the vacuum-to-three-pion amplitude, not to the ALP coefficient.
Keep the full electromagnetic charge matrix on the ALP-current side.

The executable forms the epsilon tensor and contracts its two Lorentz
indices with the metric. It obtains
$$
\bigl(\epsilon^{\mu\nu\rho\sigma}q_{1\rho}q_{2\sigma}\bigr)
\bigl(\epsilon_{\mu\nu\alpha\beta}q_1^\alpha q_2^\beta\bigr)
=\frac{\lambda(m_a^2,s_1,s_2)}2,
\qquad
\lambda(a,b,c)=a^2+b^2+c^2-2ab-2ac-2bc.
$$
Thus the scalar discontinuity follows by contracting the right side with
$\epsilon_{\mu\nu\alpha\beta}q_1^\alpha q_2^\beta$ and multiplying
by $2/\lambda(m_a^2,s_1,s_2)$. This divides by a scalar Gram factor,
never by a tensor. Both Ward contractions vanish identically. At
$\lambda=0$ the projector is degenerate; define the scalar by continuation
from nondegenerate kinematics, after any cancellations in the complete
amplitude. A projector evaluated directly at that endpoint does not
determine the scalar.

Here $\operatorname{Disc}F=F(s_1+i0,s_2)-F(s_1-i0,s_2)$.
Replacing it with $2i\operatorname{Im}F$ additionally requires Schwarz
reflection for the other arguments and real source coefficients. With a
second timelike boundary value or general complex source coefficients, keep
the discontinuity and the matrix-element product explicitly.

# Phase space and spin-one weight

The intermediate state is one $\pi^+\pi^-\pi^0$ state with each momentum
integrated once. These species are distinct: there is no $1/3!$, and no
additional sum over six relabelings. This record evaluates that charged
channel; it does not use a charge multiplicity to represent other channels.
Starting from
$$
d\Phi_3=(2\pi)^4\delta^4(q_1-p_+-p_--p_0)
\prod_{r=+,-,0}\frac{d^3p_r}{(2\pi)^3 2E_r},
$$
the script factors the charged-pair invariant $u=(p_++p_-)^2$:
$$
d\Phi_3=\frac{du}{2\pi}
d\Phi_2(q_1;p_++p_-,p_0)
d\Phi_2(p_++p_-;p_+,p_-),\qquad
\frac{d\Phi_2}{d\Omega}=\frac{\sqrt{\lambda(s,m_1^2,m_2^2)}}{32\pi^2s}.
$$
Let $x=q_1^2$, $t=(p_++p_0)^2$ and $z$ be the angle between
$\mathbf p_+$ and $\mathbf p_0$ in the charged-pair rest frame. The
executed result, after integrating the overall orientation and internal
azimuth for an orientation-independent scalar integrand, is
$$
\frac{d\Phi_3}{du,dz}=
\frac{\sqrt{\lambda(x,u,m_0^2)\lambda(u,m_+^2,m_-^2)}}{256\pi^3xu},
\quad
\left|\frac{dt}{dz}\right|=
\frac{\sqrt{\lambda(x,u,m_0^2)\lambda(u,m_+^2,m_-^2)}}{2u},
\quad
\frac{d\Phi_3}{du,dt}=\frac1{128\pi^3x}.
$$
The domain is $x\geq(m_++m_-+m_0)^2$,
$(m_++m_-)^2\leq u\leq(\sqrt{x}-m_0)^2$ and $-1\leq z\leq1$,
with $t$ obtained from the displayed Jacobian and the charged-pair energies
in the script. For a general production amplitude the orientations relative
to $q_2$ must still be integrated; the Dalitz density alone does not perform
that angular integral.

For the restricted spin-one endpoint the existing geometric current obeys
$q_1\cdot J_{3\pi}=0$. The executable obtains
$$
-J_{3\pi}^2=\det[p_r\cdot p_s]_{r,s=+,-,0}
=\frac{(1-z^2)\lambda(x,u,m_0^2)
\lambda(u,m_+^2,m_-^2)}{16u}.
$$
If the scalar weight depends on Dalitz invariants only, its orientation
average is
$$
\left\langle J_{3\pi}^{\mu}J_{3\pi}^{\alpha}\right\rangle_{\rm orientation}
=\frac{-J_{3\pi}^2}{3}
\left(-g^{\mu\alpha}+\frac{q_1^\mu q_1^\alpha}{x}\right).
$$
The factor $1/3$ follows by the explicit angular integral and transverse
trace. It is a spin-one angular weight, not an identical-particle factor.
It can be used with the declared pole residue; it does not fix regular
off-shell production terms or direct contacts.

# Single and double spectral terms

Assume temporarily that one subtraction in each variable suffices, that the
required analytic continuation is known, and that additional or anomalous
cuts have either been included or stated separately. Cauchy subtraction
produces $1/(x-s)-1/x=s/[x(x-s)]$. Applied successively, it gives the
formal representation
$$
\begin{aligned}
F(s_1,s_2)={}&F(0,0)
+\frac{s_1}{2\pi i}\int\frac{dx\,\operatorname{Disc}_{x}F(x,0)}{x(x-s_1)}
+\frac{s_2}{2\pi i}\int\frac{dy\,\operatorname{Disc}_{y}F(0,y)}{y(y-s_2)}\\
&+\frac{s_1s_2}{(2\pi i)^2}\int\!dx\,dy\,
\frac{\operatorname{Disc}_{y}\operatorname{Disc}_{x}F(x,y)}
{xy(x-s_1)(y-s_2)}.
\end{aligned}
$$
Here $F$ is the existing scalar $F_{a\gamma^*\gamma^*}$, and the integrals
run over the actual discontinuity supports. For the charged-three-pion
$x$ cut its threshold is given above. The $y$ support need not be another
three-pion threshold: the other current can have isovector contributions.
All denominators acquire their specified timelike boundary values. The
unknown subtraction data and any further subtractions remain inputs; this
formula is conditional on its convergence assumption.

The two single-variable terms exchange under $s_1\leftrightarrow s_2$.
The script verifies that the double denominator exchanges under simultaneous
$(s_1,x)\leftrightarrow(s_2,y)$. The complete double density must obey that
exchange as well. A mixed density with one three-pion support and one
different support requires the exchanged term; a symmetric full density
already integrated on the complete domain must not be doubled. In
particular, the single unitarity relation above does not calculate
$\operatorname{Disc}_{y}\operatorname{Disc}_{x}F$ or its subtraction
functions. It does not justify importing a two-pion double-density factor.

For the elementary Cauchy integrand, the boundary-value identity
$1/(x-s-i0)-1/(x-s+i0)=2\pi i\delta(x-s)$ reconstructs the input
discontinuity. Applying it once to each denominator gives the double
discontinuity, provided the spectral inputs and analytic supports satisfy
the stated assumptions. This verifies the normalization of the conditional
representation, rather than establishing its unknown physical inputs.

# Explicit poles and continuum

A pole and its spectral representation describe the same strength. The
script verifies, for a general signed pole residue,
$$
\frac{R_i}{M_i^2-s}=
\frac{R_i}{M_i^2}+
\frac{R_i s}{M_i^2(M_i^2-s)}.
$$
The second term is the once-subtracted dispersion of
$\operatorname{Disc}[R_i/(M_i^2-x)]=2\pi i R_i\delta(x-M_i^2)$.
Consequently retaining the explicit pole requires removing this same pole
distribution from the spectral input and matching the subtraction value.
Alternatively, replace the pole by the spectral contribution and its
subtraction data. A finite-width three-pion solution containing the
$\omega,\phi$ strength cannot be added to unchanged stable poles. A
physical continuum/pole separation and its uncertainty are not selected
here.

# Reproducibility and remaining work

Run
```text
.venv/bin/python scripts/derive_alp_three_pion_discontinuity.py
.venv/bin/python scripts/derive_alp_three_pion_discontinuity.py --test
```
Both completed. Five unit tests verify independent massless phase-space
recursion ($\Phi_3=x/(256\pi^3)$), the general four-vector epsilon norm and
both Ward identities, the geometric-current Gram identity and
transversality, the spin-one average, and subtracted-pole/Bose identities.
The script prints each derived coefficient with SymPy LaTeX output.

Stage 3 is complete only for these conditional normalization and dispersion
relations. A numerical three-pion discontinuity still requires the complete
source-dependent production amplitude, the compatible complex current
amplitude and relative phase, and correlated variations. The double
discontinuity, subtraction information, possible extra cuts and high-energy
continuation for the full loop remain open. These are physical inputs or
additional derivations, not adjustable normalizations of the results above.
