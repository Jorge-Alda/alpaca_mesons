---
title: Isospin breaking in ALP production and K pi rescattering
status: mixed
created: 2026-09-23
updated: 2026-09-23
bibtex_keys:
  - Alda:2026nrd
  - Kubis:2001ij
  - Nehme:2001wf
  - Nehme:2001wa
  - Kubis:2001bx
  - Cirigliano:2001mk
  - Kastner:2008ch
  - Jamin:2001zq
provenance: analytical-derivation
---

# Isospin Breaking in ALP Production and K Pi Rescattering

## Scope and inputs

Isospin breaking enters the $a\to K\bar K\pi$ calculation in four distinct places:

- the ALP source, if $\widehat c_u\ne\widehat c_d$;
- strong QCD effects proportional to $m_d-m_u$, including neutral-meson mixing and local interactions;
- charged-neutral meson mass differences and therefore split thresholds;
- electromagnetic interactions and radiative corrections.

These effects should not be represented by a single correction to the isospin-symmetric $I=1/2$ Omnes function. The Alda:2026nrd-consistent source is obtained from the full flavor traces in [the leading-order interaction record](../alp/three-pseudoscalar/alp-phi3-lo-chipt-derivation.md) [@Alda:2026nrd]. Strong and electromagnetic corrections to low-energy $\pi K$ scattering are treated systematically in chiral perturbation theory by the sources recorded in [the Kubis--Meissner extract](../extracts/Kubis:2001ij/isospin-violation-pik.md), [the Nehme--Talavera extract](../extracts/Nehme:2001wf/low-energy-pik-breaking.md), and the charged-channel extracts [Nehme:2001wa](../extracts/Nehme:2001wa/charged-pik-breaking.md) and [Kubis:2001bx](../extracts/Kubis:2001bx/charged-pik-threshold.md). [@Kubis:2001ij; @Nehme:2001wf; @Nehme:2001wa; @Kubis:2001bx]

The charge-isospin transformations and threshold-mixing formulas below are derived by [`scripts/derive_kpi_isospin_breaking.py`](../../scripts/derive_kpi_isospin_breaking.py).

## Source isospin

For any diagonal light-quark coupling,

$$
c_u\bar uu+c_d\bar dd
=c_0(\bar uu+\bar dd)+c_3(\bar uu-\bar dd),
$$

where

$$
c_0=\frac{c_u+c_d}{2},
\qquad
c_3=\frac{c_u-c_d}{2}.
$$

The $c_0$ term is isoscalar and $c_3$ is isovector. For an isoscalar source and exact isospin, coupling the bachelor kaon with $I=1/2$ to total $I=0$ selects a $K\pi$ pair with $I=1/2$. An isovector source permits total $I=1$ and generally produces both $I=1/2$ and $I=3/2$ pair amplitudes. Therefore a source with $c_u\ne c_d$ must not be dressed entirely with the BaBar $I=1/2$ wave or with an $I=1/2$ Omnes function.

In practice, evaluate the Alda:2026nrd flavor traces without setting $u=d$. Keep the distinct entries of $\widehat c_q^A$, $\widehat c_q^P$, and $B_0m_q$, and use the physical neutral-meson generators. This gives separate production kernels for every physical charge channel. The physical generators and $B_0m_q$ inputs already used in `scripts/implementation/three_meson_width.py` provide the repository conventions for this calculation.

## Charge and isospin bases

Use the Condon--Shortley convention $|\pi^+\rangle=-|1,+1\rangle$. For the $I_3=+1/2$ sector, define

$$
\mathbf F_c=
\begin{pmatrix}
F_{K^+\pi^0}\\
F_{K^0\pi^+}
\end{pmatrix},
\qquad
\mathbf F_I=
\begin{pmatrix}
F_{1/2}\\
F_{3/2}
\end{pmatrix}.
$$

The SymPy Clebsch--Gordan calculation gives

$$
\boxed{
\mathbf F_c=C\mathbf F_I,
\qquad
C=
\begin{pmatrix}
1/\sqrt3&\sqrt{2/3}\\
\sqrt{2/3}&-1/\sqrt3
\end{pmatrix},
\qquad C^{-1}=C.
}
$$

In the isospin limit,

$$
T_c=C
\begin{pmatrix}
T^{1/2}&0\\
0&T^{3/2}
\end{pmatrix}C,
$$

so that

$$
T_c=
\begin{pmatrix}
\frac13T^{1/2}+\frac23T^{3/2}
&\frac{\sqrt2}{3}(T^{1/2}-T^{3/2})\\
\frac{\sqrt2}{3}(T^{1/2}-T^{3/2})
&\frac23T^{1/2}+\frac13T^{3/2}
\end{pmatrix}.
$$

The same construction applies to the charge-conjugate $I_3=-1/2$ block after adopting consistent antikaon and pion phases.

## Split thresholds

The two thresholds in the positive-charge pair sector are

$$
\sqrt{s_1}=M_{K^+}+M_{\pi^0}\simeq0.62865\,\mathrm{GeV},
$$

and

$$
\sqrt{s_2}=M_{K^0}+M_{\pi^+}\simeq0.63718\,\mathrm{GeV}.
$$

Their separation is approximately $8.53\,\mathrm{MeV}$. Between the thresholds, only $K^+\pi^0$ is open; above the second, both channels enter unitarity. Let

$$
\rho_c(s)=\operatorname{diag}(\rho_{K^+\pi^0},\rho_{K^0\pi^+}).
$$

Rotating physical phase space to the isospin basis gives

$$
C\rho_c C=
\begin{pmatrix}
\frac13\rho_1+\frac23\rho_2
&\frac{\sqrt2}{3}(\rho_1-\rho_2)\\
\frac{\sqrt2}{3}(\rho_1-\rho_2)
&\frac23\rho_1+\frac13\rho_2
\end{pmatrix}.
$$

Thus physical threshold splitting creates off-diagonal unitarity between nominal $I=1/2$ and $I=3/2$ components even if the short-range interaction is isospin diagonal. Two independent scalar Omnes functions cannot reproduce the resulting two-threshold cusp exactly.

## Charge-basis Omnes matrix

Let $T_{ij}$ be Coulomb-subtracted charge-channel partial waves and choose the normalization

$$
\rho_i(s)=
\frac{\lambda^{1/2}(s,M_{K_i}^2,M_{\pi_i}^2)}{16\pi s}
\theta(s-s_i).
$$

The production vector obeys

$$
\operatorname{disc}\mathbf F_c(s)
=2iT_c^\dagger(s)\rho_c(s)\mathbf F_c(s+i0).
$$

Define the charge-basis Omnes matrix by

$$
\operatorname{disc}\Omega_c(s)
=2iT_c^\dagger(s)\rho_c(s)\Omega_c(s+i0),
\qquad
\Omega_c(s_0)=\mathbb 1.
$$

In the independent-isobar approximation, the dressed source is

$$
\boxed{
\mathbf F_c(s)=\Omega_c(s)\mathbf P_c(s),
\qquad
\mathbf P_c(s)=C
\begin{pmatrix}P_{1/2}(s)\\P_{3/2}(s)\end{pmatrix}.
}
$$

For $a\to K^+K^-\pi^0$, rescattering $K^+\pi^0\to K^0\pi^+$ couples that final state to $K^-K^0\pi^+$. The charge-conjugate band similarly couples to $K^+\bar K^0\pi^-$. Isospin breaking therefore couples the charge channels discussed in [the charge-channel record](kpi-s-wave-in-alp-decays.md); it is not only a correction to one Dalitz band.

Above the $K\eta$ and $K\eta'$ thresholds, enlarge $\mathbf F_c$, $T_c$, and $\Omega_c$ to include $K\eta$ and $K\eta'$ charge channels. In the exact isospin limit these channels enter the $I=1/2$ sector; isospin breaking weakly connects the sectors. The isospin-symmetric three-channel construction is documented in [the coupled-form-factor extract](../extracts/Jamin:2001zq/coupled-scalar-form-factors.md). [@Jamin:2001zq]

## Strong isospin breaking

At leading chiral order, $m_d-m_u$ induces $\pi^0$--$\eta$ mixing as well as mass and vertex corrections. At higher orders it also enters loops, wave-function renormalization, and local counterterms. Including only a mixing angle or only physical masses is not a complete fixed-order treatment; separately identified effects can cancel. [@Kubis:2001ij; @Nehme:2001wf]

For a leading-order ALP analysis, using the fitted physical neutral generators and nondegenerate $B_0m_q$ values consistently in the Alda:2026nrd flavor traces captures the source-level effects included in those inputs. Dynamical isospin violation in the rescattering matrix remains omitted unless $T_c$ is matched to the complete isospin-breaking chiral amplitude.

Charge-dependent $K_{\ell3}$ scalar form factors provide useful low-energy normalization and slope checks because they include strong and electromagnetic isospin breaking. They are not universal ALP production form factors: the source vector, subtraction constants, and left-hand cuts differ. [See the scalar-form-factor extract](../extracts/Kastner:2008ch/kl3-scalar-isospin-breaking.md) [@Kastner:2008ch]

## Electromagnetism

The channels $K^+\pi^0$ and $K^0\pi^+$ contain only one charged particle each and have no two-body Coulomb ladder. Channels such as $K^+\pi^-$ contain two oppositely charged particles and require Coulomb subtraction or resummation near threshold. Virtual electromagnetic corrections are infrared divergent by themselves; a physical rate must include real soft photons with an explicit inclusiveness prescription or photon-energy cut. [@Nehme:2001wa; @Cirigliano:2001mk]

A pairwise Coulomb factor is only an approximation in a three-body decay because photons can connect particles assigned to different isobar pairs. Electromagnetic corrections should therefore be kept separate from the hadronic Omnes matrix and attached to a clearly defined infrared-safe observable.

## Practical hierarchy

### Minimal controlled improvement

1. Evaluate separate Alda:2026nrd production kernels for all physical charge channels with nondegenerate hatted couplings, physical neutral generators, and physical masses.
2. Extract $P_{1/2}$ and $P_{3/2}$, retaining the source contribution proportional to $c_u-c_d$.
3. Construct $T_c=C\operatorname{diag}(T^{1/2},T^{3/2})C$ from isospin-symmetric dispersive amplitudes.
4. Solve the $2\times2$ charge-basis Omnes problem using physical thresholds.
5. Treat omitted dynamical $m_d-m_u$ corrections to $T_c$ as a systematic uncertainty.
6. Handle Coulomb and real-photon effects separately for charge channels containing two charged particles.

This approximation captures source isospin breaking, the $I=3/2$ wave, physical phase space, and the charged-neutral threshold cusp. It does not constitute a complete next-to-leading-order isospin-breaking scattering calculation.

### Complete fixed-order matching

Match both the ALP production kernels and the $\pi K$ amplitudes to chiral perturbation theory including all terms of the selected order in $m_d-m_u$ and $e^2$. This requires neutral mixing, physical masses in loops, vertex and wave-function corrections, local counterterms, and electromagnetic subtraction in one consistent scheme. The physical charge amplitudes, rather than separately dressed isospin amplitudes, are then the input to the Omnes matrix. [@Kubis:2001ij; @Nehme:2001wf]

### Full three-body treatment

For precision over the entire Dalitz plot, combine the charge-resolved $K\pi$, $K\eta$, and $K\eta'$ channels with crossed-channel rescattering in coupled Khuri--Treiman equations. No consulted source supplies this complete isospin-breaking ALP construction; it is a proposed extension assembled from the cited production, scattering, and form-factor ingredients.

## Verification

Run

```bash
.venv/bin/python scripts/derive_kpi_isospin_breaking.py
```

The script obtains the charge-isospin matrix from SymPy Clebsch--Gordan coefficients, derives the charge-basis scattering matrix and the off-diagonal phase-space term, decomposes the light-quark source into isoscalar and isovector parts, and computes the two physical thresholds. Assertions verify orthogonality, the scattering transformation, and restoration of diagonal isospin phase space when the two thresholds coincide.
