---
title: Conditional matching of local VPPP terms to the three-pion subtraction function
status: mixed
created: 2026-10-27
updated: 2026-10-27
bibtex_keys: [GomezDumm:2010ra, Hoferichter:2023bjm, Hoferichter:2025lcz, Wess:1971yu]
provenance: assumption
---

# Matching prescription

This record fixes a usable, explicitly conditional bridge from the local
photon-source vertex to the dispersive normalization $a(q^2)$. For this fit,
choose the overall WZW orientation so the Hoferichter current convention has
positive WZW $F_{3\pi}$, and use the same odd-parity tensor and pion phase
for the local vertex. The local full scalar is assigned to the subtraction seed before
Khuri--Treiman rescattering, with
$$
\Delta a(q^2)=\frac{\mathcal F_{\gamma,\mathrm{local}}(q^2)}{3e}.
$$
The factor $1/3$ follows the no-rescattering matching convention in which a
Dalitz-independent full scalar is the sum of three identical
single-variable functions. After rescattering, this seed is propagated
through the normalized Khuri--Treiman solution and crossed-channel average;
the bare source scalar is not added directly to the physical P wave.

This prescription, including the relative sign choice, is a model choice.
The matching script assumes it; it does not derive it from the gauged WZW
functional. Neither [@Hoferichter:2023bjm] nor
[@Hoferichter:2025lcz] derives an RChT $g_i$-to-$a(q^2)$ relation, and the
repository's [WZW record](wzw-gamma3pi-contribution.md) leaves the overall
functional orientation explicit outside this fit-specific sign choice.

# Subtraction and pole map

Use the Hoferichter subtraction point $q^2=0$. For one bare vector pole,
the photon-source response in the
[source derivation](isoscalar-vppp-photon-response.md) gives
$$
\Delta\alpha_A^{(i)}=
\frac{16\sqrt2\,t_{n,i}R_i m_\pi^2}
{M_Vf_\pi^3M_i^2}
\left[-(g_1-g_2-g_3)+(2g_4+g_5)\right],
$$
and the coefficient of the once-subtracted spectral pole is
$$
\Delta c_i=
\frac{16\sqrt2\,t_{n,i}R_i}{3M_Vf_\pi^3}
\left[(g_1-g_2-g_3)(M_i^2-3m_\pi^2)
+3m_\pi^2(2g_4+g_5)\right].
$$
The exact bare-pole decomposition is
$$
\Delta a_i(q^2)=\Delta\alpha_A^{(i)}
+\frac{q^2\Delta c_i}{M_i^2(M_i^2-q^2)}.
$$
Here $\Delta\alpha_A^{(i)}$ has units $\mathrm{GeV}^{-3}$ and $\Delta c_i$
has units $\mathrm{GeV}^{-1}$, matching the dispersive $a(q^2)$ and pole
residue conventions.

For a timelike fit, the selected continuation is to insert $\Delta c_i$ as
the shift of the corresponding resonance coefficient in the same
once-subtracted spectral term used for the other vector contributions in
[@Hoferichter:2025lcz]. Its zero-width limit is the displayed bare-pole
identity. Do not replace $M_i^2-q^2$ by an independently chosen running-width
denominator; the fit must use the cited dispersive resonance line shape.

# Identifiability and the sum-rule check

For an ideal $\omega$ pole, the map from the two VPPP combinations to
$(\Delta\alpha_A,\Delta c_\omega)$ has rank two at finite pion mass. The
pole residue alone has rank one; the subtraction shift supplies the second
row, and it vanishes in the chiral limit. With nonideal $\omega$--$\phi$
mixing, the two residue shifts alone have rank two for distinct pole masses
and nonzero source overlaps. These are ranks of the local coefficient map,
not of a data likelihood after conformal and resonance nuisances are
profiled.

The local pole contribution does not by itself obey the global sum rule
$$
\alpha_A=\frac{1}{\pi}\int_{s_{\mathrm{thr}}}^{\infty}
\frac{\operatorname{Im}\mathcal A(s')}{s'}\,ds'
+\frac{1}{\pi}\int_{s_{\mathrm{inel}}}^{\infty}
\frac{\operatorname{Im}C_p(s')}{s'}\,ds'
$$
of [@Hoferichter:2025lcz]. For each pole,
$$
\Delta\alpha_A^{(i)}-\frac{\Delta c_i}{M_i^2}
=-\frac{16\sqrt2\,t_{n,i}R_i}{3M_Vf_\pi^3}
(g_1-g_2-g_3).
$$
With only the bare $\omega$--$\phi$ spectral poles, saturation therefore
requires $g_1-g_2-g_3=0$ when the total nonstrange source overlap is nonzero.
For a nonzero value, an inelastic spectral contribution from $C_p$ or a
different high-energy completion must cancel this moment. The matching
script does not invent that completion. The coefficient $g_1-g_2-g_3$
therefore cannot yet be interpreted as independently fit by a pole-only
version of the published representation.

In the selected fit prescription, the published sum rule is retained and the
conformal/inelastic contribution must supply exactly the summed mismatch
shown above. Its moment is constrained by the VPPP combinations; its remaining
energy dependence is fitted with the declared conformal order. If that basis
cannot satisfy the constraint while describing the data, the fit requires a
larger inelastic basis or a separately justified spectral contribution.

The 2025 analysis also quotes the finite-mass estimate
$C_A=1.066(10)F_{3\pi}$. Do not add that full correction as an independent
subtraction constraint while also assigning the mass-dependent VPPP shift
above, unless the overlapping chiral contribution has first been separated.

# Reproducibility and fit use

Run:

    .venv/bin/python scripts/derive_vppp_subtraction_matching.py

The SymPy script prints the subtraction and residue shifts, their exact
once-subtracted identity, ideal and nonideal coefficient ranks, and the
sum-rule moment mismatch. Its focused checks verify the decomposition and
rank limits. For a fit, insert these shifts into the same spectral basis as
the rest of $a(q^2)$, do not add a second free subtraction or duplicate pole
residue for the same local contribution, retain the conformal coefficients
with the published sum rule, and profile their effect on both combinations.
The nonlocal remainder of each total resonance residue must be calculated
or assigned an explicit external constraint before it is added to the fit.
