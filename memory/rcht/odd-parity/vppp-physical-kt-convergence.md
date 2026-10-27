---
title: Physical-input KT numerical convergence checks
status: numerical
created: 2026-10-27
updated: 2026-10-27
bibtex_keys: [Pelaez:2019eqa, Hoferichter:2019mqg, Hoferichter:2025lcz, Colangelo:2018mtw]
provenance: numerical-computation
---

# Inputs, conventions and assumptions

The [physical pion phases](vppp-pipi-phase-input.md) provide corrected
Madrid solutions I, II and III from [@Pelaez:2019eqa]. The
[external adapter](vppp-external-kt-grid.md) documents the installed solver
revision and the unit-subtraction convention. The
[three-pion width record](vppp-three-pion-running-width.md) supplies the
phase-space measure and normalized width relation.

All pions have mass $m_\pi=0.13957039$ GeV. The full dimensionless scalar
is the sum of three single-variable unit-$a_1$ solutions; no production
subtraction function or resonance poles are included. The two comparison
energies $M_\omega=0.78266$ GeV and $M_\phi=1.019461$ GeV are explicit
numerical reference points, not new fitted mass determinations. No
on-pole widths enter this calculation: normalized width shapes are computed
directly as $(q^2/M^2)^{3/2}\sigma_\mathrm{unit}(q^2)/
\sigma_\mathrm{unit}(M^2)$.

The solver matches the supplied phase to its fixed continuation above
$s/m_\pi^2=80$. This remains a model assumption, as does using only the
phase in the single-channel equation above the kaon threshold. These runs
are Madrid alternatives and do not reproduce the Bern/Roy phase input of
[@Hoferichter:2019mqg]. Numerical convergence cannot establish the adequacy
of either approximation or agreement with measured Dalitz distributions.

# Reproducible method

The deterministic [checker](../../../scripts/check_vppp_kt_convergence.py)
performs the following operations:

1. Solve each selected phase at $\sqrt{q^2}=0.7,0.78266,0.9,1.019461,1.1$ GeV.
2. For solution I vary the dispersive upper integration limit between
   $\sqrt{s_\mathrm{max}}=2,2.5,3$ GeV, corresponding to dimensionless
   limits $205.3399289,320.8436389,462.0148400$. This conversion is essential:
   the solver takes $s_\mathrm{max}/m_\pi^2$, not an energy in GeV.
3. Construct six iterations once per physical point, then evaluate each
   stored iteration on the same $24\times16$ pair-invariant/angular
   Gauss--Legendre quadrature. The matching parameters remain the upstream
   example values $(0.1,0.05,0.3,0.2)$ in pion units.
4. At the 2.5 GeV solution-I reference, vary the external phase-space
   quadrature through $12\times8$, $24\times16$, $48\times32$, $96\times64$.
   The internal solver integration rules and interpolation grids remain
   fixed by upstream; this is not a variation of those numerical settings.
5. Repeat solution II and III at 2.5 GeV with six iterations.
6. Separately solve solution I through eight iterations at the two pole
   reference energies and compare the full complex scalar point by point
   on a $48\times32$ quadrature. Report maximum pointwise relative changes
   and Euclidean norm changes, including complex phases.
7. Check finiteness, unit subtraction at zero, $t\leftrightarrow u$
   reflection symmetry, and normalized width shapes.

```bash
.venv/bin/python scripts/check_vppp_kt_convergence.py .external/khuri_treiman_solver/build/alpaca_diagnostics/convergence
.venv/bin/python scripts/check_vppp_kt_convergence.py .external/khuri_treiman_solver/build/alpaca_diagnostics/convergence --mode amplitude
.venv/bin/python -m unittest discover -s scripts -p check_vppp_kt_convergence.py
```

Large solver outputs, amplitude arrays and detailed JSON metrics remain
under the ignored external build directory. `summary.json` contains the
integral comparisons and normalized width shapes;
`amplitude_iterations.json` contains the pointwise iteration comparisons.
The summary unit test checks an exactly known power-law integral, both
pole normalizations and comparison orientation. The adapter quadrature
unit test and four width unit tests also pass.

# Numerical results

All 25 phase/energy/cutoff solves completed without nonfinite amplitudes.
Across the scan the largest unit-subtraction error
$|\mathcal F(0)-1|$ is $3.47\times10^{-7}$; the largest reconstructed
$t\leftrightarrow u$ reflection difference divided by the maximum scalar
magnitude is $2.98\times10^{-16}$. The latter checks the symmetric
reconstruction, not an independent physical crossing relation.

The following entries are fractional changes in the integrated
$\sigma_\mathrm{unit}/\alpha^2$, quoted in percent. The reference is
solution I, 2.5 GeV cutoff, six iterations, $24\times16$ quadrature.
For the quadrature row the denominator is the corresponding
$96\times64$ result. For the iteration row it is the corresponding
six-iteration result.

| Comparison | $\omega$ reference | $\phi$ reference |
|---|---:|---:|
| One iteration / six iterations minus one | $+3.73239\%$ | $+9.86596\%$ |
| Five iterations / six iterations minus one | $+0.000344\%$ | $-0.069848\%$ |
| $24\times16$ / $96\times64$ minus one | $-0.0000181\%$ | $-0.0000886\%$ |
| $48\times32$ / $96\times64$ minus one | $-0.00000444\%$ | $-0.00000120\%$ |
| 2 GeV / 2.5 GeV cutoff minus one | $-0.185094\%$ | $-0.317436\%$ |
| 3 GeV / 2.5 GeV cutoff minus one | $+0.061382\%$ | $+0.121266\%$ |
| Phase II / phase I minus one | $-0.001838\%$ | $+0.358578\%$ |
| Phase III / phase I minus one | $-0.135091\%$ | $+0.626264\%$ |

The eight-iteration complex-amplitude comparison at $48\times32$ gives:

| Maximum pointwise relative scalar change | $\omega$ reference | $\phi$ reference |
|---|---:|---:|
| Six relative to eight iterations | $4.42\times10^{-7}$ | $1.22\times10^{-4}$ |
| Seven relative to eight iterations | $4.87\times10^{-8}$ | $3.57\times10^{-5}$ |

Over all five reference energies the $48\times32$ integral differs from
$96\times64$ by at most $4.44\times10^{-8}$ relative. The baseline
five-to-six iteration difference is at most $6.98\times10^{-4}$ relative.
Thus external phase-space integration is much more accurate than the
observed iteration, cutoff and phase dependences on these points. These
finite comparisons do not bound the error of the fixed internal solver
grids or an extrapolation to infinite dispersive cutoff.

## Normalized width shapes

For solution I at 2.5 GeV with six iterations and $24\times16$ quadrature,
the following ratios have exact pole normalization by construction. They
are conditional shapes, not widths in GeV or branching fractions.

| $\sqrt{q^2}$ [GeV] | $\Gamma_{3\pi}(q^2)/\Gamma_{3\pi}(M_\omega^2)$ | $\Gamma_{3\pi}(q^2)/\Gamma_{3\pi}(M_\phi^2)$ |
|---|---:|---:|
| 0.7 | 0.2280147 | 0.01088233 |
| 0.78266 | 1 | 0.04772642 |
| 0.9 | 5.716964 | 0.2728502 |
| 1.019461 | 20.95276 | 1 |
| 1.1 | 38.75463 | 1.849619 |

Over these points, changing the cutoff from 2.5 GeV to 2 or 3 GeV changes
the $\omega$-normalized shape by at most $0.1921\%$ or $0.0686\%$,
respectively; for the $\phi$-normalized shape the corresponding maxima
are $0.1869\%$ and $0.0664\%$. At the same cutoff, phase II changes these
shapes by at most $0.5143\%$ and $0.3644\%$, while phase III changes them
by at most $0.9435\%$ and $0.8031\%$. These are central-model
sensitivities, not confidence intervals and not a comparison with the
unavailable exact Bern phase.

# Consequences and remaining scope

A single KT iteration was insufficient for percent-level shape work.
The demonstrated iteration and external quadrature behavior supports
using eight iterations and $48\times32$ quadrature for the next grids,
retaining cutoff and phase alternatives as explicit variations. The
present calculation covers five virtualities and cannot replace a grid
or interpolation-convergence study across every experimental bin. It also
does not check the upstream internal grids, match a published Dalitz
distribution, determine absolute vector partial widths or establish
production-subtraction matching. None of those issues makes the exact
Bern phase a prerequisite for a fit using a declared alternative phase
model; it is needed to reproduce that specific literature baseline.


# Direct comparison of phase matching points

The checker also linearly interpolates each prepared phase table at
$\sqrt{s}=0.8,1.15$ GeV. Values below are in degrees; the solver itself
always receives radians. The comparison input is documented in the
[Roy/VFF source extract](../../extracts/Colangelo:2018mtw/roy-p-wave-phase-input.md)
[@Colangelo:2018mtw].

| Phase prescription | $0.8$ GeV | $1.15$ GeV |
|---|---:|---:|
| Prepared Madrid solution I | 107.3943 | 163.8513 |
| Prepared Madrid solution II | 107.4265 | 163.6314 |
| Prepared Madrid solution III | 107.5341 | 163.7817 |
| Earlier Roy input estimates | $108.9(2.0)$ | $166.5(2.0)$ |
| Roy/VFF fitted central values | $110.4(0.7)$ | $165.7(2.4)$ |

Solution I therefore differs from the Roy/VFF fitted central values by
$-3.0057^\circ$ and $-1.8487^\circ$ at these points. The earlier Roy
estimates and the fitted phase points must not be conflated. In particular,
the Madrid paper's compatibility statement within CFD uncertainties is
not a bound on point-by-point differences from the final fitted Roy/VFF
central curve. The small I/II/III spread seen here cannot quantify the
Madrid-versus-Roy/VFF variation in the fit. That comparison requires the
full alternate phase curve, not only its two matching values.
