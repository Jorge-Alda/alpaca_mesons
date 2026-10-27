---
title: External one-subtraction KT amplitude grid and software diagnostic
status: numerical
created: 2026-10-27
updated: 2026-10-27
bibtex_keys: [Niecknig:2012sj, Hoferichter:2019mqg, Hoferichter:2025lcz, Pelaez:2019eqa]
provenance: numerical-computation
---

# Inputs, conventions and scope

The adapter [export_vppp_kt_grid.py](../../../scripts/export_vppp_kt_grid.py)
provides the full unit-subtraction scalar required by the
[three-pion width calculator](vppp-three-pion-running-width.md). The physics
prescription and exact-isospin approximation are those of that record
[@Hoferichter:2019mqg; @Hoferichter:2025lcz; @Niecknig:2012sj].

Software provenance: the optional external installation is
[HISKP-ph/khuri_treiman_solver](https://github.com/HISKP-ph/khuri_treiman_solver),
commit `7e32f96455bdccf49d6f532b5d148c039da891e0`, in the ignored directory
`.external/khuri_treiman_solver`. Its GPL source is not copied into this
project. The upstream requested article citation is arXiv:2212.11767;
its INSPIRE record was not retrieved during this diagnostic, so article-level
claims are not inferred from it here. The following interface observations
are from the installed code, not a literature result.

**Unit correction:** the upstream `cpp/example/KT.cpp` explicitly uses
$m_\pi=1$, and the phase continuation in `cpp/src/phase.cpp` matches at
$s=80$. Thus this adapter passes $M/m_\pi$ and $s/m_\pi^2$ to the solver,
and requires the phase table's first column in $m_\pi^2$ units and its second
column in radians. The presence of GeV mass constants elsewhere in the
upstream code does not establish the example's units. The continuation is
hard-coded, so supplying a phase file is not enough to select an arbitrary
high-energy prescription. Its matching requires a nonzero phase derivative;
an exactly constant zero-phase file is unsuitable.

# Prepared physical phase alternatives

The [physical phase-input record](vppp-pipi-phase-input.md) now supplies
three Madrid/GKPY central phase tables from [@Pelaez:2019eqa], with its
2024 erratum included. They cover $s/m_\pi^2=4\ldots80$ in radians and
use the same common pion mass as this adapter. They are documented
alternatives to the exact Bern phase of [@Hoferichter:2019mqg], which has
not been obtained. The solver's fixed continuation above 80 remains an
additional model assumption. Above the common kaon threshold, using only
the phase in a single-channel KT equation also neglects inelastic channels.

The historical calculation below used a synthetic phase and remains a
software diagnostic. Preparing physical phase tables does not turn that
calculation into a physical result. Each selected physical phase requires
its own KT solution and numerical convergence checks.

# Method

1. Import `_khuri_pinocchio` from the explicitly supplied build directory.
2. Instantiate `IterationV3Pi` with `mass_1=1`, `mass_4=sqrt(q_squared)/m_pi`,
   `max_subs=1`, and caller-selected iterations/cutoff. Construction solves
   immediately; do not call `solve_kt` again.
3. Evaluate `SubtractionConstant.a1`, isospin argument 1, `Setting.above`.
   The upstream initial basis is the Omnès function, with unit constant
   subtraction. This is the single-variable basis, not the full scalar.
4. Construct physical Gauss--Legendre quadrature in GeV units, then sum
   $\mathcal F(s)+\mathcal F(t)+\mathcal F(u)$ after converting each invariant
   to $m_\pi^2$ units. No factor $1/3$ is applied. With zero rescattering,
   the full unit-subtraction scalar would be 3, consistently with
   $\mathcal F(s,t,u)=a(q^2)[\mathcal F(s)+\mathcal F(t)+\mathcal F(u)]$.
5. Save the NPZ contract, phase hash, upstream commit and module hash;
   integrate using the existing width script. The common shape is repeated
   along the resonance dimension under the explicitly assumed common phase
   prescription. Distinct vector-specific rescattering inputs require
   separate calculations.

# Synthetic diagnostic, not physical input

No exact physical Bern phase table was available for this historical diagnostic. The numerical check used a
**synthetic**, smooth phase
$\delta(s)=\pi[1-\exp(-(s/m_\pi^2-4)/25)]$,
1921 equally spaced nodes in $s/m_\pi^2\in[4,100]$. This is a software
fixture, not a fitted scattering phase or a substitute for one.

Inputs: $m_\pi=0.13957039$ GeV, $\sqrt{q^2}=0.7,0.78266$ GeV,
one iteration, cutoff 1000 in pion mass units, 8 pair nodes, 6 angular nodes.
The matching parameters are the upstream example values
$(\epsilon,\mathrm{validity},\widetilde\epsilon,
\widetilde{\mathrm{validity}})=(0.1,0.05,0.3,0.2)$ in the solver units.
The pole mass is 0.78266 GeV; the **arbitrary diagnostic** pole width
0.001 GeV tests normalization only. No physical width is claimed.

Results: single-variable basis at zero is $1+6.60820863\times10^{-7}i$
at both virtualities; the spurious imaginary part measures numerical error.
Reflection $z\to-z$ changes the reconstructed full scalar by at most
$9.16\times10^{-16}$. The width program gives 0.00044734 and 0.001 GeV,
with exact pole normalization. These numbers establish the file handoff,
not convergence or agreement with a physical KT solution.

# Reproduction and checks

Prepare a JSON with the fields documented in the script; the diagnostic's
local copy is `build/alpaca_diagnostics/synthetic_configuration.json` under
the external checkout. For reproducibility, its phase is generated with
`numpy.linspace(4,100,1921)` and the explicit formula above.

```bash
.venv/bin/python scripts/export_vppp_kt_grid.py .external/khuri_treiman_solver/build/alpaca_diagnostics/synthetic_configuration.json .external/khuri_treiman_solver/build/alpaca_diagnostics/synthetic_grid.npz
.venv/bin/python scripts/build_vppp_three_pion_width.py .external/khuri_treiman_solver/build/alpaca_diagnostics/synthetic_grid.npz .external/khuri_treiman_solver/build/alpaca_diagnostics/synthetic_width.npz
.venv/bin/python -m unittest discover -s scripts -p export_vppp_kt_grid.py
.venv/bin/python -m unittest discover -s scripts -p build_vppp_three_pion_width.py
```

The new quadrature unit test checks invariant conservation, reflection,
pair integration of unity, and angular integration of $1-z^2$; it passes.
The four width unit tests also pass. The solver imports under Python 3.14.8.
Before physics use, supply the selected physical phase, verify its conversion
and high-energy matching, demonstrate iteration/quadrature/cutoff convergence,
and compare with a published Dalitz distribution or reference amplitude.
One synthetic iteration is not a converged KT solution. The adapter does
not implement uncertainty propagation, radiative corrections, isospin
breaking or production normalization $a(q^2)$.
