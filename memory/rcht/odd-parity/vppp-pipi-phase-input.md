---
title: Physical pion-scattering phase tables for the external KT calculation
status: numerical
created: 2026-10-27
updated: 2026-10-27
bibtex_keys: [Pelaez:2019eqa, ParticleDataGroup:2024cfk, Hoferichter:2019mqg]
provenance: numerical-computation
---

# Conventions, inputs and assumptions

The [literature extract](../../extracts/Pelaez:2019eqa/p-wave-phase-parameterization.md)
supplies the full $I=1$, $P$-wave scattering amplitude, including the small
inelastic contribution, and all corrected central parameters from
[@Pelaez:2019eqa] and its 2024 erratum. The three solutions are published
model alternatives. Their spread is not a statistical confidence interval.

Use exact-isospin kinematics with $m_\pi=0.13957039$ GeV, matching the
[KT adapter](vppp-external-kt-grid.md). **Mass approximation:** all pion
states have the charged-pion mass. The one kaon threshold uses
$m_K=(0.493677+0.497611)/2=0.495644$ GeV, the arithmetic mean of charged
and neutral kaon masses from [@ParticleDataGroup:2024cfk]. We do not combine
charged and neutral kaon thresholds within the amplitude. This is an
explicit exact-isospin input convention, not a refit of the article's
parameters or a determination of isospin breaking.

The output has dimensionless first column $s/m_\pi^2$ and second column
$\delta^1_1$ in radians. The domain $[4,80]$ ends at approximately
$\sqrt{s}=1.248$ GeV, below the parameterization's 1.4 GeV matching energy.
Coefficients controlling the article's separate region above 1.4 GeV
are therefore not needed. Table marginal errors are stored as metadata;
no correlated uncertainty is invented or propagated.

# Reproducible method

The deterministic script
[prepare_vppp_pipi_phases.py](../../../scripts/prepare_vppp_pipi_phases.py)
performs these steps:

1. Set the documented masses and corrected coefficients for I, II, III.
2. Evaluate the conformal amplitude and the physical upper-lip kaon-loop
   inelastic amplitude on a dense real grid.
3. Assemble the full amplitude from Eq. (20), including its interference
   term. Form $S=1+2i\sigma_\pi t^1_1$ and compute
   $\delta^1_1=\operatorname{unwrap}(\arg S)/2$, starting from zero at
   threshold. The threshold endpoint uses the limiting value $S=1$.
4. Write separate header-free two-column tables and a JSON record of the
   exact input coefficients, masses, source references and marginal errors.
5. Check threshold behavior, elastic unitarity below the kaon threshold,
   continuity, monotonic phase over the supplied domain, phase range,
   equivalence with the factorized $S$ matrix, and round-trip readability.

Invocation:

```bash
.venv/bin/python scripts/prepare_vppp_pipi_phases.py
.venv/bin/python -m unittest discover -s scripts -p prepare_vppp_pipi_phases.py
```

# Outputs and scope

- [Solution I](../../../data/odd_parity/pipi_pwave_pelaez_2019_I.dat).
- [Solution II](../../../data/odd_parity/pipi_pwave_pelaez_2019_II.dat).
- [Solution III](../../../data/odd_parity/pipi_pwave_pelaez_2019_III.dat).

These are **Madrid/GKPY alternatives**, not the exact Bern phase input
of [@Hoferichter:2019mqg]. Above the kaon threshold the physical scattering
amplitude is inelastic; the single-channel KT solver retains its phase
only. This approximation must accompany any ensuing three-pion result.
The solver further fixes a continuation toward $\pi$ above $s/m_\pi^2=80$;
that continuation is an additional model assumption, not the published
Madrid parameterization above this point. Different phase inputs require
new KT amplitude grids, not substitution into previously solved grids.
Iteration, quadrature and continuation sensitivity of the physical KT
solution remain separate from the table preparation checks.

The [input manifest](../../../data/odd_parity/pipi_pwave_pelaez_2019_manifest.json)
records all coefficients, marginal errors and SHA-256 hashes. Each table
contains 3801 numeric rows with spacing 0.02 in $s/m_\pi^2$ and no header,
as required by the upstream stream-based reader.

## Verification result

The three focused unit-test methods pass. All solutions start at zero,
are continuous across the common kaon threshold and increase monotonically
on the supplied grid. The endpoint diagnostics are:

| Solution | $\delta^1_1(80)$ [rad] | $d\delta^1_1/d(s/m_\pi^2)$ at 80 |
|---|---|---|
| I | 2.9218643241 | 0.004330648987 |
| II | 2.9150444091 | 0.004063830378 |
| III | 2.9209761474 | 0.004394348257 |

The endpoint derivative is estimated by a three-point backward difference;
it is nonzero, as required for the installed continuation's matching.
These checks establish phase-table preparation, not convergence or physical
validation of a subsequent KT solution. No GPL solver source was modified.
