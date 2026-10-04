---
title: Why the BaBar ISR three-pion table cannot constrain the one-multiplet couplings
status: derived
created: 2026-10-04
updated: 2026-10-04
bibtex_keys: [BaBar:2004ytv, Niecknig:2012sj]
provenance: assumption
---

# Why the BaBar ISR $3\pi$ Table Is Not a First-Release $d_i$ Constraint

## Decision

Do **not** build an $e^+e^-\to\pi^+\pi^-\pi^0$ Born cross-section predictor for the
current first release. The tabulated range does not contain the resonances where
the VVP sensitivity lives, and the region it does cover is one our
single-multiplet model is documented to fail. Recorded so the table is not
re-litigated as "the obvious first target".

## Three independent reasons

**1. Neither peak is in the data, and no restricted window helps.** The 78 rows of
[`babar_3pi_born_spectrum.csv`](../../../data/odd_parity/babar_3pi_born_spectrum.csv)
are bin centres on a uniform 25 MeV grid spanning $\sqrt{s'}=1.0625$ to
$2.9875$ GeV, i.e. 1.05--3.00 GeV. The $\omega$ at 782.7 MeV and the $\phi$ at
1019.5 MeV both lie **below** the range. The measured cross section peaks at
$\sqrt{s'}=1.3125$ GeV, in the $\omega'(1350)$ region. The $d_i$ sensitivity
lives in $\omega\to3\pi$ and $\phi\to3\pi$, whose peaks are absent.

Restricting to a low-mass window does not rescue it: $\omega'$ sits only
$0.73\,\Gamma_{\omega'}$ above the $\phi$, its Lorentzian floor is 11% of its
peak value, and $\omega'+\omega''$ own at least **72.6%** of the coherent sum at
every one of the 78 bins. The quantitative derivation, including the error-budget
trap of a benign-looking but nuisance-dominated systematic, is in the
[separability record](three-pion-resonance-separability.md).

**2. Our model cannot describe this window.** BaBar's own fit needs four
isoscalar isoscalar resonances, $\omega,\phi,\omega',\omega''$, with
$\mathcal B(V\to e^+e^-)\mathcal B(V\to3\pi)$ as *fit parameters*, a $\rho\pi$
intermediate state, and a fixed relative phase $\phi_\phi=163^\circ$; the paper
claims a good description only **up to 1.8 GeV** [@BaBar:2004ytv]. Our baseline
is a single multiplet ($\rho,\omega,\phi$) with no $\omega'$ or $\omega''$, so
this is precisely the range where the
[one-multiplet assumption](odd-parity-rcht-implementation-decisions.md) is
weakest. See the [line-shape extract](../../extracts/BaBar:2004ytv/mass-binning-and-fitted-products.md).

**3. The measured line shape contains no direct VVP term.** BaBar assumes
$V\to3\pi$ proceeds *exclusively* through $\rho\pi$; the fitted products absorb
the total. Their remark that VMD gives $\Gamma(\omega\to3\pi)\approx1.4$ MeV
against the measured $7.52$ MeV is Introduction motivation, not part of the fit
model [@BaBar:2004ytv]. Their analysis therefore constrains a product of a
$\rho\pi$-mediated amplitude and fitted propagators, not the VVP operators.

**4. Fixed-width isobars are the wrong tool here.** Crossed-channel three-body
cuts shift and attenuate the $\rho$ band in $\phi\to3\pi$ and shift and enhance
it in $\omega\to3\pi$, and the same source calls fixed-width isobar models
"doomed to fail" for these two decays [@Niecknig:2012sj]. A tree-level
single-multiplet amplitude mispredicts exactly the shape measured. The remedy is
scoped separately in [the dispersive treatment plan](dispersive-three-pion-treatment-plan.md).

## What it would still require

A predictor would have to (a) restrict to the low bins near 1.05--1.2 GeV, (b)
carry $\omega'$ and $\omega''$ as **fixed** BaBar inputs rather than fitted
states, and (c) dress the vector coupling with $\Pi(m_V^2)^{-1}$, because BaBar
absorb vacuum polarization into the resonance electronic width rather than
applying it to the cross section. Point (c) is the one genuinely reusable piece
of machinery: at the $\phi$ the hadronic $\Pi$ is $O(10^3)$, so omitting the
dressing is not a small correction. See
[the Born-definition extract](../../extracts/BaBar:2004ytv/born-cross-section-definition.md).

## Consequences

The table stays a checksum-validated loader with no `Observation`, exactly as
recorded in [the data-implementation record](odd-parity-data-implementation.md).
Its additive `sigma_syst_nb` nuisance and the $1.2\%$ luminosity nuisance remain
declared and unused. No predictor was written and no fit was run against it.