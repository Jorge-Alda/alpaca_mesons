---
title: Scoped plan for a dispersive three-pion final-state treatment
status: mixed
created: 2026-10-04
updated: 2026-10-04
bibtex_keys: [Niecknig:2012sj, GomezDumm:2013sib, KLOE:2003kas, BaBar:2004ytv]
provenance: assumption
---

# Dispersive $V\to3\pi$ Treatment: Scoped Project and Staged Alternative

## Why this is on the table

Three separate blockers now point at the same root cause. The multi-body final
state is the problem:

- [BaBar $3\pi$](babar-3pi-spectrum-fit-blocker.md) is unusable because
  $\omega',\omega''$ own the whole range, with no clean window
  ([separability record](three-pion-resonance-separability.md)).
- [KLOE $\phi\to3\pi$](kloe-phi-three-pion-fit-blocker.md) is unusable because
  the held data is a four-number amplitude decomposition that a tree-level
  amplitude mispredicts: crossed-channel three-body cuts *shift and attenuate*
  the $\rho$ band in $\phi\to3\pi$ and *shift and enhance* it in $\omega\to3\pi$
  [@Niecknig:2012sj].

Our baseline is a single multiplet at tree level, which is the regime the
literature describes as inadequate for these decays.

## (c) — the full dispersive treatment

Implement $V\to3\pi$ in the Omnès framework [@Niecknig:2012sj]:

$$
\mathcal M(s,t,u)=i\epsilon_{\mu\nu\alpha\beta}n^\mu p_+^\nu p_-^\alpha p_0^\beta\,
\mathcal F(s,t,u),\qquad
\mathcal F(s,t,u)=\mathcal F(s)+\mathcal F(t)+\mathcal F(u),
$$

$$
\mathcal F(s)=\Omega(s)\left\{a+\frac{s}{\pi}\int_{4m_\pi^2}^\infty
\frac{ds'}{s'}\frac{\sin\delta(s')\hat{\mathcal F}(s')}{|\Omega(s')|(s'-s)}\right\},
\qquad
\Omega(s)=\exp\!\left\{\frac{s}{\pi}\int_{4m_\pi^2}^\infty
\frac{\delta(s')}{s'(s'-s)}ds'\right\}.
$$

| Component | Content | Size |
|---|---|---|
| Omnès $\Omega(s)$ | **already implemented** — `scripts/implementation/omnes.py`, `PelaezRodas2016ElasticPhase` [@Pelaez:2016tgi] | done |
| $\hat{\mathcal F}(s)$ inhomogeneity | $\pi\pi$ P- **and** S-wave crossed-channel angular averages, including $\rho\pi$ and $\pi\pi\pi$ crossed topologies | 1--2 weeks |
| Unitarity integral equation | solve on an $s$-grid for a complex function; iterative relaxation | 3--5 days |
| Dalitz density | 2-D Dalitz integration, normalized to $\Gamma_{V\to3\pi}$ | 1 week |
| WZW/anomaly in the same formalism | anomalous inhomogeneity and $L_\mu^{KZ}$ handling, so the subtraction constant $a$ is the RChT quantity we want to fit | 1--2 weeks |
| Validation | Watson final-state-interaction test, $\omega$ enhancement versus $\phi$ attenuation, reproduce published Dalitz plots and $\|\alpha\|,\|\beta\|$ | 1--2 weeks |

**Estimate: 6--10 weeks for one channel**, and it does not stop there. The
$\omega$ and $\omega''$ channels and the radiative set all need it for internal
consistency before any coupling extraction is meaningful. Six months to a year
for a publishable treatment. This is a project in its own right, larger than the
coupling fit it would serve.

### The part that keeps expanding

The integral is not the hard bit. The hard bit is that the subtraction constant
$a$ — the object carrying the $d_i$ information — is defined *relative to* a
specific $\hat{\mathcal F}$ and a specific $\pi\pi$ amplitude with specific
isospin conventions. Mapping $a$ onto RChT inherits all of it: chiral loops in
the $\rho\pi$ propagator to avoid the $\Im A_P$ double count that
[@GomezDumm:2013sib] warn violates analyticity at $O(p^6)$; $\eta$-$\pi$ mixing
in the crossed channels; and unitarity of the *coupled* $\rho\pi/\omega\pi$
system. That is a second amplitude-modelling project stacked on the first.

## (c') — the staged alternative

Capture most of the physics for days rather than months:

1. Replace the fixed-width intermediate propagator by the **Omnès function**
   $\Omega_{11}(s)$ in the tree-level RChT amplitude. This is already available.
2. Omit the crossed-channel inhomogeneity $\hat{\mathcal F}$ entirely.
3. Carry the residual as an explicit systematic, bounded by the published
   numbers: once-subtracted with crossed channels reaches
   $\chi^2/\mathrm{ndof}\approx1.2$, versus $\approx2.0$ with the Omnès function
   alone [@Niecknig:2012sj].

(c′) does **not** give a defensible Dalitz *shape*. It does give a defensible
$\Gamma_{V\to3\pi}$ constraint and an honest price tag for the shape. Note that
the total width is far more robust than the density precisely because a
shift-and-attenuation of the $\rho$ band perturbs an integral much less than a
pointwise comparison.

## Operator-basis caveat: (c) may not be sufficient even if built

Before any of the effort above, one must confront a limitation of the coupling
basis itself. From the
[VVP extract](../../extracts/Ruiz-Femenia:2003jdx/vvp-three-pseudoscalar.md),
Ruiz--Femenia et al. state that their four VVP operators are

> a complete basis for constructing vertices with only **one-pseudoscalar**; for
> a larger number of pseudoscalars additional operators may emerge

[@Ruiz-Femenia:2003jdx]. Consequences:

- **The direct $V\to3\pi$ contact operators are not written down.** The paper
  argues a direct term is necessary — VMD gives $\Gamma(\omega\to3\pi)\approx1.4$
  MeV against the measured $7.52$ MeV [@BaBar:2004ytv] — but does not supply the
  operators, because they lie outside its basis.
- **The published $\rho$-exchange amplitude is entirely OPE-determined.** It
  depends only on $d_1+8d_2-d_3$ and $d_3$ (Eq. (30)), and the short-distance
  matching fixes $d_1+8d_2$ and $d_3$ *individually*. Hence $d_1+8d_2-d_3$ is
  fixed too, and the whole $\rho$-exchange three-pion amplitude carries **zero
  free-$d_i$ information** beyond the OPE prior.
- **$d_4$ enters no three-pion observable**: not the $\rho$-exchange amplitude,
  and not the $\langle VVP\rangle$ matching.

So the three-pion route carries no free-$d_i$ information unless new
three-pseudoscalar VVP operators are constructed and matched — an EFT-derivation
project in its own right, upstream of (c). This must be settled **before**
funding (c), because a dispersive treatment of an OPE-fixed amplitude would
reproduce the prior rather than test it.

By contrast, $V\to P\gamma$ is a **one-pseudoscalar** vertex, which is exactly
the regime where the published basis is stated to be complete. That is a
structural argument for the radiative route, independent of the saturation
question in [the probe record](radiative-channel-saturation-probe.md).

## Recommendation

**Defer (c).** Run the cheap
[radiative saturation probe](radiative-channel-saturation-probe.md) first: it
costs days and decides whether the radiative route works at all. If the
radiative datums turn out to be saturated too, then (c) stops being optional and
becomes the only route — at which point it should be scoped and staffed as its
own project with its own timeline, not appended to this one.

## Explicitly not claimed

No dispersive predictor exists in this repository. The Omnès class is an
operator-level kernel used elsewhere in the project; it has not been wired to
any $V\to3\pi$ amplitude, no Dalitz integration has been written, and no fit has
been run. All sizes above are estimates, not measurements.