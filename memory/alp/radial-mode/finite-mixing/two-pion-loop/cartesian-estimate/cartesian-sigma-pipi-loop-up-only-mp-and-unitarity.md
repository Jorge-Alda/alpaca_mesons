---
title: Loop scale m_P in the up-only model and the epsilon at which R reaches one, compared with tree-level unitarity
status: mixed
created: 2026-11-28
updated: 2026-11-28
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh, Murayama:2026ioh, ParticleDataGroup:2024cfk]
provenance: numerical-computation
---

# Scale $m_P$ and the value of $\epsilon$ at $R=1$ in the up-only model

No loop integral, scalar loop function, loop coefficient or loop amplitude is computed. $R$ is a size ratio, not a computed correction. Units: GeV; $\epsilon=f_\pi/f_a=f_\pi\kappa_u/(\sqrt2m_u)$.

## Scale $m_P$

In the up-only model the scale $m_P$ of rule D3 (largest pole mass of the nine meson-like states) is taken to be the mass of the $\eta'$-like pole, $m_P=m_{\eta'}$. This is an input from the user, supported by the check below: the $\eta'$-like pole is the heaviest meson-like pole and its mass changes little with $\kappa_u$.

Check: [script](../../../../../../scripts/check_up_only_eta_prime_heaviest.py), tests (4, pass) in [`test_check_up_only_eta_prime_heaviest.py`](../../../../../../scripts/test_check_up_only_eta_prime_heaviest.py), inputs of [the only-up inputs record](../../../../meson-mixing/finite-mixing/spectrum-and-poles/quadratic-mixing/only-up-finite-mixing-inputs.md), $\epsilon$ from 0.1 to 150 in 25 logarithmic steps.

- The $\eta'$-like pole is the heaviest of the nine meson-like poles at every grid point; the ALP-like pole is excluded and exceeds it for $\epsilon\gtrsim13$.
- $m_{\eta'}$ is not $0.95778$ GeV: the input fixes only the singlet diagonal entry, and $\eta$–$\eta'$ mixing raises the pole to $1.00235$ GeV at small $\epsilon$ (shift $+4.65\%$). Over the grid it stays between $0.9916$ and $1.0037$ GeV (shift $+3.5\%$ to $+4.8\%$); the smallest value is at $\epsilon\simeq13$, where the ALP-like pole crosses it and the overlap label of the solver is conventional.
- For $\epsilon\gtrsim13$, $m_a>m_P$, as assumed in the [up-only estimate record](cartesian-sigma-pipi-loop-up-only-estimate.md).

## $\epsilon$ at $R=1$

From the leading term of the up-only record, $R\to\dfrac{B_0^2m_u^2}{8\pi^2f_\pi^2m_P^2}\epsilon^2$, so
$$
\epsilon_{R=1}=\frac{2\sqrt2\,\pi f_\pi m_{\eta'}}{B_0m_u}\simeq132\quad(m_{\eta'}\simeq1.002\ {\rm GeV}),
$$
with $m_u=2.20$ MeV, $f_\pi=92.3$ MeV, $B_0=m_{\pi^\pm}^2/(m_u+m_d)$ of the inputs record. This is the extrapolation of the large-$\epsilon$ term, and it is outside the range where the estimate is reliable ([limits record](cartesian-sigma-pipi-loop-limits.md)). For the polar estimate, $R\to\frac{21B_0m_u}{64\pi^2f_\pi^2}\epsilon^4$ gives $\epsilon_{R=1}\simeq2.5$ with the same inputs; the two estimates are not equivalent and neither is authoritative ([comparison record](cartesian-sigma-pipi-loop-polar-comparison.md)).

## Comparison with tree-level unitarity

The record [`murayama-alp-unitarity.md`](../../../../radial-mode/couplings-and-decays/murayama-alp-unitarity.md) gives, for the up-only model:

- $\epsilon\le52.5$ (primary, $|{\rm Re}\,a_0|\le1/2$) and $\epsilon\le111$ (secondary, $|a_0|\le1$) from the constant $aP_iP_jP_k$ contacts, with $m_u=2.16$ MeV and $f_\pi=92.07$ MeV (with the inputs of the inputs record above these become $48.7$ and $103.3$);
- $\lambda\le2\pi$ (primary) from $\{aa,\sigma\sigma\}$, independent of $\epsilon$ at large $\epsilon$ because $m_\sigma^2-m_a^2=2\lambda f_a^2$.

## Interpretation (assumption)

- **Unitarity as the bound.** The tree-level partial-wave bounds on $\epsilon$ above are taken as the physical range of $\epsilon$ in the up-only model.
- **$R(\epsilon)$ as the error estimate.** Within that range, $R(\epsilon)=(\epsilon/\epsilon_{R=1})^2$ is taken as the order-of-magnitude relative uncertainty of the tree amplitude of $\sigma\to\pi^+\pi^-$ (and about twice that for a rate, which scales as $|M|^2$). It is a size ratio, not a computed correction, and it is not a bound on the loop.
- **Only the Cartesian estimate is used.** The polar $R\propto\epsilon^4$ is far above one in this range, and since the exact amplitude is the same in both parametrizations the polar $\epsilon^4$ graphs must cancel in the exact loop if both calculations are correct. The polar $R$ is therefore an artefact of the size rule and is not used. The Cartesian $\epsilon^2$ term rests on the single graph $C\sigma9$ at scale $m_P$ and may itself miss cancellations or be an artefact, so $R(\epsilon)$ is a rough scale.
- **Limits.** The unitarity bounds are tree-level, and loops of relative size $R$ shift them by a comparable amount, so the two are consistent only for $R\ll1$. At the bounds, $m_a$ is a few GeV and $f_a$ a few MeV, outside the range where the leading-order chiral Lagrangian and the extrapolation of the large-$\epsilon$ term are reliable ([limits record](cartesian-sigma-pipi-loop-limits.md)).
