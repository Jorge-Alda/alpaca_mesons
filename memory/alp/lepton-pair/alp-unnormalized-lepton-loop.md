---
title: Unnormalized ALP dilepton loop and pointlike ultraviolet subtraction
status: derived
created: 2026-10-25
updated: 2026-10-25
bibtex_keys: [Messerli:2025rnv, Alda:2026nrd]
provenance: analytical-derivation
---

# Inputs and conventions

Use the scalar $F_{a\gamma^*\gamma^*}$ in the [tree TFF record](alp-doubly-virtual-tff.md), with $i\mathcal M=i e^2R F$, metric $(+---)$, $\epsilon^{0123}=+1$, and no-mixing sources in Alda:2026nrd's convention [@Alda:2026nrd]. The loop integrand is the [Messerli:2025rnv functional](../../extracts/Messerli:2025rnv/dilepton-loop-formalism.md); its bibliography entry was inspected [@Messerli:2025rnv]. The [normalization record](alp-lepton-pair-normalization.md) supplies the independently derived spin sum and two-real-photon cut. Natural units are used; masses are in GeV and $F$ in $\mathrm{GeV}^{-1}$.

Let $q=p_++p_-$, $q^2=s=m_a^2$, and $p=p_-$ in the lepton propagator. All propagators carry the Feynman $+i0$ prescription. Work at leading QED order and $s>4m_\ell^2$. The user specifies $2m_e<m_a<2m_\tau$ and full loop virtualities. No numerical TFF, local coefficient, matching scale, or subtraction constant is selected.

# Removing the real-photon normalization

For a nonzero real-photon value, the existing reduced amplitude obeys the loop-coefficient convention $\left.g_\ell\right|_{\mathrm{TFF}}=2\alpha^2m_\ell F_{a\gamma\gamma}\mathcal A_\ell$. The executable [derivation](../../../scripts/derive_alp_unnormalized_lepton_loop.py) substitutes the normalized TFF into that relation and cancels $F_{a\gamma\gamma}$ before integration. The resulting definition is
$$
\left.g_\ell\right|_{\mathrm{TFF}}=
\frac{4i\alpha^2m_\ell}{\pi^2s}\int d^4k\,
\frac{sk^2-(q\cdot k)^2}
{[k^2+i0][(q-k)^2+i0][(p-k)^2-m_\ell^2+i0]}
F_{a\gamma^*\gamma^*}(k^2,(q-k)^2).
$$
This unnormalized integral also defines the coefficient when $F_{a\gamma\gamma}=0$. Its dimension is zero, as required for the coefficient of $\bar u i\gamma_5v$. The equality fixes the coefficient convention used here; it does not independently match a UV local ALP--lepton Lagrangian or its relative sign to that convention.

SymPy expands the integrand for two independent source responses and obtains exactly the corresponding linear combination of integrands. Integration preserves this linearity when the same regularization, subtraction, and matching are used for each response. Two responses with equal values at $(0,0)$ can enter with opposite coefficients, so their real-photon values cancel while their off-shell difference remains nonzero. The script imposes that equality only at the origin and leaves the two functions independent elsewhere. There is no division by the cancelling value and no conclusion that the loop vanishes.

# Coherent local contribution and width

Express any separately matched local contribution in the same coefficient convention. With $\mathcal M=(\left.g_\ell\right|_{\mathrm{TFF}}+\left.g_\ell\right|_{\mathrm{local}})\bar u i\gamma_5v$, the executed spin sum gives
$$
\Gamma(a\to\ell^+\ell^-)=
\frac{\sqrt{s}\,\sigma_\ell}{8\pi}
\left|\left.g_\ell\right|_{\mathrm{TFF}}+
\left.g_\ell\right|_{\mathrm{local}}\right|^2,
\qquad \sigma_\ell=\sqrt{1-\frac{4m_\ell^2}{s}}.
$$
The modulus contains interference. The local term is not determined by the hadronic TFF or by this width identity. Its scale, scheme, and conversion from a UV Wilson coefficient remain matching inputs. A consistently matched local contribution can absorb ultraviolet dependence assigned to a pointlike photon term; it cannot be added independently to an already matched contribution a second time.

# Two-real-photon discontinuity

The script reintegrates the existing cut angular function $-\pi/[2(1-\sigma_\ell z)]$ over $-1<z<1$. Multiplication by the unnormalized coefficient gives
$$
\frac{1}{2i}\operatorname{Disc}_{\gamma\gamma}
\left.g_\ell\right|_{\mathrm{TFF}}
=\frac{\pi\alpha^2m_\ell F_{a\gamma\gamma}}{\sigma_\ell}
\log\!\left(\frac{1-\sigma_\ell}{1+\sigma_\ell}\right).
$$
For real coefficients and compatible Schwarz reflection this is the corresponding imaginary part. For general complex source coefficients retain the discontinuity statement. The reduced cut tends to $-\pi$ at threshold, while the unnormalized coefficient remains finite and the two-body width vanishes with $\sigma_\ell$. Its massless-lepton behavior is proportional to $m_\ell\log(m_\ell^2/s)$ and vanishes, verifying helicity suppression. If source contributions cancel $F_{a\gamma\gamma}$, this particular cut vanishes; the off-shell dispersive loop and other cuts need not vanish. Additional hadronic cuts prevent interpreting this contribution alone as the full imaginary part or a universal lower bound.

# Pointlike ultraviolet behavior

This subsection is a declared approximation used to diagnose the subtraction, not a physical full-range ALP TFF. Replace $F$ by a constant $F_{a\gamma\gamma}$ and take the hard loop momentum while the external momenta and masses remain fixed. SymPy performs the Lorentz-symmetric angular replacement $(q\cdot k)^2\to sk^2/d$ and obtains $(1-1/d)sk^2$ for the numerator. At $d=4$, the leading Minkowski integrand for the coefficient is $3i\alpha^2m_\ell F_{a\gamma\gamma}/(\pi^2(k^2)^2)$.

Under $k^0=iK^0$, $d^4k=i\,d^4K$ and $k^2=-K^2$, the Euclidean angular area is $2\pi^2$. The executable radial integral between positive scales $\mu$ and $\Lambda$ is therefore
$$
\left.\Delta g_\ell\right|_{\mathrm{TFF,hard}}
=-6\alpha^2m_\ell F_{a\gamma\gamma}
\log\frac{\Lambda}{\mu}
=-3\alpha^2m_\ell F_{a\gamma\gamma}
\log\frac{\Lambda^2}{\mu^2}.
$$
This is the leading logarithm in the stated Euclidean radial-cutoff and loop-coefficient convention, with fixed $\alpha$ and constant $F$. It is not a finite pointlike prediction or an identification of a minimal-subtraction constant. Cutoff independence requires the opposite local dependence,
$$
\frac{d}{d\log\Lambda}\left.g_\ell\right|_{\mathrm{local}}
=+6\alpha^2m_\ell F_{a\gamma\gamma},
$$
for this assigned pointlike contribution. The finite local part is not fixed. A hadronic component falling as $K^{-n}$ with $n>0$ instead has a convergent hard radial integral $\int_\mu^\infty dK\,K^{-1-n}=\mu^{-n}/n$. Establishing that falloff, its source-specific coefficient and finite-mass corrections requires the independent high-energy matching in the [isoscalar plan](alp-isoscalar-tff-plan.md); a retained stable-pole expression cannot establish the needed QCD continuation by itself.

# Reproducibility and remaining work

Executed:

```text
.venv/bin/python scripts/derive_alp_unnormalized_lepton_loop.py
.venv/bin/python scripts/derive_alp_unnormalized_lepton_loop.py --test
```

Four unit tests pass. They check removal of the real-photon denominator and source linearity, the integrated cut and threshold limit, the coherent width and helicity limit, and the ultraviolet logarithm/local cancellation plus convergent-tail integral. The script derives expressions with SymPy and prints intermediate and final results with its LaTeX printer.

Stage 6 is complete only for this unnormalized functional, cut normalization, and ultraviolet bookkeeping. A numerical full-loop result over the user's mass range still requires the complete source-resolved TFF with all represented cuts and analytic continuation, a controlled high-energy part, and a matched local contribution in one scale and scheme. No numerical integration or full ALP width has been produced.
