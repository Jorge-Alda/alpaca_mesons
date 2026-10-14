---
title: Pion single- and double-virtual short-distance inputs
status: literature
created: 2026-10-13
updated: 2026-10-15
bibtex_keys: [Knecht:2001qf]
provenance: literature
---

# Pion TFF targets

Source `pionpole.tex`, lines 742–747, states that the single-virtual
form factor falls as $1/Q^2$ and that the LMD ansatz fails this test,
whereas LMD+V can satisfy it with $h_1=0$ [@Knecht:2001qf].

Appendix A, equation labelled `OPE_FF`, gives

$$
\mathcal F(\lambda^2q^2,(p-\lambda q)^2)
=\frac{2F_\pi}{3q^2}\left[
\lambda^{-2}+\lambda^{-3}\frac{q\cdot p}{q^2}
+O(\lambda^{-4})\right],
$$

up to $O(\alpha_s)$ corrections (source lines 1847–1862)
[@Knecht:2001qf]. The anomaly phase in this paper is negative,
$\mathcal F(0,0)=-N_C/(12\pi^2F_\pi)$. The repository now uses
$\epsilon^{0123}=+1$, $\epsilon_{0123}=-1$ in the metric $(+---)$,
retaining its external-meson phase. Its scalar anomaly consequently has
this negative sign, and the corresponding symmetric Euclidean target is
$\lim_{Q^2\to\infty}Q^2F_\pi(-Q^2,-Q^2)=-2F_\pi/3$.
The formerly documented positive target belonged to the opposite common
tensor orientation. This scalar translation preserves the physical tensor
amplitude and relative resonance/anomaly signs; see the
[canonical convention record](../../rcht/odd-parity/levi-civita-convention.md)
and its [executable evidence](../../../data/vector_sector/levi_civita_convention.md).

The paper permits replacing the chiral-limit decay constant by $F_\pi$
for these stated leading pion limits. It supplies no physical eta/eta-prime
axial-decay-constant matrix or singlet-source OPE for the current fit.
No corresponding numerical eta target is inferred here.

The INSPIRE-provided BibTeX record is stored unmodified in
`bibliography.bib`; [source navigation](../../sources/Knecht:2001qf/INDEX.md)
identifies the archived arXiv files.
