---
title: Tensor omega-phi mixing source audit
status: literature
created: 2026-10-13
updated: 2026-10-13
bibtex_keys: [Kucukarslan:2006wk]
provenance: literature
---

# Tensor mixing audit

Source: arXiv:hep-ph/0603061v1, downloaded 2026-10-13.
The archive and extracted source are indexed under
[the source folder](../../sources/Kucukarslan:2006wk/INDEX.md).
INSPIRE-HEP record 711934 supplies the unmodified BibTeX entry.

Sections 2–3 use the Ecker tensor kinetic terms and the leading $F_V,G_V$
operators, with a separate omega–phi tensor bilinear and photon exchange.
The source adopts a particular quark-counting mass prescription; it is not a
complete flavor-resolved SU(3)-breaking nonet action for all the repository's
vertices. [@Kucukarslan:2006wk]

Section 6 defines vector interpolating fields by

$$
V_\mu=\partial^\nu V_{\mu\nu}/M_V.
$$

For their mixed transverse two-point function it gives

$$
\Theta(k^2)=\left[2\sqrt2M_\rho(\widehat m-m_s)
+\frac{\sqrt2}{9}e^2F_V^2\right]\frac{k^2}{M_\phi M_\omega}.
$$

This expression belongs to the explicitly derivative-defined vector fields;
its factor $k^2$ does not establish that a tensor-source response vanishes at
$k^2=0$. The paper also states that a physical exchange amplitude has a more
complicated momentum dependence than this correlator numerator.
[@Kucukarslan:2006wk]

The paper reports a small admixture $0.059\pm0.005$ and a departure from ideal
mixing of $3.4\pm0.3$ degrees in its conventions and data treatment (Section 4).
These dated, conditional results are comparison inputs, not a new independent
precision prior. [@Kucukarslan:2006wk]

## Application decisions

The following are modeling decisions, not claims established by this paper:

- Derive our own invariant tensor mass action and its matrix response; retain
  the repository's positive-strange phi convention and source sign.
- Keep strong mass mixing separate from electromagnetic mixing. Explicit
  photon propagation already generates the latter; importing an effective
  electromagnetic bilinear as well would double count it unless the action
  is reorganized consistently.
- A minimal quark-spurion plus singlet-mass truncation is an assumption and
  does not establish a complete tensor/Proca interaction matching.
