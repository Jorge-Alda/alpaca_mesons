---
title: Alda neutral-meson generators from mass diagonalization
status: literature
created: 2026-10-13
updated: 2026-10-13
bibtex_keys: [Alda:2026nrd]
provenance: literature
---

# Physical neutral generators

Verified on 2026-10-13 in the downloaded **arXiv:2610.13671v1** source,
submitted 2026-10-09 14:37:51 UTC. The arXiv abstract page identifies v1
as the current version: <https://arxiv.org/abs/2610.13671v1>.
The source archive was obtained from <https://arxiv.org/src/2610.13671>
and stored at `../../sources/Alda:2026nrd/arxiv-2610.13671-source.tar`.
Its SHA256 is
`3322b89e81dabd86e4e6cd03d26bc940668694c12714e197377bd32972e67339`.
The extracted `main_arXiv.tex`, appendix `app:mesonmass`, Eqs.
`eq:fit_mesonic` and `eq:matrices_physmesons`, lines 2230–2323,
is the evidence used below. [@Alda:2026nrd]

This supersedes the earlier local-checkout provenance. All central values
and limitations previously transcribed in this record are also present in
arXiv v1; no numerical correction to those entries was needed. This is a
comparison of this record with the downloaded source, not a claim that the
entire local manuscript is identical to arXiv v1.

The displayed central generators are

$$
\lambda_{\pi^0}=\operatorname{diag}(1.01707,-0.98257,-0.01140),\quad
\lambda_\eta=\operatorname{diag}(0.79999,0.83749,-0.81156),\quad
\lambda_{\eta'}=\operatorname{diag}(0.57061,0.57721,1.15812).
$$

Use the full diagonalization/eigenvectors for precision; these entries are
rounded. Preserve the isospin admixtures and phases. The fit supplies
$B_0m_u=6568$, $B_0m_d=11740$, $B_0m_s=235921$ in MeV squared,
$m_0^2=0.67124$ GeV squared and their correlation matrix.
The appendix describes an effective LO calibration, excludes the eta mass
from the likelihood, and predicts $m_\eta=494.10$ MeV; the quoted parameter
errors exclude chiral truncation effects. [@Alda:2026nrd]

Alda uses $\Phi=\sum\lambda_P P$, $u=\exp(i\Phi/(2f_\pi))$,
$\operatorname{tr}(\lambda_i\lambda_j)=2\delta_{ij}$, $f_\pi=93$ MeV.
The existing RChT kernel uses $\Phi_R=\sum T_P P$ and
$u=\exp(i\Phi_R/(\sqrt2F))$. For the same fields and $F=f_\pi$, their
explicit convention conversion is $T_P=\lambda_P/\sqrt2$. This conversion
must not be confused with independently renormalizing physical states.
Supply $\chi=2\operatorname{diag}(B_0m_u,B_0m_d,B_0m_s)$ to the recorded
trace definitions; the anomalous $m_0^2$ is not part of this quark spurion.
This is the user-selected input prescription; numerical propagation and
predictor integration remain to be performed.
