---
title: Nearly symmetric high-energy ALP transition form factor
status: mixed
created: 2026-11-15
updated: 2026-11-15
bibtex_keys: [Alda:2026nrd, Messerli:2025rnv, Ludtke:2020moa, Larin:1993tq, Bass:2002mv, Grozin:2012ic]
provenance: analytical-derivation
---

# Scope and inputs

Use the conventions of [Alda:2026nrd and Messerli:2025rnv](alp-high-energy-tff-conventions.md), natural units, $p=q_1+q_2$, $p^2=m_a^2$, $q=(q_1-q_2)/2$, and spacelike $q^2$ with $|q^2|\gg m_a^2,|q\cdot p|$. The calculation combines the [local current expansion](alp-hard-current-ope.md), [LO chiral matrix element](alp-axial-current-matrix-element.md), [singlet evolution](alp-singlet-axial-current-evolution.md), and the [heavy-flavor triangle](alp-heavy-quark-triangle.md). It keeps the tree-level chiral current, the one-loop heavy triangle and the $O(\alpha_s)$ coefficient of the light axial current. Corrections to the heavy triangle and mixed heavy-light gluon graphs are not calculated. The current matrix element is a named renormalized input; identifying it with the LO chiral expression at a low scale is an approximation.

# Assembly

At leading chiral order for the light sector, with the trace over $u,d,s$ and $Q=\operatorname{diag}(2/3,-1/3,-1/3)$,
$$
\begin{aligned}
F_{a\gamma^*\gamma^*}={}&\frac{c_\gamma+2N_c\langle Q^2c_q^A\rangle}{4\pi^2f_a}
+\sum_{h=c,b}\frac{2N_cQ_h^2}{4\pi^2f_a}
\left[c_{hh}^P+(c_{hh}^A-c_{hh}^P)B_4\right]\\
&-\frac{4f_\pi^2}{f_aq^2}\left(1-\frac{\alpha_s(\mu)}\pi\right)
\left\{\left\langle\left(Q^2-\frac29\mathbf1\right)\mathbb M_a^{-1}(\bm{\mathcal C})\right\rangle_{\mu_0}
+\frac29\left[\text{singlet evolution from }\mu_0\text{ to }\mu\right]\left\langle\mathbb M_a^{-1}(\bm{\mathcal C})\right\rangle_{\mu_0}\right\}+\cdots .
\end{aligned}
$$
Here $B_4$ has the arguments $4m_h^2/m_a^2,4m_h^2/q_1^2,4m_h^2/q_2^2$ of its [definition](alp-loop-function-b4.md), $\mu=\sqrt{-q^2}$, and the ellipsis includes terms beyond the stated chiral, perturbative and local-expansion accuracy. The scale-dependent singlet factor and the convention for the renormalized current are specified in the singlet evolution record; the Wilson coefficient assumes its $O(\alpha_s)$ form applies also to the singlet channel. The displayed singlet decomposition applies directly while the effective theory has three active flavors. Above charm, the light-flavor singlet is not the full active-flavor singlet: its evolution mixes with the charm current, so the three-flavor factor cannot simply be continued with $n_f=4$ or $5$. The active-flavor evolution is derived below. The required threshold values of all active-flavor current matrix elements and their finite matching corrections are named inputs; the [threshold record](alp-heavy-flavor-thresholds.md) supplies the running-coupling convention but its simple singlet-continuity assumption does not determine their individual matching. The expression uses no $\eta$ or $\eta'$ decay constant as an ALP input.

The first line combines the photon contact and direct quark triangles; no heavy axial-current matrix element is added to the heavy triangle. For a common axial rotation, the shifts of $c_\gamma$ cancel those of the direct derivative couplings, and the chiral matrix element depends on the invariant sources. The light mass-dependent hard remainder is proportional to $m_q^2$ and therefore beyond LO chiral counting; its subtraction and cutoff-scheme calculation are given in [heavy-gluon threshold matching](alp-heavy-gluon-threshold-matching.md). It is not claimed to vanish at fixed quark masses. The direct gluonic coefficient is absent among massless dimension-three operators in the symmetric expansion at the order stated by Ludtke:2020moa; this does not determine higher-dimensional source terms or mixed heavy-light gluon contributions.

Below a heavy threshold in all soft invariants, $c_G$ entering $\bm{\mathcal C}$ is matched by $c_G\to c_G+c_{hh}^P$, as derived in the threshold matching record. This local matching is not justified near or above the charm-pair threshold. There the renormalized vacuum-to-ALP light-current matrix elements, including charm-induced effects, are named nonperturbative inputs; the free-charm triangle also omits charmonium effects. Consequently the displayed three-flavor chiral formula is conditional in that region and is not a determination of those inputs. The high-energy expansion is not a result for a singly virtual photon or a fixed unequal ratio of virtualities.

# Evolution in each active-flavor interval

Let $\mu_1$ and $\mu_2$ lie in the same interval with $n_f$ active flavors. The singlet and traceless flavor currents evolve independently. Applying their projectors gives, for each active flavor,
$$
\begin{aligned}
\langle0|\bar q_i\gamma^\mu\gamma_5q_i|a\rangle_{\mu_2}
={}&\langle0|\bar q_i\gamma^\mu\gamma_5q_i|a\rangle_{\mu_1}\\
&+\frac{1}{n_f}\left\{\exp\left[\frac{2n_f}{\pi(11-2n_f/3)}
\left(\alpha_s(\mu_2)-\alpha_s(\mu_1)\right)\right]-1\right\}
\sum_{j=1}^{n_f}\langle0|\bar q_j\gamma^\mu\gamma_5q_j|a\rangle_{\mu_1}.
\end{aligned}
$$
This is the solution with the leading singlet anomalous dimension and one-loop coupling, as in the [singlet record](alp-singlet-axial-current-evolution.md). It applies separately in each interval; at a threshold use the matched current values in the new theory as inputs. At $O(\alpha_s)$ fixed-order accuracy the current matching is tree level; any finite higher-order matching is not supplied here. In a four-flavor interval, the evolution of a light current contains the charm-current matrix element with coefficient one quarter of the curly bracket. A light singlet alone therefore does not evolve multiplicatively. Weighting only the $u,d,s$ currents by their charge squares gives a change equal to $2/(3n_f)$ times the curly bracket times the sum of all active-flavor matrix elements. The script checks this charge-weighted result and the projectors for $n_f=3,4,5$, their composition, and the four-flavor charm coefficient. Mixed heavy-light terms associated with this evolution should be kept or omitted consistently with the perturbative accuracy; they are not supplied by the free heavy triangle.

# Heavy-current threshold matching at two loops

Bass:2002mv gives the following operator matching at leading power in the heavy mass and through its first nonzero perturbative order, in its $\overline{\rm MS}$ current convention (see the [extract](../../extracts/Bass:2002mv/heavy-singlet-axial-current-matching.md)):
$$
\left(\bar h\gamma^\mu\gamma_5h\right)_{n_f}=
\left(\frac{\alpha_s^{(n_f)}}{\pi}\right)^2
\left[\log\frac{m_h}{\mu}+\frac18\right]
\sum_{i=1}^{n_f-1}\left(\bar q_i\gamma^\mu\gamma_5q_i\right)_{n_f-1}
+O(\alpha_s^3,1/m_h).
$$
Thus in the charge-weighted electromagnetic-current expansion this term contributes the stated coefficient multiplied by $Q_h^2$ to each light axial current. At $\mu=m_h$ the coefficient is $(\alpha_s/\pi)^2/8$. The script derives this charge-weighted insertion and verifies the threshold value. The relation applies only when all soft momenta are small compared with $m_h$; it does not describe $m_a$ near or above the charm-pair threshold. It is beyond the fixed $O(\alpha_s)$ Wilson-coefficient accuracy retained above, and the corresponding mixed-gluon terms must be included consistently if it is used. In particular, do not add this contribution on top of a heavy triangle without a matching calculation that subtracts their common perturbative part. The relation does not fix finite matching of the light-singlet current or the nonperturbative charm-region matrix elements. Grozin:2012ic discusses the singlet and nonsinglet light-current matching separately; its [extract](../../extracts/Grozin:2012ic/axial-current-decoupling.md) records the perturbative diagram orders and the absence of a finite light-singlet coefficient in that source.

# Verification

Run `.venv/bin/python scripts/derive_alp_symmetric_tff.py --test` and `.venv/bin/python scripts/derive_alp_hard_current_ope.py --test`. The first script checks light and heavy source rotations, triangle limits, the soft/hard subtraction, photon exchange and both Ward contractions. The second independently checks the local Wick contraction, Clifford algebra and scalar projection. All results refer to $m_a>0$; the scalar matrix element cancels the projection's apparent $m_a^{-2}$.
