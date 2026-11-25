---
title: Loop-size estimate for a to pi+ pi- pi0 in the up-only radial model: diagrams, vertices and result
status: numerical
created: 2026-11-25
updated: 2026-11-25
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh, Manohar:1983md]
provenance: numerical-computation
---

# Loop-size estimate for $a\to\pi^+\pi^-\pi^0$ in the up-only radial model

**Conventions.** Alda:2026nrd and DiLuzio:2026poh, polar fields, connection-free frame, $\epsilon=f_\pi/f_a$, $\langle\lambda_i\lambda_j\rangle=2\delta_{ij}$. Inputs are in [only-up-finite-mixing-inputs.md](only-up-finite-mixing-inputs.md). $\lambda=10^{-3}$ is a labelled default. Kinematics are the symmetric Dalitz point $s_0=s_{+-}=s_{+0}=s_{-0}=(m_a^2+2m_{\pi^+}^2+m_{\pi^0}^2)/3$, a labelled choice. No loop integral is evaluated.

## Estimate

$$
M_1^{\rm est}=\sum_d S_d\,\frac{(m_d^2)^{2-r_d}}{16\pi^2}\;T_d,\qquad R=\frac{|M_1^{\rm est}|}{|M_0|},
$$

where $T_d$ is the value of the diagram $d$ (defined in the next section), and:

- the sum runs over the allowed diagrams $d$ listed below;
- $r_d$ is the number of internal lines of $d$;
- $S_d=1/|{\rm Aut}_d|$: $\tfrac12$ for a closed internal line or for two identical internal lines, $1$ otherwise;
- $m_d=m_a$ if $d$ contains an A line, and $m_d=m_\sigma$ otherwise;
- every product of two momenta that contains the loop momentum is set to $m_d^2$;
- the factors $\epsilon$ and $N$ are inside the vertices and the pole vectors;
- $M_0$ is the exact tree of [the criterion record](murayama-radial-crossing-implicit.md).

The rule $(m_d^2)^{2-r_d}/(16\pi^2)$ comes from the loop measure $d^4k\to m^4/(16\pi^2)$ and from $1/m^2$ for each propagator. The result is dimensionless. It was checked by rescaling every dimensionful input by a common factor $c\in\{3,\,0.37,\,7.1\}$: $R$ is unchanged to $10^{-13}$.

## Legs, vertices and internal lines

**Legs.** A leg is a pair $(\lambda_x,\xi_x)$, as in the [finite-mixing three-pion record](../three-pseudoscalar/alp-finite-mixing-three-pion-amplitudes.md): $\lambda_x$ is the meson matrix, $\xi_x$ the ALP component. A radial leg is written $\sigma$ and carries no meson matrix. Momenta are $k_x$, all incoming.

- Charged pions: $(\lambda_{\pi^+},0)$ and $(\lambda_{\pi^-},0)$, with $\lambda_{\pi^\pm}$ as defined in the [finite-mixing three-pion record](../three-pseudoscalar/alp-finite-mixing-three-pion-amplitudes.md).
- Neutral pion: $(\lambda_{\pi^0},\xi_{\pi^0})$ with $\lambda_{\pi^0}$ the neutral-pion matrix and $\xi_{\pi^0}=1/\sqrt{N(m_{\pi^0}^2)}$, from the pole vector below with $m_{\pi^0}^2$.
- External ALP: the pole vector $v_a=\big(2\epsilon\,\mathbb M^{-1}_{m_a^2}(\mathcal C)/\sqrt{N(m_a^2)},\ 1/\sqrt{N(m_a^2)}\big)$, written as a leg. The label $a$ in the diagrams below means $v_a$.
- Radial exchanges: $\sigma_{\pi^+\pi^-}$ denotes the radial cubic vertex $V(\sigma;\cdot,\cdot)$ whose $\sigma$ leg is the tree exchange into $\pi^+\pi^-$; it carries the factor $V(\sigma;\pi^+,\pi^-)/(m_\sigma^2-s_0)$. $\sigma_{a\pi^0}$ is defined in the same way, with the factor $V(\sigma;a,\pi^0)/(m_\sigma^2-s_0)$ and the exchange into $a\pi^0$.

Here $\mathcal C=B_0\widehat P_{01}$, $\mathbb M_s(X)=B_0\{\widehat S_{00},X\}+\tfrac{m_1^2}{3}\operatorname{tr}(X)\mathbf 1-sX$, and $\mathbb M^{-1}_{s}(\mathcal C)$ is the inverse of $\mathbb M_s$ applied to $\mathcal C$. $N(s)=1+2\epsilon^2\operatorname{tr}\big[(\mathbb M^{-1}_{s}(\mathcal C))^2\big]$.

**Internal lines.** Each internal line places a leg at each of the two vertices it joins.

- ALP line: the leg $v_a$ at each of its two vertices.
- Meson line: the sum over the nine generators $i=0,\dots,8$ of the leg $(\lambda_i,0)$ at each of its two vertices, minus the same with the leg $v_a$ at both vertices. The sum is the $U(3)$ completeness relation; the subtracted term removes the ALP pole, which the ALP line already carries.
- Radial line: the leg $\sigma$ at each of its two vertices.

The diagram value $T_d$ is the sum over the internal generator labels of the product of the vertex functions $V$ of its blocks, with each internal line contributing its legs as above.

**Potential part.** For a vertex with legs $x_1,\dots,x_n$, let $R$ be the set of radial legs, $A$ the set of legs contributing their ALP component, and $P$ the set of legs contributing their meson matrix; a leg can have only the roles it has. Let $r=|R|$, $s=|A|$, $m=|P|$ and $X_{rs}=\widehat S_{rs}$ for even $m$, $X_{rs}=\widehat P_{rs}$ for odd $m$:
$$
V^{(n)}_{\rm pot}=B_0\sum_{(R,A,P)}r!\,s!\,\epsilon^{r+s}\,\frac{(-1)^{\lfloor m/2\rfloor}}{m!}\prod_{x\in A}\xi_x\sum_{\pi\in S_m}\big\langle X_{rs}\,\lambda_{\pi_1}\cdots\lambda_{\pi_m}\big\rangle .
$$
In the up-only model, the only nonzero sources are, with $r\le1$ and $(r,s)\neq(0,0)$,
$$
\widehat S_{rs}+i\widehat P_{rs}=\frac{m_u(-i)^s}{s!}\Big(\tfrac13\mathbf 1+\tfrac12\lambda_3+\tfrac1{2\sqrt3}\lambda_8\Big),
$$
and
$$
\widehat S_{00}=\tfrac13(m_u+m_d+m_s)\,\mathbf 1+\tfrac12(m_u-m_d)\,\lambda_3+\tfrac1{2\sqrt3}(m_u+m_d-2m_s)\,\lambda_8 .
$$
No source with $r\ge2$ is nonzero.

**Four-meson kinetic part** ($n=4$, all legs with meson matrices):
$$
V^{(4)}_{\rm kin}=\frac1{48}\sum_{\pi\in S_4}\big(-k_{\pi_2}\cdot k_{\pi_4}\big)\,\big\langle[\lambda_{\pi_1},\lambda_{\pi_2}][\lambda_{\pi_3},\lambda_{\pi_4}]\big\rangle ,
$$
with every product $k_i\cdot k_j$ that contains the loop momentum replaced by $m_d^2$.

**ALP kinetic part** ($n\le4$):
$$
V_{\rm kin,a}=-\sum_{i<j}\xi_i\xi_j\,k_i\cdot k_j\;C(\text{remaining legs}),\qquad C(\varnothing)=1,\quad C(\sigma)=2\epsilon,\quad C(\sigma,\sigma)=2\epsilon^2,
$$
where $C$ is evaluated on the remaining legs, and $C(\cdot)=1$ for a remaining leg that is not radial.

**Radial cubic vertex.** For legs $x,y$ without a radial component:
$$
V(\sigma;x,y)=\frac1{f_\pi}\Big\{B_0\Big[-\tfrac{\epsilon}{2}\langle\widehat S_{10}\{\lambda_x,\lambda_y\}\rangle+\epsilon^2\big(\xi_x\langle\widehat P_{11}\lambda_y\rangle+\xi_y\langle\widehat P_{11}\lambda_x\rangle\big)+2\epsilon^3\xi_x\xi_y\langle\widehat S_{12}\rangle\Big]-2\epsilon\,\xi_x\xi_y\,k_x\cdot k_y\Big\}.
$$

**Sigma-squared vertex.** $V(\sigma,\sigma;x,y)=-2\epsilon^2\xi_x\xi_y\,k_x\cdot k_y/f_\pi^2$. It has no potential part, because no $\widehat S_{2s}$ source is nonzero.

**Sigma-cubed vertex.** $V_{\sigma^3}=-6\lambda f_a$, from $\mathcal L\supset-\lambda f_a\sigma^3$.

## Diagrams

The 50 diagrams are drawn as graphs in [murayama_loop_diagrams.tex](../../../latex/murayama_loop_diagrams.tex): vertices with their external legs, and internal lines with their types. Each diagram is labelled $k.i$, where $k$ is its class and $i$ its index within the class. The classes are defined by the internal lines they contain.

| class | internal lines | diagrams | scale $m_d$ |
|---|---|---|---|
| 1 | one ALP line, closed on one vertex | 1.1–1.3 | $m_a$ |
| 2 | one ALP line and one meson line | 2.1–2.5 | $m_a$ |
| 3 | two ALP lines | 3.1–3.3 | $m_a$ |
| 4 | one radial line and one meson line (bubbles 4.1–4.6, triangles 4.7–4.14) | 4.1–4.14 | $m_\sigma$ |
| 5 | one radial line and one ALP line (bubbles 5.1–5.4, triangles 5.5–5.12) | 5.1–5.12 | $m_a$ |
| 6 | two radial lines (bubbles 6.1 and 6.7, triangles 6.2–6.4, boxes 6.5–6.6) | 6.1–6.7 | $m_\sigma$ |
| 7 | ALP self-energy on an external leg | 7.1–7.4 | $m_a$ |
| 8 | radial self-energy on an external leg | 8.1–8.2 | $m_\sigma$ |

Each diagram has symmetry factor $S_d=1/|{\rm Aut}_d|$ and $r_d$ internal lines. Diagrams in classes 1–6 have $r_d=1$ to $4$; classes 7 and 8 are two-point corrections with $r_d=1$.

**Radial exchanges.** $\sigma_{\pi^+\pi^-}$ and $\sigma_{a\pi^0}$ denote a radial vertex whose $\sigma$ leg is the tree exchange into the named pair. Its factor is $V(\sigma;\pi^+,\pi^-)/(m_\sigma^2-s_0)$, or $V(\sigma;a,\pi^0)/(m_\sigma^2-s_0)$ for $\sigma_{a\pi^0}$. The diagrams of classes 1 to 6 with such a vertex are drawn with the exchange as an external leg of that vertex.

**Sigma-cubed vertex.** Diagrams 6.1, 6.3, 6.4 and 6.7 contain the vertex $V_{\sigma^3}$. They are marked in the atlas.

**External residues.** Classes 7 and 8 are the self-energy corrections of the external legs. Class 7 is the closed ALP line, nonzero on $\pi^\pm$ and zero on $a$ and $\pi^0$. Class 8 is the closed radial line on the external $a$ and $\pi^0$ legs.

## Truncation

Diagrams with three or more ALP or radial lines are excluded. The selection rules allow five of them: two boxes with two radial, one ALP and one meson line; one triangle with two radial and one ALP line; and two triangles with two ALP and one radial line. The sigma self-energies with pseudoscalar lines are also excluded: they have one $\sigma_{\pi^+\pi^-}$ and one $\sigma_{a\pi^0}$ vertex joined by two ALP, one ALP and one meson, or two meson lines. Their size is not estimated. The retained class-5 bubble is the largest term, so the sum is reliable only if the excluded terms do not cancel it.

## Result ($\lambda=10^{-3}$, symmetric Dalitz point)

| $\kappa_u$ | $m_a$ (GeV) | $m_\sigma$ (GeV) | $M_0$ | $R$ | channel |
|---|---|---|---|---|---|
| 0.10 | 0.24882 | 0.23397 | $-1.670797$ | $1.27\times10^{1}$ | closed |
| 0.15 | 0.35619 | 0.35096 | $-1.550431$ | $7.62\times10^{1}$ | closed |
| 0.20 | 0.45916 | 0.46794 | $-1.610290$ | $7.19\times10^{1}$ | open |
| 0.30 | 0.71104 | 0.70191 | $-2.565123$ | $1.50\times10^{3}$ | open |
| 0.50 | 1.17770 | 1.16985 | $-4.116116$ | $1.21\times10^{4}$ | open |
| 1.00 | 2.34249 | 2.33970 | $-8.130230$ | $2.03\times10^{5}$ | open |

The decay channel opens at $\kappa_u=0.16533$. The first two rows are below threshold and their values are formal. On the grid $\kappa_u\in\{0.1,0.15,0.16533,0.2,0.205,0.21,0.22,0.23,0.25,0.3,0.4,0.5,0.7,1.0\}$, $R\ge12$.

**Contributions** $C_k/|M_0|$ by class from the run (positive values have the opposite sign to $M_0$, which is negative):

| class | diagrams | $\kappa_u=0.3$ | $\kappa_u=1.0$ |
|---|---|---|---|
| 1 | 1.1–1.3 | $+1.321\times10^{1}$ | $+1.781\times10^{3}$ |
| 7 | 7.1–7.4 | $+1.497\times10^{-3}$ | $+1.540\times10^{-3}$ |
| 2 | 2.1–2.5 | $-7.058$ | $-2.202\times10^{1}$ |
| 3 | 3.1–3.3 | $-1.082\times10^{2}$ | $-1.488\times10^{4}$ |
| 4, bubbles | 4.1–4.6 | $-6.612\times10^{1}$ | $-3.542\times10^{2}$ |
| 4, triangles | 4.7–4.14 | $+3.310$ | $-1.198\times10^{1}$ |
| 5, bubbles | 5.1–5.4 | $-1.534\times10^{3}$ | $-1.892\times10^{5}$ |
| 5, triangles | 5.5–5.12 | $+2.122\times10^{2}$ | $+1.008\times10^{3}$ |
| 6.1 | 6.1 | $+5.381\times10^{-5}$ | $+5.299\times10^{-5}$ |
| 6, triangles | 6.2–6.4 | $+6.463\times10^{-1}$ | $+7.522$ |
| 6, boxes | 6.5–6.6 | $+4.832\times10^{-1}$ | $+2.333\times10^{-1}$ |
| 8 | 8.1–8.2 | $-1.399\times10^{1}$ | $-1.788\times10^{3}$ |
| 6.7 | 6.7 | $-1.087\times10^{-10}$ | $-8.454\times10^{-13}$ |

The class-5 bubbles dominate at both points.

**Accidental zeros.** $R<1$ occurs only in two windows, each around a zero of $M_1^{\rm est}$ where the contributions cancel. The zeros are at $\kappa_u=0.207690$ and $0.225347$, with $M_1^{\rm est}=+2.3\times10^{-3}$ and $-1.9\times10^{-3}$ respectively, and $M_0\approx-1.54$ and $-1.75$. $M_1^{\rm est}$ changes sign across a window of width about $10^{-4}$, and $R\gg1$ on both sides. These are not physical crossings. No $\kappa^*$ with $R(\kappa^*)=1$ exists in the open-channel range apart from these windows.

## Checks

- The unit tests confirm: the code's list equals the allowed diagrams, each appearing once, with the stated $S_d$; the vertex closed forms; the charged residue slope against the closed form above; unit rescaling with $R$ unchanged to $10^{-9}$; and frame invariance ([frame-invariance record](murayama-radial-frame-invariance.md)).

The criterion, the exact tree and the caveats are in [murayama-radial-crossing-implicit.md](murayama-radial-crossing-implicit.md). The up-only parameters are in [murayama-radial-up-only-crossing.md](../three-pseudoscalar/murayama-radial-up-only-crossing.md).
