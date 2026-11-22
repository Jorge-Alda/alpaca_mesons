---
title: One-loop power counting for ALP lines in the LO ALP-meson theory at finite f_pi/f_a
status: mixed
created: 2026-11-21
updated: 2026-11-22
bibtex_keys: [Alda:2026nrd, Kaiser:1998ds, Gasser:1983yg, Manohar:1983md]
provenance: analytical-derivation
---

# One-loop power counting for ALP lines in the LO ALP-meson theory at finite $f_\pi/f_a$

## Conventions and units

Units GeV. $x=a/f_a$, $U=u^2=\exp(i\Phi/f_\pi)$, $\Phi=\sum_i\lambda_iP_i$, $\operatorname{tr}(\lambda_i\lambda_j)=2\delta_{ij}$, $N_f=3$, notation of [@Alda:2026nrd]. Neutral fields in the flavour basis: $\lambda_q=\sqrt2\,E_{qq}$ ($P_u,P_d,P_s$); charged pairs $\Phi_{qq'}=\sqrt2\,\phi_{qq'}$. Vertices $V^L$ are derivatives of the Lagrangian (amputated, without the factor $i$), all momenta incoming.

Loop function (dimensional regularization; $\int d^dk\,1=0$, so quadratic divergences vanish):
$$
A_0(m^2)=-\mu^{4-d}\!\int\!\frac{d^dk}{(2\pi)^d}\frac{i}{k^2-m^2+i0}
\;\to\;\frac{m^2}{16\pi^2}\Big[1-\ln\frac{m^2}{\mu^2}\Big],
$$
after $\overline{\rm MS}$ subtraction of $\frac{m^2}{16\pi^2}[\frac{2}{4-d}-\gamma_E+\ln4\pi]$. Then $\langle\phi^2\rangle=-A_0$. Other subtraction constants change only the polynomial part, which is absorbed by the local counterterms of the next order [@Gasser:1983yg]. For size estimates $|A_0(m^2)|\to m^2/(16\pi^2)$ (approximation: the logarithm is dropped, i.e. $\mu$ near the mass in the loop); this is the loop factor of [@Manohar:1983md].

## Inputs

- **Dressed frame (input hypothesis, tested in a separate tree-unitarity calculation).** Since $x$ is a single scalar, the ALP currents are pure gauge and are removed by an $x$-dependent chiral rotation ([general redefinition record](alp-general-chiral-redefinition-invariants.md)). The LO Lagrangian is then
$$
\mathcal L=\frac{f_\pi^2}{4}\langle u_\mu u^\mu+\chi_+\rangle-\frac{f_\pi^2m_0^2}{4N_f}X_{TC}^2+\frac12(\partial a)^2-\frac12m_{a,0}^2a^2,
$$
with $\chi=2B_0\widehat M(x)$, $u_\mu$ free of sources, $X_{TC}=-\operatorname{tr}\Phi/f_\pi-2\widehat\Theta(x)$ ([dressed topological record](alp-dressed-topological-invariant.md)). Kinetic terms are canonical; ALP-meson mixing is only in $M^2$. No bare $a^4$ term is included (assumption).
- **Source convention (choice, stated).** $M(x)=m_q+ix\{c_q^P,m_q\}$, constant diagonal $c_q^A=\operatorname{diag}(c^A_u,c^A_d,c^A_s)$, $c_q^P=\operatorname{diag}(c^P_u,c^P_d,c^P_s)$, $C_V=0$. Then $W_R=e^{ixc_q^A}$, $W_L=e^{-ixc_q^A}$ and
$$
\widehat M(x)=e^{-2ixc_q^A}(1+2ixc_q^P)\,m_q,\qquad
\widehat M_n=(-2i)^n\,m_q\,(c_q^A)^{n-1}(c_q^A-n\,c_q^P)\ (n\ge1),\qquad
\widehat\Theta(x)=c_{GA}\,x .
$$
**Correction to the task statement:** $\widehat M_n$ for $n\ge2$ is *not* a function of $c_q^{AP}=c_q^A-c_q^P$ alone; e.g. $\widehat M_2=-4m_q(c_q^{A})^2+8m_qc_q^Ac_q^P=-4m_q[(c_q^{AP})^2-(c_q^P)^2]$. Only $\widehat M_1=-2i\,m_qc_q^{AP}$ ([first-order record](alp-dressed-mass-first-order-invariants.md)) depends on $c_q^{AP}$ alone. Two representatives with the same $c_q^{AP}$ differ at second order: $c_q^P=0$ ($\widehat M=m_qe^{-2ixc_q^A}$, a pure phase) and $c_q^A=0$ ($\widehat M$ linear, $\widehat M_{n\ge2}=0$). The latter has no $a^2$ quark-mass term and becomes unstable at large $f_\pi/f_a$ unless $m_{a,0}^2$ or the topological term compensates ([stability-bound record](alp-lo-bare-alp-mass-bound.md)). The illustrative numbers below use $c_q^P=0$, $c_q^A=c\,\mathbf 1$.
- Quadratic blocks and pole solver: [finite-mixing record](finite-alp-meson-mixing.md), [pole solver](alp-lo-pole-solver.md), [joint vertex expansion](alp-lo-joint-vertex-expansion.md); loop structure of mixed two-point functions: [loop LSZ record](../../nlo/loop-corrections-to-alp-meson-lsz-mixing.md).

## L1. Vertex scaling

For neutral (commuting) fields the dressed potential is, exactly,
$$
\mathcal L\supset f_\pi^2B_0\sum_q m_q\big[\cos\vartheta_q+2c^P_q x\sin\vartheta_q\big]-\frac{f_\pi^2m_0^2}{4N_f}\Big(\frac{\sqrt2\sum_qP_q}{f_\pi}+2c_{GA}x\Big)^2,\qquad
\vartheta_q=\frac{\sqrt2P_q}{f_\pi}+2c^A_qx .
$$
In general (any flavour structure) the term of $\frac{f_\pi^2}{4}\langle\chi_+\rangle$ with $n_a$ ALPs and $n_P$ mesons gives
$$
V^L_{n_a;i_1\dots i_{n_P}}=\frac{B_0\,f_\pi^{2-n_P}f_a^{-n_a}}{2^{1+n_P}}\sum_{\sigma\in S_{n_P}}\sum_{k=0}^{n_P}\frac{(-i)^{n_P}\langle\lambda_{\sigma_1}\cdots\lambda_{\sigma_k}\widehat M_{n_a}\lambda_{\sigma_{k+1}}\cdots\lambda_{\sigma_{n_P}}\rangle+i^{n_P}\langle\lambda_{\sigma_1}\cdots\lambda_{\sigma_k}\widehat M_{n_a}^\dagger\lambda_{\sigma_{k+1}}\cdots\lambda_{\sigma_{n_P}}\rangle}{k!\,(n_P-k)!},
$$
momentum independent. For neutral legs this reduces to $V^L_{q^{n_P}a^{n_a}}=2^{n_P/2}B_0f_\pi^{2-n_P}f_a^{-n_a}\operatorname{Re}[(-i)^{n_P}(\widehat M_{n_a})_{qq}]$. In this convention $\widehat M_n$ is real for even $n$ and imaginary for odd $n$, so only $n_a+n_P$ even survives: there are no cubic vertices, and the one-loop two-point functions are pure tadpoles. The kinetic term $\frac{f_\pi^2}{4}\langle u_\mu u^\mu\rangle$ contains no ALP ($n_a=0$ only; its quartic vertices vanish for commuting neutral fields). With $\widehat\Theta$ linear in $x$, $X_{TC}^2$ is quadratic in the fields: it contributes to $M^2$ only. A nonzero $\widehat\Theta_2$ would add cubic and higher vertices.

**Superseded (2026-11-22):** the parity statement and the absence of cubic vertices above hold only in the representative with constant diagonal $c_q^A$, $c_q^P$. The vertex formula in the invariants $\widehat S_n$, $\widehat P_n$, $\widehat\Theta_n$, with cubic vertices for general sources, is in the [dressed-vertex power-counting record](alp-dressed-vertex-power-counting-invariants.md).

## L2. One-loop tadpoles

Combinatorics (checked against the second derivative of the Coleman-Weinberg potential): with $D^{-1}(s)=\sum_n v^{(n)}v^{(n)\mathsf T}/(s-m_n^2)$,
$$
\delta M^2_{IJ}=\frac12\sum_nV^L_{IJKL}\,v^{(n)}_Kv^{(n)}_L\,A_0(m_n^2),\qquad
\Gamma(s)=s-M^2+\Sigma(s),\quad \Sigma_{IJ}(s)=-\frac12\sum_nV^L_{IJKL}(s,k^2\!=\!m_n^2)\,v^{(n)}_Kv^{(n)}_L A_0(m_n^2),
$$
where for momentum-dependent vertices the legs are $(p,-p,k,-k)$, $k\cdot k\to m_n^2$ and $p\cdot k\to0$. A complex charged field contributes $-V^L_{IJ\phi\phi^*}A_0$. Check: for SU(2) this gives $m_\pi^2=m^2[1+\frac{m^2}{32\pi^2f_\pi^2}\ln\frac{m^2}{\mu^2}]+\dots$ on shell.

**Charged-meson lines (kinetic vertices included, exact; charged mesons do not mix with the ALP):**
$$
\Sigma_{qq}(s)=\frac{s-2B_0m_q}{3f_\pi^2}\sum_{r\ne q}A_0^{qr},\quad
\Sigma_{qr}(s)=-\frac{s+2B_0(m_q+m_r)}{3f_\pi^2}A_0^{qr},\quad
\delta M^2_{qa}=\frac{2\sqrt2B_0}{3f_\pi f_a}\sum_{r\ne q}\big(2m_qc^{AP}_q+m_rc^{AP}_r\big)A_0^{qr},
$$
with $A_0^{qr}=A_0(B_0(m_q+m_r))$; the mixed entry depends on $c_q^{AP}$ only (it uses $\widehat M_1$).

**Unmixed approximation (labelled: ALP and meson propagators diagonal in $(P_u,P_d,P_s,a)$; $m_0$ neglected in internal meson lines; $m_a^2\simeq M^2_{aa}$).**

(a) ALP line in a meson entry ($q,q'$ neutral or charged):
$$
\frac{\delta M^2_{qq'}}{B_0(m_q+m_{q'})}=\frac{2\big[m_qc^A_q(c^A_q-2c^P_q)+m_{q'}c^A_{q'}(c^A_{q'}-2c^P_{q'})\big]}{m_q+m_{q'}}\,\frac{A_0(m_a^2)}{f_a^2}\;\xrightarrow{c^P=0,\,c^A=c}\;2c^2\frac{A_0(m_a^2)}{f_a^2}.
$$
(b) Meson lines in $M^2_{aa}$:
$$
\delta M^2_{aa}=\frac{4B_0}{f_a^2}\sum_{q,r}m_qc^A_q(c^A_q-2c^P_q)\,A_0\big(B_0(m_q+m_r)\big),\qquad
\frac{\delta M^2_{aa}}{M^2_{aa}|_{\rm quark}}=\frac{\sum_{q,r}m_qc^A_q(c^A_q-2c^P_q)A_0(B_0(m_q+m_r))}{f_\pi^2\sum_qm_qc^A_q(c^A_q-2c^P_q)},
$$
with $M^2_{aa}|_{\rm quark}=\frac{4B_0f_\pi^2}{f_a^2}\sum_qm_qc^A_q(c^A_q-2c^P_q)=-\frac{f_\pi^2}{2f_a^2}\operatorname{Re}\operatorname{tr}\chi_2$.
(c) ALP line in $M^2_{aa}$:
$$
\delta M^2_{aa}=\frac{f_\pi^2B_0}{2f_a^4}\operatorname{Re}\operatorname{tr}\widehat M_4\,A_0(m_a^2)=\frac{8f_\pi^2B_0}{f_a^4}\sum_qm_q(c^A_q)^3(c^A_q-4c^P_q)\,A_0(m_a^2).
$$
**Topological contribution.** $M^2_{aa}\supset 2f_\pi^2m_0^2c_{GA}^2/(N_ff_a^2)$ gives the ALP mass but no self-coupling: $V^L_{aaaa}=(f_\pi^2B_0/f_a^4)\operatorname{Re}\operatorname{tr}\widehat M_4$ is independent of $c_{GA}$, $m_0$ and $m_{a,0}$ (unit-tested). So the ALP mass and the ALP self-coupling are independent; in particular
$$
\frac{\delta M^2_{aa}|_{\rm ALP\ line}}{m_a^2}=\frac{V^L_{aaaa}}{32\pi^2}\ \text{(with }|A_0|\to m_a^2/16\pi^2),
$$
whatever generates $m_a$.

**Exact in mixing.** The script evaluates $\delta M^2$ by the pole sum and separates the ALP-like pole (largest $|v_a^{(n)}|$) from meson-like poles. For large $f_\pi/f_a$ the ALP-like pole is heavy, $m_a^2\to M^2_{aa}$, its meson components fall like $f_a/f_\pi$, and its contribution to meson entries approaches (a) (agreement to $2\times10^{-3}$ at $f_\pi/f_a=30$). Tree-level remark: although the mixing angle falls like $f_a/f_\pi$, the meson-like poles shift by $O(1)$ (independent of $f_\pi/f_a$) through $-M^2_{qa}M^2_{q'a}/M^2_{aa}$; e.g. with $c_q^A=\mathbf1$, $c_q^P=0$, $c_{GA}=0$, $m_{a,0}=0$ and the illustrative inputs below, the meson-like poles are $(0.0195,0.250,0.968)\ \mathrm{GeV}^2$ at $f_\pi/f_a=0.01$ and $(0.0195,0.054,0.730)\ \mathrm{GeV}^2$ at $f_\pi/f_a=30$. With $c_q^P=0$, $c_G=0$, $m_{a,0}=0$ the Lagrangian is shift symmetric and one pole is exactly massless.

## L3. S-wave bubbles and unitarity

For a constant $P P\to aa$ vertex $V_{PPaa}$, the $s$-channel bubble is $\mathcal M=\frac12V_{PPaa}^2B_0(s)/(16\pi^2)$ with $\operatorname{Im}B_0=\pi\sqrt{1-4m_a^2/s}$ (from the Feynman-parameter integral). With $a_0=\mathcal M/(16\pi)$ for an angle-independent amplitude,
$$
\operatorname{Im}a_0(PP\to PP)\big|_{aa}=\sqrt{1-\tfrac{4m_a^2}{s}}\;\tfrac12\,|a_0(PP\to aa)|^2,
$$
the factor $\frac12$ being the identical-particle factor (verified symbolically, residual $0$).

Loop-to-tree ratios for $P_qP_q\to P_qP_q$ (flavour basis, unmixed; tree $V^L_{qqqq}=4B_0m_q/f_\pi^2$; vertices $V^L_{qqqa}=4\sqrt2B_0m_qc^{AP}_q/(f_\pi f_a)$, $V^L_{qqaa}=8B_0m_qc^A_q(c^A_q-2c^P_q)/f_a^2$), with $m_{qq}^2=2B_0m_q$ and $B_0(s)$ the dimensionless bubble function:

| intermediate state | loop / tree | loop factor |
|---|---|---|
| $PP$ | $\dfrac{m_{qq}^2}{16\pi^2f_\pi^2}B_0$ | $1/(16\pi^2f_\pi^2)$ |
| $aP$ | $4(c^{AP}_q)^2\Big(\dfrac{f_\pi}{f_a}\Big)^2\dfrac{m_{qq}^2}{16\pi^2f_\pi^2}B_0=4(c^{AP}_q)^2\dfrac{m_{qq}^2}{16\pi^2f_a^2}B_0$ | $1/(16\pi^2f_a^2)$ |
| $aa$ | $4[c^A_q(c^A_q-2c^P_q)]^2\Big(\dfrac{f_\pi}{f_a}\Big)^4\dfrac{m_{qq}^2}{16\pi^2f_\pi^2}B_0$ | $\dfrac{1}{16\pi^2f_a^2}\cdot\dfrac{f_\pi^2}{f_a^2}$ |

Rule: since $V_{n_a,n_P}\propto f_\pi^{2}\,f_\pi^{-n_P}f_a^{-n_a}$, a one-loop diagram relative to tree carries one extra $f_\pi^2$ (extra vertex) and $1/f_X^2$ for each end-pair of internal line $X$, times $1/16\pi^2$: bubble ratio $\sim m^2f_\pi^2/(16\pi^2f_{X_1}^2f_{X_2}^2)$; tadpole ratio $\sim m_X^2/(16\pi^2f_X^2)$ with $m_X$ the mass in the loop. The mass scale in the bubbles comes from the vertices ($B_0m_q$), not from $m_a$: tuning $m_a$ does not reduce the $aP$ and $aa$ bubbles.

## Result table (tadpoles: $|A_0|\to m^2/16\pi^2$; bubbles: $|B_0|\to1$)

**Superseded (2026-11-22):** the closed forms with $c_q^A$, $c_q^P$ in this table and in the L3 table are representative dependent. The loop classes with their powers of $\epsilon$ and invariant contractions are in the [dressed-vertex power-counting record](alp-dressed-vertex-power-counting-invariants.md); the powers of $f_\pi/f_a$ below agree with it.

| loop class | loop factor | power of $f_\pi/f_a$ | closed form (ratio to tree) |
|---|---|---|---|
| meson line, meson entry | $1/(16\pi^2f_\pi^2)$ | 0 | $\sim m_P^2/(16\pi^2f_\pi^2)$; charged lines above |
| ALP line, meson entry | $1/(16\pi^2f_a^2)$ | 2 at fixed $m_a$; 4 if $m_a^2\propto(f_\pi/f_a)^2$ | $2c^A_q(c^A_q-2c^P_q)\,m_a^2/(16\pi^2f_a^2)$ |
| meson line, ALP entry | $1/(16\pi^2f_\pi^2)$ | 0 | (b) |
| ALP line, ALP entry | $1/(16\pi^2f_a^2)$ | 4 (for any origin of $m_a$) | $V^L_{aaaa}/(32\pi^2)=\dfrac{4(f_\pi/f_a)^4\sum_q(c^A_q)^3(c^A_q-4c^P_q)m_{qq}^2}{16\pi^2f_\pi^2}$ |
| $aP$ bubble | $1/(16\pi^2f_a^2)$ | 2 | table of L3 |
| $aa$ bubble | $f_\pi^2/(16\pi^2f_a^4)$ | 4 | table of L3 |

Expected structures: meson loops $\sim m_P^2/(16\pi^2f_\pi^2)$ — confirmed. One ALP + one meson line $\sim(f_\pi/f_a)^2c^2m^2/(16\pi^2f_\pi^2)$ — confirmed for the $aP$ bubble (coefficient $4(c^{AP})^2$); it does not occur in the two-point functions (no cubic vertices). Two ALP lines up to $(f_\pi/f_a)^4c^4$ — confirmed ($aa$ bubble, ALP line in $M^2_{aa}$), but the coefficient is $[c^A(c^A-2c^P)]^2$, not $(c^{AP})^4$. ALP line in a meson entry $\sim c^2m_a^2/(16\pi^2f_a^2)$ — confirmed, coefficient $2c^A_q(c^A_q-2c^P_q)$, which vanishes in the representative $c_q^A=0$.

## L4. Regime map

Solve ratio $=1$ for $f_\pi/f_a$ (closed forms; $c^P=0$, $c^A_q=c$):

- (i) untuned $m_{a,0}$, $m_a^2\simeq 2c^2(f_\pi/f_a)^2\sum_qm_{qq}^2=2c^2(f_\pi/f_a)^2(m_\pi^2+2m_K^2)$: ALP line in a meson entry, $\big(f_\pi/f_a\big)_*=\big[16\pi^2f_\pi^2/(4c^4(m_\pi^2+2m_K^2))\big]^{1/4}$.
- (ii) $m_a^2\simeq\frac{2}{3}m_0^2c_{GA}^2(f_\pi/f_a)^2$: $\big(f_\pi/f_a\big)_*=\big[3\cdot16\pi^2f_\pi^2/(4c^2c_{GA}^2m_0^2)\big]^{1/4}$.
- (iii) $m_a$ fixed (bare mass tuned, $m_{a,0}^2=m_a^2-M^2_{aa}|_{f_\pi/f_a}$): $\big(f_\pi/f_a\big)_*=4\pi f_\pi/(\sqrt2\,c\,m_a)$.
- ALP line in $M^2_{aa}$ (all cases): same as (i). $aP$ bubble: $4\pi f_\pi/(2c\,m_{qq})$. $aa$ bubble: $[16\pi^2f_\pi^2/(4c^4m_{qq}^2)]^{1/4}$. Meson lines: no crossing (needs $m_P\sim4\pi f_\pi$).

**Illustrative numbers** (inputs: $m_\pi=0.13957$, $m_K=0.49368$, $f_\pi=0.09207$ GeV, $c=1$, $c_{GA}=3$, $m_0^2=m_\eta^2+m_{\eta'}^2-2m_K^2=0.730\ \mathrm{GeV}^2$ from the LO trace relation, $m_a=1$ GeV in (iii), $m_{ss}^2=2m_K^2-m_\pi^2$): (i) $0.90$; (ii) $0.63$; (iii) $0.82$; ALP line in $M^2_{aa}$ $0.90$; $aP$ bubble $4.1$ ($u\bar u$), $0.85$ ($s\bar s$); $aa$ bubble $2.0$ ($u\bar u$), $0.92$ ($s\bar s$). These are illustrative, not predictions. Recomputation with physical inputs ($m_{\pi^0}$, $B_0=m_{\pi^0}^2/(m_u+m_d)$, $\mu=0.77$ GeV, with and without $\overline{\rm MS}$ logarithms): [physical-input crossing record](alp-loop-crossing-physical-inputs.md); the numbers in this paragraph were illustrative only. Conclusion: with $c=O(1)$ ALP loops reach the size of tree terms at $f_\pi/f_a\approx0.6$–$2$; for $f_\pi/f_a\gg1$ the ALP-loop expansion is controlled by $c\,f_\pi/f_a$ (and by $c_{GA}f_\pi/f_a$ for the mass in case (ii)), not by $f_\pi/f_a$ alone.

The tadpole-ratio coefficient $2c^A_q(c^A_q-2c^P_q)$ of (a) holds only in the linear-spurion representative and for an off-shell entry; the frame-invariant pole-mass version, with the dressed invariants, is in the [invariant crossing record](alp-invariant-loop-crossing-meson-poles.md).

## L6. Derivative frame versus dressed frame

Toy model (one meson + ALP, derivative frame $\frac12(\partial P-v\partial a)^2$, $v=\sqrt2cf_\pi/f_a$; dressed $P=\hat P+va$): the on-shell one-loop pole shifts $v^{(n)\mathsf T}\delta M^2v^{(n)}$ agree in both frames, symbolically at $O(v^2)$ and numerically exactly; the off-shell $aa$ entry of $\delta M^2$ differs. In this toy the map is linear, so the check is a congruence. For the full nonlinear chiral rotation the equality of on-shell results is an **assumption**: S-matrix invariance under field redefinitions, trivial Jacobian in dimensional regularization, anomaly already contained in $\widehat\Theta$. Off-shell Green functions and the counterterms differ between frames.

## Reproducible calculation

[`derive_alp_loop_power_counting.py`](../../../scripts/derive_alp_loop_power_counting.py) (imports `derive_alp_lo_joint_vertices.vertex` and `solve_alp_meson_spectrum.solve_spectrum_numerical`). Tests in [`test_derive_alp_loop_power_counting.py`](../../../scripts/test_derive_alp_loop_power_counting.py): $\widehat M_n$ against the Taylor series; closed-form vertex against the general joint-vertex code; momentum independence; no ALP vertex from the kinetic term; linear topological source only quadratic; neutral closed Lagrangian against general vertices; factor $\frac12$ against Coleman-Weinberg; SU(2) chiral logarithm; unmixed (a)-(c) against tensor contraction plus charged loops; pole sum at large $f_\pi/f_a$; topological term without self-coupling; massless pole at the shift-symmetric point; bubble unitarity; regime crossings; L5 amplitudes; L6 frame comparison.

```sh
.venv/bin/python scripts/derive_alp_loop_power_counting.py
.venv/bin/python -m unittest discover -s scripts -p test_derive_alp_loop_power_counting.py
```
