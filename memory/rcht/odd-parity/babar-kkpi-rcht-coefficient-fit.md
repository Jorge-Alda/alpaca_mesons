---
title: RChT coefficient profiles from the BaBar isovector K star K spectrum
status: numerical
created: 2026-10-30
updated: 2026-10-30
bibtex_keys: [BaBar:2007ceh, GomezDumm:2010ra, Alda:2026nrd, Ruiz-Femenia:2003jdx, ParticleDataGroup:2024cfk]
provenance: numerical-computation
---

# Observable and coefficient coordinates

This record fits RChT coefficients to BaBar's Table 7 isovector
$e^+e^-\to K^*(892)K$ component [@BaBar:2007ceh]. It uses the crossed
$K^*\to K\gamma^*$ vertex with the isovector source
$$
Q_I=\operatorname{diag}(1/2,-1/2,0).
$$
The stable-$K^*$ quasi-two-body phase space is the one used by BaBar below
$2\,\mathrm{GeV}$. The component is modeled directly; no inclusive $KK\pi$
rate or tau endpoint conversion is inserted.

For the charged $K^*K$ pair, the current form factor is built from
$$
\begin{aligned}
C_I(s)={}&\frac12\Big[(c_1+c_2+8c_3-c_5)m_K^2
+(-c_1+c_2+c_5-2c_6)M_{K^*}^2\\
&\hspace{23mm}+(c_1-c_2+c_5)s+8c_4(m_K^2-m_\pi^2)\Big],\\
H_I(s)={}&\frac12\Big[(d_1+8d_2-d_3)m_K^2+d_3(M_{K^*}^2+s)\Big],\\
F_{VP}(s)={}&\frac{2\sqrt2}{F M_{K^*}M_V}C_I(s)
-\frac{4F_\ell}{F M_{K^*}}\frac{H_I(s)}{D_\rho(s)}
+r_{\rho'}\left[D_{\rho'}(s)-D_{\rho'}(0)\right],\\
F_\ell={}&F_V+2z_V\chi_{+,\ell},\\
\sigma_1(s)={}&\frac{\pi\alpha^2}{3s^3}
\lambda(s,M_{K^*}^2,m_K^2)^{3/2}|F_{VP}(s)|^2.
\end{aligned}
$$
Here $D_\rho$ is the upper-rim rho inverse propagator and $D_{\rho'}$ is
BaBar's four-pion running-width shape with fixed
$\mathcal B_{4\pi}=0.64$. The factor in $\sigma_1$ includes the two
charge-conjugate final states. The stable-$K^*$ normalization was checked
against the charged-minus-neutral isovector form-factor difference and a
direct spin/orientation integral.

The preliminary parameter representative is from
[odd-preliminary-fit.md](odd-preliminary-fit.md). All remaining coefficients
and even/source inputs are held at that representative in the one-coordinate
fits. The new exact real-photon-null coordinate is
$$
c_1\mapsto c_1+\frac{\delta}{2},\qquad
c_5\mapsto c_5+\frac{\delta}{2},\qquad
c_1-c_2+c_5\mapsto c_1-c_2+c_5+\delta.
$$
It leaves the two real-photon VJP numerators unchanged. The $\rho'$ residue
is complex and fitted. Subtracting $D_{\rho'}(0)$ is an explicit matching
assumption that preserves the preliminary real-photon amplitudes; an
unsubtracted variant is also retained.

# Inputs and likelihood assumptions

The coefficient fit reads the 19 bins in
[the tau-overlap CSV](../../../data/vector_sector/babar_2007_isovector_kstar_k_tau_region.csv)
and the 29 bins through $2\,\mathrm{GeV}$ in
[the extended CSV](../../../data/vector_sector/babar_2007_isovector_kstar_k_below_2gev.csv).
The [fit script](../../../scripts/fit_babar_kkpi_rcht_coefficients.py)
records the corresponding preliminary central coefficients and fixed
inputs. It uniformly averages over each published energy bin because the
ISR/acceptance response is unavailable. Statistical errors are diagonal.
BaBar provides no covariance for the additional 15% Dalitz-fit systematic
or the 5.1% total cross-section systematic, so both are varied as common
normalization nuisances and in a diagonal sensitivity model.

The $c_4$ fit uses the first six bins, with all other coefficients fixed to
the preliminary representative and no $\rho'$ pole. The $\delta$ fit uses
all 19 or all 29 bins, fixes $c_4$ to its preliminary value, and profiles the
complex $\rho'$ pole. Separate joint fits vary both $c_4$ and $\delta$.
These fits report conditional likelihoods; the original low-energy
measurement likelihood is not re-evaluated simultaneously.

# Results

In the first-six-bin fit, canonical $c_4$ has the best minimum
$$
c_4=0.0234567,\qquad \chi^2/\mathrm{dof}=1.4444/5,
$$
with local statistical error $0.00224$. The profile is disconnected. The
statistical 68% intervals are
$$
[-0.05321,-0.04951]\;\cup\;[0.02115,0.02563],
$$
and the 95% intervals are
$$
[-0.05540,-0.04693]\;\cup\;[0.01877,0.02762].
$$
The second local minimum is $c_4=-0.0514089$ with
$\chi^2=1.76387$. Under the common 15% and 5.1% normalization nuisances,
the 68% intervals widen to
$[-0.05477,-0.04854]\cup[0.01998,0.02754]$. The sign-related branch must
be retained. This $c_4$ result is conditional on the other preliminary
coefficients and is not a unique global fit over the full coefficient
space.

For 19 bins with $c_4$ fixed at its preliminary value, the complex-$\rho'$
fit gives
$$
c_1-c_2+c_5=0.0024790,\qquad \chi^2/\mathrm{dof}=19.7333/14.
$$
The profile's statistical 68% interval is
$[-0.01179,0.00852]$ and its 95% interval is
$[-0.05132,0.05267]$. The initial preliminary value is $0.02407$. Including
the common systematics barely changes this profile; using diagonal
systematics broadens the 68% interval to approximately
$[-0.02148,0.02009]$.

When $c_4$ and $\delta$ are varied together over the same 19 bins, the best
point is $c_4=0.09485$ and $c_1-c_2+c_5=-0.03830$, with
$\chi^2/\mathrm{dof}=15.481/13$. Their local correlation is $-0.9957$.
The fitted $\rho'$ width is $0.0554\,\mathrm{GeV}$, much narrower than
BaBar's broad global-fit benchmark. This shows that the individual
coefficient values depend strongly on allowing the coefficient directions
and excited pole to move together. Several other joint-fit minima lie within
$\Delta\chi^2<0.2$, including one with negative $c_4$ and positive
$c_1-c_2+c_5$.

Extending the $\delta$-only fit through $2\,\mathrm{GeV}$ gives
$c_1-c_2+c_5=-0.01206$, $\chi^2/\mathrm{dof}=42.169/24$, and puts the
fitted $\rho'$ mass at its lower search bound. The joint 29-bin fit gives
$c_4=-0.06397$, $c_1-c_2+c_5=0.03842$, and
$\chi^2/\mathrm{dof}=31.442/23$, with local correlation $-0.9920$.
Several distinct minima are present. The line-shape and fit-window changes
therefore remain material model dependence; do not collapse these results
to a single coefficient covariance. The 29-bin second joint minimum has
$c_4=0.11176$, $c_1-c_2+c_5=-0.04476$, and
$\chi^2=36.469/23$.

# Remaining null directions and reproduction

This component cannot determine $c_7$, $d_4$, the common shift of
$c_2,c_5,c_6$, or any $g_i$. Other original coefficient directions remain
fixed by the preliminary representative rather than newly measured by
BaBar. The coefficient ranks and nulls are summarized in the
[KK pi rank audit](kkpi-vector-current-coefficient-rank.md). The
nonlinear multistart search does not prove global minimality.

Recreate the two data files and run the deterministic coefficient analysis:
```text
.venv/bin/python scripts/extract_babar_kkpi_isovector_data.py
.venv/bin/python scripts/fit_babar_kkpi_rcht_coefficients.py --output data/vector_sector/babar_kkpi_rcht_coefficient_fit.json
```
The [saved fit output](../../../data/vector_sector/babar_kkpi_rcht_coefficient_fit.json)
contains all minima, profiles, covariance diagnostics, fitted coefficient
vectors, bin predictions and normalization/systematic assumptions.
