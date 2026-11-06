---
title: Conditional rho-to-two-pion width and self-energy conventions
status: derived
created: 2026-11-06
updated: 2026-11-06
bibtex_keys: [Ecker:1988te, Ecker:1989yg, Jamin:2000wn, ParticleDataGroup:2024cfk, GomezDumm:2013sib, Watson:1952ji, Omnes:1958hv]
provenance: analytical-derivation
---

# Inputs, conventions, and assumptions

This derivation obtains the charged two-pion P-wave width from the leading
even-parity RChT $VPP$ interaction and compares distinct off-shell
continuations. The Lagrangian and antisymmetric-tensor conventions are those
of [the RChT basics](../../rcht/foundations/rcht-basics.md) and [the generic
$VPP$ vertex record](../../rcht/foundations/rcht-vpp-vpa-photon-pa-expansion.md)
[@Ecker:1988te; @Ecker:1989yg]. Use
$$
\mathcal L_{VPP}=\frac{iG_V}{\sqrt2}\langle V_{\mu\nu}u^\mu u^\nu\rangle,
\qquad U=\exp(i\sqrt2\Phi/F),
$$
with the unit-normalized generators
$$
T_\rho=\frac{\operatorname{diag}(1,-1,0)}{\sqrt2},\qquad
\lambda_+=|u\rangle\langle d|,\qquad
\lambda_-=|d\rangle\langle u|.
$$
Natural units, metric $(+---)$, relativistic states, $S=1+iT$, and the upper
boundary $s+i0$ are used. $F$, $G_V$, masses, and widths are in GeV; $s$ and
self-energies are in GeV$^2$. The physical pole mass $M_\rho$ and the tensor
mass $M_{\mathrm{tensor}}$ are kept distinct, as in the [on-shell $G_V$
width record](../../rcht/foundations/vector-gv-fit.md).

The cut includes only $\pi^+\pi^-$ in exact light-quark isospin and retains
the leading $G_V$ vertex. The first continuation below uses a unit transverse
spin-one projector,
$$
\frac{i(q_\mu\varepsilon_\nu-q_\nu\varepsilon_\mu)}{\sqrt{s}},
\qquad s=q^2>0.
$$
This is a stated virtual-projector convention. It is not the matrix element
of an unstable rho as an external LSZ state, and it does not by itself
define the off-shell tensor-field self-energy.

# Unit-projector width from the cut

The executable
[`derive_alp_rho_pipi_self_energy.py`](../../../scripts/derive_alp_rho_pipi_self_energy.py)
first evaluates the charged-pion flavor trace,
$$
\operatorname{tr}\!\left(T_\rho[\lambda_+,\lambda_-]\right)=\sqrt2,
$$
and expands the tensor vertex. Apart from its common phase, the tensor
coefficient is $2G_V/F^2$. In the rho rest frame the tensor contraction
with a unit transverse polarization gives the current coefficient
$$
-\frac{G_V\sqrt{s}}{F^2}
$$
of $\varepsilon\cdot(p_+-p_-)$. The transverse polarization sum is
$$
- (p_+-p_-)^2+\frac{[q\cdot(p_+-p_-)]^2}{s}=4q_{\pi\pi}^2,
\qquad
q_{\pi\pi}(s)=\frac12\sqrt{s-4M_\pi^2}.
$$
The integrated two-body phase space is $q_{\pi\pi}/(4\pi\sqrt{s})$.
Combining it with the spin average and flux factor yields
$$
\Gamma_{\rho\to\pi\pi}^{\mathrm{unit}}(s)
=\frac{G_V^2 q_{\pi\pi}^3(s)}{6\pi F^4}
=\frac{G_V^2}{48\pi F^4}(s-4M_\pi^2)^{3/2}
\theta(s-4M_\pi^2).
$$
Thus the cut turns on as $(s-4M_\pi^2)^{3/2}$, as required for a two-body
P wave [@Jamin:2000wn]. No identical-particle factor is present for the
distinct charged pions.

# Off-shell continuations are not unique

The on-shell width does not choose the continuation away from $s=M_\rho^2$.
The script derives and compares the following conventions:

1. **Unit transverse spin-one projector:** the width above, proportional to
   $q_{\pi\pi}^3(s)$.
2. **Fixed tensor-mass continuation:** continuing the tensor-to-Proca
   normalization with $M_{\mathrm{tensor}}$ fixed multiplies the unit-projector
   width by $s/M_{\mathrm{tensor}}^2$. At the rho pole this reproduces the
   stored on-shell RChT formula
   $$
   \frac{\Gamma_{\rho\to\pi\pi}(M_\rho^2)}{G_V^2}
   =\left|\operatorname{tr}\!\left(T_\rho[\lambda_+,\lambda_-]\right)\right|^2
   \frac{M_\rho^2q_{\pi\pi}^3(M_\rho^2)}
   {12\pi F^4M_{\mathrm{tensor}}^2}.
   $$
3. **Pole-normalized phenomenological P-wave shape:** the line-shape record
   gives the alternative
   $$
   \Gamma_{\rho\to\pi\pi}^{\mathrm{phen}}(s)
   =\Gamma_{\rho\to\pi\pi}(M_\rho^2)\frac{M_\rho^2}{s}
   \left[\frac{q_{\pi\pi}(s)}{q_{\pi\pi}(M_\rho^2)}\right]^3,
   $$
   for $s\geq4M_\pi^2$ [the width and line-shape comparison](../../rcht/foundations/vector-width-lineshape-prescriptions.md)
   [@Jamin:2000wn; @ParticleDataGroup:2024cfk].

All three share the P-wave threshold and match their own pole normalization,
but their off-shell dependence differs. The existing line-shape record
proposes the pole-normalized expression as a practical model; the tensor
vertex does not prove it equivalent to either RChT continuation. No unique
off-shell width is selected here.

# Self-energy sign, dispersion, and subtraction conditions

Choose the inverse response convention
$$
D_\rho(s)=M_{0\rho}^2-s-\Delta_\rho^{\mathrm{self}}(s).
$$
With $S=1+iT$, the two-pion cut obeys
$$
2\operatorname{Im}\Delta_\rho^{\mathrm{self}}(s+i0)
=\sum_{\pi\pi}\int d\Phi_2\,\overline{|\mathcal M|^2},
\qquad
\operatorname{Im}\Delta_\rho^{\mathrm{self}}(s+i0)
=+\sqrt{s}\,\Gamma_{\rho\to\pi\pi}(s).
$$
For the unit-projector convention this is
$$
\operatorname{Im}\Delta_\rho^{\mathrm{self}}(s+i0)
=\frac{G_V^2\sqrt{s}}{48\pi F^4}
(s-4M_\pi^2)^{3/2}\theta(s-4M_\pi^2).
$$
The upper-rim inverse response therefore has imaginary part
$-\sqrt{s}\,\Gamma_{\rho\to\pi\pi}(s)$. Its response-level spectral density is
$$
\rho_\rho^{\mathrm{vec}}(s)
=\frac{\operatorname{Im}\Delta_\rho^{\mathrm{self}}(s+i0)}
{\pi\left[
\left(M_{0\rho}^2-s-\operatorname{Re}\Delta_\rho^{\mathrm{self}}(s)\right)^2
+\left(\operatorname{Im}\Delta_\rho^{\mathrm{self}}(s+i0)\right)^2
\right]}\geq0
$$
above threshold. This is the density of the dressed response; it is not the
complete source-specific two-pion TFF density, which also contains
$F_\pi^V$ and the ALP production partial wave.

For $N$ subtractions at $s=0$, below the two-pion threshold, the real part is
$$
\operatorname{Re}\Delta_\rho^{\mathrm{self}}(s)
=\sum_{k=0}^{N-1}c_{\rho k}s^k
+\frac{s^N}{\pi}\,\mathrm{PV}\!\int_{4M_\pi^2}^{\infty}
\frac{\operatorname{Im}\Delta_\rho^{\mathrm{self}}(x+i0)}
{x^N(x-s)}\,dx.
$$
For the unregulated unit-projector cut, $\operatorname{Im}\Delta(x)\sim x^2$;
$N=2$ is logarithmically divergent and $N=3$ converges. A form factor or
ultraviolet completion can change this count. The three-subtraction
polynomial contains $c_{\rho0},c_{\rho1},c_{\rho2}$.

One real-axis mass and residue convention is
$$
M_{0\rho}^2=M_\rho^2+\operatorname{Re}\Delta_\rho^{\mathrm{self}}(M_\rho^2),
$$
This is paired with the on-shell real-part subtraction
$$
\operatorname{Re}\Delta_\rho^{\mathrm{self}}(s)
-\operatorname{Re}\Delta_\rho^{\mathrm{self}}(M_\rho^2)
-(s-M_\rho^2)\left.\frac{d}{ds}\operatorname{Re}\Delta_\rho^{\mathrm{self}}(s)\right|_{s=M_\rho^2}.
$$
The subtracted real part has zero value and slope at $M_\rho^2$, giving a
unit real-axis inverse-response slope. These two conditions do not fix the remaining
quadratic subtraction $c_{\rho,2}$. The complex pole must instead be found by
continuing $D_\rho^{-1}$ to the second sheet; its residue depends on
$\left.\frac{d}{ds}\Delta_\rho^{\mathrm{self,II}}(s)\right|_{s=s_p}$ and is
not fixed by the real-axis mass condition alone. The numerator of a dressed
response must respect whichever residue convention is used.

# Stable-pole limit and verification

With a positive regulator,
$$
\rho_\rho^{\mathrm{stable}}(s)
=\frac{\varepsilon}{\pi[(M_\rho^2-s)^2+\varepsilon^2]}
\xrightarrow{\varepsilon\to0^+}\delta(s-M_\rho^2),
\qquad \int_0^\infty\rho_\rho^{\mathrm{stable}}(s)\,ds\to1
$$
for $M_\rho^2>0$, provided local mass and residue terms are absorbed
consistently. The script prints the intermediate tensor contraction, spin
sum, phase space, width, optical-theorem sign, dispersion remainder,
ultraviolet subtraction count, real-axis conditions, and stable spectral
weight using SymPy's LaTeX printer.

Run:

```text
.venv/bin/python scripts/derive_alp_rho_pipi_self_energy.py
.venv/bin/python scripts/derive_alp_rho_pipi_self_energy.py --test
```

The five tests passed. They check the flavor trace and on-shell tensor-mass
normalization, P-wave threshold and positive cut, optical-theorem sign and
spectral-density form, Cauchy subtractions and real-axis conditions, and
stable-pole weight plus the differences among the three off-shell shapes.

# Limits and relation to the TFF

This derivation supplies conditional response-level width and subtraction
conventions; it chooses no numerical $G_V$, $F$, rho mass, or subtraction
coefficient. The remaining third subtraction and off-shell continuation
must be stated before a numerical mixed $F^{(\rho I)}$ term is used. Applying
one of these widths to the crossed $t/u$ poles, to $W_\gamma(y)$, or on top
of the elastic phase already represented by $F_\pi^V$ requires a separate
analytic justification. The source-specific overlap between this response
and the two-pion double-spectral reconstruction remains open; see the
[rho occurrence map](alp-rho-occurrence-map.md).
