---
title: Normalized leading single/double-trace scalar radiative and photon kernels
status: derived
created: 2026-10-14
updated: 2026-10-14
bibtex_keys: [Alda:2026nrd, Ruiz-Femenia:2003jdx, Ecker:1988te, Ecker:1989yg]
provenance: analytical-derivation
---

# First stage-B normalized kernels

## Inputs and scope

Use [stage-4 normalization and response](vector-stage4-amplitudes.md),
[trace enumeration](odd-breaking-flavor-enumeration.md),
[Lorentz projections](odd-breaking-lorentz-reduction.md), and the
[source/EOM audit](odd-trace-source-eom-audit.md).
Retain the leading derivative trace families, the leading $i\chi_-$ mass
trace families and the original $c_4$ commutator. No additional $\chi_+$
derivative dressing, two-derivative mass/source extension, loop, OPE relation
or independently fitted topological function is included in these kernels.
The existing vector action/source response remains unchanged.

These are normalized one-scalar-SM-pseudoscalar kernels, not a full independent
axial/pseudoscalar-source functional or a completed normalized ALP predictor.
Use actual external-state conjugations for charged flavor matrices, the
positive-strange phi phase, $\epsilon_{0123}=+1$, $(+---)$ and GeV units.

## Trace conventions and covariant representatives

For any flavor matrices $A,B,Z$, use the stage-A functions

$$
F_0=\tfrac12\operatorname{tr}((AB+BA)Z),\quad
F_1=\operatorname{tr}A\operatorname{tr}(BZ),\quad
F_2=\operatorname{tr}B\operatorname{tr}(AZ),\quad
F_3=\operatorname{tr}(AB)\operatorname{tr}Z.
$$

For VJP, replace the original trace in each $c_1,c_2,c_3,c_5,c_6,c_7$
covariant seed by $q^J_\ell=2F_\ell$, with four independent coefficient
families. The original coefficient family is $\ell=0$. There is no
additional double-trace $c_4$ because its C-odd projection vanishes.

For symmetric VVP, replace the trace in the $d_1,d_2,d_3,d_4$ seeds by
$(q^V_0,q^V_1,q^V_2)=(2F_0,F_1+F_2,2F_3)$ and symmetrize the two
vector assignments. Define $q_-=F_1-F_2$. A concrete covariant representative
of the scalar-visible antisymmetric family is

$$
\mathcal O_-=
\epsilon_{\mu\nu\rho\sigma}\left[
\operatorname{tr}(\nabla_\alpha V^{\mu\alpha})
\operatorname{tr}(V^{\rho\sigma}u^\nu)
-\operatorname{tr}((\nabla_\alpha V^{\mu\alpha})u^\nu)
\operatorname{tr}(V^{\rho\sigma})\right].
$$

Its real dimensionless coefficient is $d_-$. Both field assignments give
$q_-$ multiplying the antisymmetric Lorentz combination. The other two
source-resolved antisymmetric VVP directions lie in the scalar projection's
kernel; this scalar API does not claim to implement their axial-source
vertices. Connection terms and nonlinear expansions of these representatives
remain in the covariant action; they are not discarded by EOM conversion.

In the SM expansion set $Z=T_P$ for derivative traces and
$Z=\{T_P,\chi\}$ for mass traces, defining $h^J_\ell$ and
$h^V_\ell$ by the same functions. No neutral mass is substituted for an
independent mass-spurion trace.

## Executable contraction and result

Run [the source/contraction script](../../../scripts/audit_odd_trace_sources.py).
It reuses the actual Fourier-index contractions of the original operators,
weights their trace placements, and retains both vector/photon assignments.
The output is derived by SymPy reductions, not inserted from a guessed flavor
rescaling. It gives, separately for each symmetric trace family,

$$
C_\ell=q^J_\ell[(c_{1\ell}+c_{2\ell}-c_{5\ell})m_P^2
+(-c_{1\ell}+c_{2\ell}+c_{5\ell}-2c_{6\ell})s
+(c_{1\ell}-c_{2\ell}+c_{5\ell})t]+4c_{3\ell}h^J_\ell,
$$

$$
H_\ell=q^V_\ell[(d_{1\ell}-d_{3\ell})m_P^2+d_{3\ell}(s+t)]
+4d_{2\ell}h^V_\ell,
\qquad H_-=-d_-q_-(s-t).
$$

Write $C=\sum_{\ell=0}^3C_\ell$ and
$H=\sum_{\ell=0}^2H_\ell+H_-$ and insert them into exactly the existing
stage-4 $W(s)$ formulas:

$$
F_P(s_1,s_2)=\frac{\sqrt2N_Cq(Q,Q)}{8\pi^2F}
+\frac4{FM_V}[C(W_1,Q;s_1,s_2)+C(W_2,Q;s_2,s_1)]
-\frac{4\sqrt2}{F}H(W_1,W_2;s_1,s_2).
$$

For the radiative scalar coupling, with $F_K=\operatorname{tr}(T_V[[T_P,Q],\chi])$,

$$
G=\frac{2\sqrt2e}{FM_{\rm ext}M_V}
[C(T_V,Q;M_{\rm ext}^2,0)-4c_{40}F_K]
-\frac{4e}{FM_{\rm ext}}H(T_V,W_0;M_{\rm ext}^2,0).
$$

The external $1/M_{\rm ext}$ is canonical tensor normalization, not an
internal pole. The WZW coefficient and the established crossed/identical-
photon width factors are unchanged. The anti family is Bose-even because
both $q_-$ and $s-t$ change sign; it vanishes for both photons real.

## Implementation and checks

[`odd_trace_amplitudes.py`](../../../scripts/implementation/odd_trace_amplitudes.py)
supplies validated immutable `OddTraceCouplings`, trace functions, $C,H$,
normalized TFF/radiative components and widths. `from_baseline` sets all new
families to zero and recovers the existing kernel. Its stored $c_7,d_4$ slots
preserve baseline correspondence but are scalar-blind; double-trace $c_4$
values are rejected. No fit/data factory is silently switched to the new model.

Run [the verification runner](../../../scripts/verify_odd_trace_amplitudes.py):

```text
.venv/bin/python scripts/verify_odd_trace_amplitudes.py
```

All 32 new/established tests pass: every baseline coefficient, charged $c_4$,
double-trace normalization, photon singleton-trace blindness, anti-family
sign, Bose exchange, crossing, internal phase invariance, invalid inputs
and poles, plus the existing vector/amplitude/data regressions. The symbolic
script checks the unreduced photon Ward contractions before projection.

## What is counted now

The leading catalogue has 40 source-resolved tensors, 29 beyond the original
eleven, before a complete EOM/source-contact quotient. This scalar API stores
38 coefficient slots, omitting two anti source-only directions; the stored
slots include known scalar/photon-blind directions and are not 38 fit parameters.

The script constructs the exact formal scalar EM polynomial design matrix,
using $\operatorname{tr}Q=0$ but no physical-mass EOM substitution. Its rank
is 23, compared with baseline rank 8: **15 additional scalar electromagnetic
kernel directions**. This counts generic flavor/virtuality structures with
independent mass traces, not full-action operators, actual channels with
fixed central inputs or a measurement-weighted/marginalized fit rank.

The full new action-operator count will follow the covariant/source closure
and EOM audit described in [the source record](odd-trace-source-eom-audit.md),
before calibration. One-spurion normalized kernels and complete ALP/source
predictors remain to be derived; a zero-spurion count is not the final count
for the full selected extension.
