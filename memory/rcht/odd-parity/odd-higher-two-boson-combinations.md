---
title: Higher odd-action combinations in stripped ALP two-boson vertices
status: derived
created: 2026-10-14
updated: 2026-10-15
bibtex_keys: [Alda:2026nrd, Ruiz-Femenia:2003jdx]
provenance: analytical-derivation
---

Canonical tensor and odd-coupling conventions follow [the Levi-Civita audit](levi-civita-convention.md). The linked original scripts and saved artifacts use the explicitly identified reference backend; apply that record's signed coefficient and tensor adapter when interpreting them in the canonical convention.

# Inputs and assumptions

Use the [certified higher action basis](odd-higher-covariant-operator-count.md),
[physical independent-source ALP insertions](two-source-alp-rcht-odd-parity-extension.md),
[canonical tensor audit](alp-two-boson-blind-direction-audit.md), and
[neutral generator conventions](radiative-sm-flavor-algebra.md).
This executes the projection proposed in the
[operator relevance assessment](odd-operator-alp-relevance.md).

Tree level, one physical ALP and two boson legs, metric $(+---)$ and
$\epsilon^{0123}=+1, \epsilon_{0123}=-1$. Both boson tensors have canonical wedge form.
The scalar vacuum $X=\chi_+$ is constant. Its term linear in a pseudoscalar
vanishes at a real scalar vacuum; connection/source corrections are retained
in the full action but do not enter this cubic vertex through an extra $X$
insertion. No transverse polarization or on-shell mass substitution is used.

The derivative and mass/source weights remain independent:

$$
u_\mu|_a=-\frac{2\partial_\mu a}{f_a}\mathbb M_a^{-1}(C),\qquad
(i\chi_-)|_a=\frac{4B_0a}{f_a}
\{c_q^{AP}+\mathbb M_a^{-1}(C),m_q\}.
$$

The projection uses abstract flavour matrices first. The second projection
uses the actual saved diagonal $\rho,\omega,\phi$ generator span and
$Q=\operatorname{diag}(2,-1,-1)/3$, with symbolic coefficients in that span.
$X=2\chi$, $\chi=2B_0m_q$, keeps the three quark masses symbolic; no partial
isospin replacement, fitted mixing angles or numerical rank tolerance is used.
Diagonal derivative and mass/source weights are independent in this stage.

# Executable calculation

Run [the script](../../../scripts/project_odd_higher_two_boson.py):

```text
.venv/bin/python scripts/project_odd_higher_two_boson.py
```

It reconstructs the actual epsilon/metric contractions of every selected
representative, replaces tensor legs by wedges and $u$'s Lorentz vector by
$p_a=-q_A-q_B$, and replaces ordered derivative momenta by their commuting
cubic-vertex values. Derivatives on constant $X$ vanish. For VVP it includes
both vector assignments with the matching symmetric/antisymmetric flavour
parity. Flavour words remain abstract in this first stage.

Let

$$
E=\epsilon_{\mu\nu\rho\sigma}q_A^\mu\varepsilon_A^\nu
q_B^\rho\varepsilon_B^\sigma,\qquad
s=q_A^2,\quad t=q_B^2,\quad h=q_A\cdot q_B,
\qquad m_a^2=s+t+2h.
$$

SymPy coefficient reduction expresses every retained image exactly as
$E$ times a linear polynomial in $s,t,h$ with the specified flavour traces.
It then reduces the joint trace/invariant columns by exact rref. Each saved
combination has the form $K_i=\sum_j R_{ij}g_j$, with $R$ computed by rref,
so $\sum_jg_j\mathcal P_j=\sum_i K_i\mathcal P_{p_i}$.
Both coefficient maps and every parent image are exported.

# Results

| Block | Action representatives | Abstract trace/invariant rank | Neutral rank |
|---|---:|---:|---:|
| $X$:VJP | 89 | 30 | 21 |
| $X$:VVP | 43 | 17 | 14 |
| $Ydd$:VJP | 31 | 12 | 9 |
| $Ydd$:VVP | 18 | 7 | 7 |
| Total above | 181 | 66 | 51 |

The abstract projection has 101 individually zero images. Neutral restriction
has 110 individually zero images among those 181; other dependencies relate
nonzero images. VJP in the neutral restriction has $B=Q$. VVP retains two
arbitrary vectors in the neutral generator span, so internal photon response
matrices in that span remain representable; it is not replaced by $Q$ in
both slots.

For electromagnetic-compatible diagonal ALP sources, the script derives
$[U,Q]=[Y,Q]=0$. The linear scalar-photon $f_-$ insertion is proportional to
a charge commutator, so all seven $fmdd$ representatives vanish in this neutral
cubic projection, including their ordered derivatives. Thus the 188 higher
representatives have rank 51 under the stated neutral, stripped projection.
This does not apply to charged calibration channels or arbitrary
flavour-changing axial/pseudoscalar sources.

# Scope and remaining work

The [joint power-counted source/channel calculation](odd-power-counted-channel-maps.md)
now combines all leading and higher raw parents, adds Fourier/source/canonical
normalization and checks every action-relation projection. It shows why the
independent-source stripped map below cannot directly descend to the EOM
quotient: three leading relations need the LO-compatible source equation
with its topological contact. Its compatible primitive ranks and on-shell
channel limits supersede using 51 alone for a truncation count; complete
ALP predictor and data-fit interpretation remain separate.

These are formal polynomial vertex combinations, **not 51 fit parameters or
a completed physical ALP sensitivity rank**. They exclude Fourier/source,
electric-charge and canonical LSZ normalization from the displayed coefficient
maps. They retain symbolic masses and arbitrary independent diagonal ALP
weights rather than one chosen ALP model. The leading 40 representatives have
not yet been combined with these columns. Fixed EOM-induced topological
contacts, anomaly contributions and full source/response assembly must enter
the normalized predictor before interpreting any null direction physically.

A $Ydd$ insertion with fixed scalar source can overlap a leading derivative
momentum polynomial; a fixed vacuum spurion can make an $X$ dressing overlap
existing flavour combinations. Neither is removed using independent-source
EOM by this vertex calculation. Additional reductions must be made jointly,
with the relevant fixed source/mass assumptions stated.

For $a\to VV$, use the VVP kernel with canonical vector legs. For
$a\to V\gamma$, combine direct VJP and VVP with one photon response. For
$a\to\gamma\gamma$, include both VJP photon assignments and the VVP
kernel with two responses, plus the anomaly/contact sector. Exact response
coefficients and on-shell kinematics may further reduce ranks. Hard-limit
matching has not been solved by this projection.

# Artifacts and verification

- [Machine-readable images and combinations](../../../data/vector_sector/odd_higher_two_boson_projection.json).
- [Step-by-step LaTeX contractions and coefficient maps](../../../data/vector_sector/odd_higher_two_boson_projection.md).

Every epsilon polynomial is reconstructed from its three invariant
coefficients. Both abstract and neutral matrices satisfy exact reconstruction
$M_{\rm pivot}R=M$. These checks use no sampled kinematics or floating-point
flavour ranks.
