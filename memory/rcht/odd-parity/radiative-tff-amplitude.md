---
title: Normalized neutral-pseudoscalar transition amplitude in tensor RChT
status: derived
created: 2026-10-13
updated: 2026-10-14
bibtex_keys: [Ruiz-Femenia:2003jdx, Ecker:1988te, Ecker:1989yg, CLEO:1997fho, BaBar:2009rrj]
provenance: analytical-derivation
---

# Normalized neutral-pseudoscalar transition amplitude

## Inputs and assumptions

The VJP basis, WZW action and their normalizations are the literature inputs in
[the operator digest](../../extracts/Ruiz-Femenia:2003jdx/vector-odd-parity.md)
and [the two-photon digest](../../extracts/Ruiz-Femenia:2003jdx/pseudoscalar-two-photon-form-factor.md).
The VVP basis is in [its digest](../../extracts/Ruiz-Femenia:2003jdx/vvp-three-pseudoscalar.md).
Photon-vector mixing and the tensor quadratic action are in
[the diagonalization derivation](../foundations/rcht-photon-vector-diagonalization.md).
The neutral meson field and $Q$ follow [the flavor record](radiative-tff-flavor.md).

Assume tree level, one common vector mass $M_V$, common $F_V$, leading order
in the electromagnetic coupling and the single-trace operator action already
recorded in this repository. Neutral meson flavor matrices commute with $Q$.
The formulas also define an assumed single-trace extension to the existing
nonet fields; the cited octet paper does not establish singlet matching or
additional OZI operators. No meson matrix is rescaled here. A numerical
$\eta,\eta'$ application must supply its mass-spurion traces and field/decay
constant normalization explicitly; its physical mass alone is insufficient.
Natural units, masses in GeV, metric $(+---)$, source $\epsilon_{0123}=+1$.

Define $s_1=q_1^2$, $s_2=q_2^2$, $p^2=m_P^2$, $p+q_1+q_2=0$ for the
all-incoming vertex, and $D_i=M_V^2-s_i$. Set

$$
F_Q=\operatorname{tr}(\{Q,Q\}T_P),\qquad
F_H=\operatorname{tr}(\{Q,Q\}\{T_P,\chi\}),\qquad \chi=2B_0m_q.
$$

$F_Q$ is dimensionless and $F_H$ has dimension GeV$^2$; $F_H$ includes the
normalization of the $\chi=2B_0m_q$ building block.
For $\pi^0$ with $\chi_u=\chi_d=m_\pi^2$ at leading chiral order,
$F_Q=\sqrt2/3$, $F_H=2m_\pi^2F_Q$. Treating this relation as a universal
identity for the physical $\eta'$ would omit its anomalous mass contribution.

## Derivation

Reproduce the actual contractions and reductions with
[`scripts/derive_radiative_tff_amplitude.py`](../../../scripts/derive_radiative_tff_amplitude.py):

```text
.venv/bin/python scripts/derive_radiative_tff_amplitude.py
```

1. Restrict the quadratic tensor action to the exact electromagnetic tensor
   subspace $V=rF_\gamma$. After removing the factor $F_\gamma^2/4$, it is
   $(M_V^2-s)r^2-2g_\gamma r$. Differentiation gives
   $r=g_\gamma/(M_V^2-s)$ with $g_\gamma=\sqrt2 eF_VQ$ in flavor space.
   This obtains the internal response directly from the action; the external
   on-shell $1/M_V$ wavefunction is not used for an internal photon line.
2. Expand $u_\mu=-\sqrt2\partial_\mu\Phi/F$,
   $f_+^{\mu\nu}=-2eQF_\gamma^{\mu\nu}$,
   $i\chi_-=\sqrt2\{\Phi,\chi\}/F$.
3. Insert Fourier derivatives $\partial\to iq$ in each operator. Contract all
   Lorentz indices with the explicit metric and Levi-Civita tensor. The
   derivative and mass blocks have opposite building-block signs.
4. Add both photon assignments for VJP and both vector Wick assignments for
   VVP.
5. Expand the WZW trace with both photon assignments. Choose the global meson
   phase so the pion WZW form factor is positive. This flips every term
   relative to the negative WZW coefficient in the source amplitude and leaves
   all relative signs intact.

The signed form factor is defined by
$i\mathcal M=i e^2\epsilon_{\alpha\beta\rho\sigma}
\epsilon_1^\alpha\epsilon_2^\beta q_1^\rho q_2^\sigma
F_P(s_1,s_2)$, with the chosen global phase.

For source-level ALP insertion, the
[shared photon/action phase dictionary](../../alp/photon/alp-photon-odd-action-phase-matching.md)
now gives $E=\epsilon(q_1,\varepsilon_1,q_2,\varepsilon_2)=-R$ with
$R=\epsilon(\varepsilon_1,\varepsilon_2,q_1,q_2)$ used above. In the
unrephased source basis the action is $i e^2E F_P$. An ALP LSZ overlap
computed for that source field changes sign when used with the rephased
meson amplitude $i e^2R F_P$. This flips every term together and fixes
the relative anomaly–resonance phase; the scalar formula below is unchanged.

$$
\begin{aligned}
C(s,t)&=F_Q\big[(c_1+c_2-c_5)m_P^2
 +(-c_1+c_2+c_5-2c_6)s+(c_1-c_2+c_5)t\big]+4c_3F_H,\\
H(s,t)&=F_Q\big[(d_1-d_3)m_P^2+d_3(s+t)\big]+4d_2F_H.
\end{aligned}
$$

The executed script yields

$$
\boxed{
F_P(s_1,s_2)=\frac{\sqrt2 N_CF_Q}{8\pi^2F}
 +\frac{4F_V}{FM_V}\left[\frac{C(s_1,s_2)}{D_1}
                         +\frac{C(s_2,s_1)}{D_2}\right]
 -\frac{4\sqrt2 F_V^2}{F}\frac{H(s_1,s_2)}{D_1D_2}.}
$$

Every term has dimension GeV$^{-1}$. The formula is symmetric under
$s_1\leftrightarrow s_2$; the spacelike observable is
$Q^2|F_P(-Q^2,0)|$ in GeV. The antisymmetric Lorentz structure obeys both
electromagnetic Ward identities. Widths, loop effects and resonance mass
splittings are outside this derivation.

## Published normalization and rank consequences

For the pion with both photons real, the script obtains

$$
F^{VJP}_\pi(0,0)=\frac{8\sqrt2 F_Vm_\pi^2}{3FM_V^3}
(c_1+c_2+8c_3-c_5),\qquad
F^{VVP}_\pi(0,0)=-\frac{8F_V^2m_\pi^2}{3FM_V^4}
(d_1+8d_2-d_3).
$$

These reproduce both amplitudes in Eq. (34) of Ruiz-Femenia:2003jdx with the
declared common phase flip. The integer 8 arises from the mass building block
and the two vector assignments, rather than from an octet dimension factor.
All five OPE substitutions give, for the explicitly assumed relation
$F_H=2m_P^2F_Q$,

$$
F_P(s_1,s_2)=\frac{\sqrt2 F_Q}{8\pi^2F}
\frac{N_CM_V^4-4\pi^2F^2(m_P^2+s_1+s_2)}{D_1D_2}.
$$

In particular $F_\pi(0,0)=1/(4\pi^2F)-Fm_\pi^2/(3M_V^4)$ at $N_C=3$,
which reproduces the published $(1-\Delta)$ correction. The script verifies
that the full TFF sensitivity annihilates **all six** homogeneous OPE null
directions in this limit: adding these TFF rows to the OPE block preserves
rank 5, and adding them to the action-derived common-pole radiative baseline
preserves rank 7.
With independent $F_H$, additional sensitivity is restricted to $c_3,d_2$;
$c_4,c_7,d_4$ and the common $c_2,c_5,c_6$ shift remain absent from neutral TFFs.
The latter is explicitly $\delta c_2=\delta c_5=\delta c_6$, with every
other coupling unchanged. The script also verifies those four invisible
directions for **arbitrary** $F_H$. Thus OPE plus neutral TFF data cannot
exceed rank 7 in this truncation. The equal-spurion rank-5 statement is a
stronger special case. The subsequent
[SM radiative action reconstruction](radiative-sm-amplitude.md) finds that
its common-pole radiative baseline spans only one correlated $c_3,d_2$
combination and also reaches $c_4$ through charged kaons. Adding TFFs with
the same **independent** exploratory $\chi$ inputs can then raise rank
$7\to8$ by the other $c_3,d_2$ combination. For equal-spurion TFF rows, the
zero-added-rank result remains exact.

The constrained single-virtual expression approaches a nonzero constant at
large $Q^2$, $\sqrt2 F_QF/(2M_V^2)$, under these assumptions. Thus this
single-resonance OPE approximation does not enforce the single-virtual
$1/Q^2$ asymptotic behavior. Connecting the data enables a diagnostic fit;
it establishes no precision validity over the CLEO/BaBar range.

## Implementation and data

[`pseudoscalar_tff.py`](../../../scripts/implementation/pseudoscalar_tff.py)
implements the full uneliminated amplitude. `TFFInputs` requires explicit
$F,F_V,M_V,m_P,F_Q,F_H$; no $\eta$ or $\eta'$ mass-spurion value is inferred.
The factories in
[`odd_parity_data.py`](../../../scripts/implementation/odd_parity_data.py)
connect selected CLEO mesons and the BaBar pion table to the existing
`Observation` interface. They retain the checksum-pinned data and authorized
likelihood policy: CLEO diagonal statistical plus systematic variance;
BaBar mean absolute asymmetric point error converted MeV to GeV, with a
separate 2.3% form-factor normalization nuisance. The 3% cross-section
normalization belongs to the cross-section observable. A multi-block fit must
declare block independence through the existing `FitProblem` interface.
No fit, coupling estimate, pooled experimental result or $\chi^2$ is reported.
