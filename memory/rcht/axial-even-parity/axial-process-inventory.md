---
title: Axial-vector inventory for existing ALP processes
status: derived
created: 2026-11-10
updated: 2026-11-10
bibtex_keys: [ParticleDataGroup:2024cfk]
provenance: analytical-derivation
---

# Axial-vector inventory for existing ALP processes

This record lists, for every ALP process class in the [ALP index](../../alp/INDEX.md), whether an axial-vector state can appear as an external particle, and whether a single axial propagator can be exchanged. Each exclusion gives its reason. It is the input that tells the vertex work (tasks 3 and 4) which channels need an axial pole and which do not.

## Inputs and assumptions

- The ALP $a$ has $J^{PC}=0^{-+}$. Strong and electromagnetic subprocesses conserve $P$ and $C$; the ALP vertices are treated as sources with $a$ at rest carrying $J=0$.
- The axial states are the $1^{++}$ nonet: $a_1(1260)$, $f_1(1285)$, $f_1(1420)$ and the $1^{++}$ strange member $K_1$ (strange mixing is deferred, see [design open questions](../../../openspec/changes/define-axial-meson-rcht-sector/design.md)). The $1^{+-}$ nonet ($b_1$, $h_1$, $K_{1B}$) is outside this sector.
- Tree level, one axial propagator per diagram, even-sector interactions only. Loops are not inventoried.
- $J^{PC}$ assignments of the resonances follow [@ParticleDataGroup:2024cfk].

## Quantum-number rules

Let $P_f$ be the parity of a final state of $n$ particles with orbital angular momentum $L$. For two particles $P_f=P_1P_2(-1)^L$.

- **R1.** $0^{-+}\to 1^{++}+0^{-+}$ at $J=0$: $P_f=-(-1)^L$ must equal $-1$, so $L$ is even. Total spin $S=1$ (from the axial state) with $J=0$ needs $L=1$. Forbidden.
- **R2.** $1^{++}\to0^{-+}0^{-+}$: $P_f=(-1)^L$ must equal $+1$, so $L$ is even. $J=1$ with $S=0$ needs $L=1$. Forbidden.
- **R3.** $1^{++}\to 1^{--}+0^{-+}$ (for example $a_1\to\rho\pi$): $P_f=(-1)^L=+1$, so $L$ is even; $S=1$ from the vector, and $L=0,2$ give $J=1$. Allowed.
- **R4.** $1^{++}\to\gamma+0^{-+}$ (for example $a_1\to\pi\gamma$): same counting as R3. Allowed.
- **R5.** $1^{++}\to\gamma\gamma$ is forbidden for a single spin-1 state (Landau–Yang).
- **R6.** $C(1^{++})=+1$ and $C(\gamma)=-1$: no single-photon transition $A\leftrightarrow\gamma$, and no $A\to\ell^+\ell^-$ through one virtual photon.
- **R7.** $1^{++}\to 3P$ (for example $a_1\to3\pi$) is allowed by parity for odd orbital $L$ among the three pseudoscalars.

## Classification

| ALP process class (index) | Final state | External axial state | Single axial exchange | Classification |
|---|---|---|---|---|
| Two pseudoscalars ([index](../../alp/two-pseudoscalar/INDEX.md)) | $PP$ | $a\to AP$ is closed by R1 | $A\to PP$ is closed by R2 | Excluded (both). |
| Three pseudoscalars ([index](../../alp/three-pseudoscalar/INDEX.md)) | $3P$ | $a\to AP$ is closed by R1; $a\to A+2P$ is a four-body channel with no single-particle pole | Topology $a\to A^*P$, $A^*\to PP$ is closed by R2; no other one-$A$ tree topology | Excluded at single-exchange level. Contact-type axial terms are checked in task 3.1. |
| Photon plus two pseudoscalars ([index](../../alp/photon-two-pseudoscalar/INDEX.md)) | $\gamma PP$ | Not a two-body channel; not applicable | $a\to A^*P$, $A^*\to\gamma P$ is allowed by R4 | Included as an exchange channel. The $A\gamma P$ vertex is fixed in task 4.1. |
| Vector plus two pseudoscalars ([index](../../alp/vector-two-pseudoscalar/INDEX.md)) | $VPP$ | Not applicable | $a\to A^*P$, $A^*\to VP$ is allowed by R3 | Included. This is the leading axial exchange channel for $a\to VPP$. |
| Four pseudoscalars ([index](../../alp/four-pseudoscalar/INDEX.md)) | $4P$ | Not applicable | $a\to A^*P$, $A^*\to3P$ is allowed by R7 | Included as an exchange channel. |
| Baryons ([index](../../alp/baryon/INDEX.md)) | $B\bar B$ | Not applicable | The one-ALP baryon coupling is a local vertex without a meson propagator; an axial pole needs a baryon–axial vertex beyond the present order | Not included at the present order. Revisit only if a baryon–axial vertex is introduced. |
| Photon ([index](../../alp/photon/INDEX.md)) | $\gamma\gamma$ | $a\to A\gamma$ is possible only above the axial threshold | $a\to A^*\gamma$, $A^*\to\gamma$ is closed by R6 | Excluded for the two-photon channel. The vector-conversion contributions belong to the vector records. |
| Lepton pairs ([index](../../alp/lepton-pair/INDEX.md)) | $\ell^+\ell^-$ | $A\to\ell^+\ell^-$ is closed by R6 | No single axial pole in the $\gamma^*\gamma^*$ tree topology | Excluded at tree level. Hadronic intermediate axial states inside the $\gamma^*\gamma^*$ loop are not classified here. |
| Shared decay notes ([index](../../alp/shared/INDEX.md)) | not a process | not applicable | not applicable | Conventions only; no classification needed. |

## Consequences

- Axial exchange enters the ALP program at $a\to\gamma PP$, $a\to VPP$ and $a\to4P$. These are the channels for which tasks 4.1–4.2 must produce vertices.
- The two-pseudoscalar and three-pseudoscalar classes receive no single-axial pole. Any axial effect there must come from contact-type operators, which task 3.1 enumerates.
- External axial states are kinematically closed unless $m_a$ exceeds the axial mass threshold. The ALP mass range used by each process record decides whether an external axial state is needed. This is recorded here as a condition, not a result.

## Verification

- Coverage check: each of the nine entries of the [ALP index](../../alp/INDEX.md) appears in the table above, and each has a link to its subfolder index. The check was done by reading the index files; no script was run.
- The rules R1–R7 were derived by hand from $P_f$ and $C$ counting. They are not yet checked by a script. Their use is limited to the classification above.
