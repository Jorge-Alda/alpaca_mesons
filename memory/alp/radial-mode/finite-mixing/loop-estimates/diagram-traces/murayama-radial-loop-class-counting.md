---
title: Loop classes with the radial mode in the Murayama model, by counting
status: mixed
created: 2026-11-23
updated: 2026-11-23
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh]
provenance: assumption
---

# Loop classes with the radial mode, by counting


This record is task 3.1 of the OpenSpec change `phi-dependent-loop-tree-crossing-murayama`. It lists the loop classes of the estimate once the radial mode $\sigma$ is added to the ALP-only classes. It gives for each class the power of $\epsilon=f_\pi/f_a$ relative to the leading tree term, the scale of the loop, and the contraction left symbolic.

**No loop is evaluated here.** The powers are labels from end counting, used to organise the estimate. They are not an expansion of the crossing, which is implicit and exact in $\epsilon$ ([crossing record](../../../../loop-corrections/finite-mixing/power-counting-and-validity/crossing-and-validity/alp-loop-crossing-exact-tree-condition.md)). No loop integral, scalar loop function, loop coefficient or loop amplitude is computed.

## Inputs

- ALP-only classes 1 to 6, with their powers and contractions: [power-counting record, Result 4](../../../../decay-amplitudes/vertices-and-kinematics/tree-vertices/finite-mixing/alp-dressed-vertex-power-counting-invariants.md#result-4-loop-classes-of-the-amplitude).
- Radial couplings $\sigma PP$, $\sigma aP$ and $\sigma aa$, from the Lagrangian with $\widehat S_{10}$, $\widehat P_{11}$ and $\widehat S_{20}$: [only-up radial record](../../../couplings-and-decays/finite-mixing/murayama-sigma-two-pseudoscalar-decays.md) and [general-source record](../../../couplings-and-decays/finite-mixing/radial-decay-polar-cartesian-general-sources.md).
- Exact tree with radial exchange, which gives the $\sigma$-exchange denominators: [radial-exchange record](../../../../three-pseudoscalar/finite-mixing/alp-three-pseudoscalar-radial-exchange.md).
- Scale prescription per class: [naturalness record](../../../../loop-corrections/finite-mixing/power-counting-and-validity/loop-counting/alp-loop-naturalness-assumption.md).
- Radial mass $m_\sigma^2 = 2\lambda f_a^2 + (B_0f_\pi^2/f_a^2)(\langle\widehat S_{10}\rangle-2\langle\widehat S_{20}\rangle)$: [only-up radial record](../../../couplings-and-decays/finite-mixing/murayama-sigma-two-pseudoscalar-decays.md).

## Counting rule

1. The leading tree term carries one external ALP end, with power $\epsilon$ and the pole factor $1/\sqrt{N(m_a^2)}$ ([tree record](../../../../decay-amplitudes/vertices-and-kinematics/tree-vertices/finite-mixing/alp-finite-mixing-three-meson-tree.md)).
2. Each internal ALP-like end carries $\epsilon$. An internal ALP line carries two ends and the factor $1/N(m_A^2)$ from its pole normalization. This is the rule of the ALP-only classes.
3. Each internal $\sigma$ end carries $\epsilon$ as well, because the $\sigma PP$ coupling is $B_0\widehat S_{10}/f_a$ and $1/f_a=\epsilon/f_\pi$. Per end, $\sigma$ and the ALP are of the same order.
4. The radial mode has no admixture of mesons at leading order ($\widehat P_{10}=0$, $\langle\widehat S_{11}\rangle=0$), so its propagator carries no $N$ factor. $\sigma$ lines therefore give $\epsilon^2$ per line, without $1/N$.
5. Mixed couplings such as $\sigma aP$ carry one $\sigma$ end and one ALP end. They are counted as two ends.

## Classes

The classes are the cycles listed in [the estimate record](../../../couplings-and-decays/finite-mixing/murayama-radial-up-only-estimate-result.md). For each class the number of internal lines $r$ and the loop scale $m_d$ are as follows.

| class | internal lines | diagrams | scale $m_d$ |
|---|---|---|---|
| 1 | one ALP line, closed on one vertex | 1.1–1.3 | $m_a$ |
| 2 | one ALP line and one meson line | 2.1–2.5 | $m_a$ |
| 3 | two ALP lines | 3.1–3.3 | $m_a$ |
| 4 | one radial line and one meson line | 4.1–4.14 | $m_\sigma$ |
| 5 | one radial line and one ALP line | 5.1–5.12 | $m_a$ |
| 6 | two radial lines | 6.1–6.7 | $m_\sigma$ |
| 7 | ALP self-energy on an external leg | 7.1–7.4 | $m_a$ |
| 8 | radial self-energy on an external leg | 8.1–8.2 | $m_\sigma$ |

The diagrams are drawn in [the atlas](../../../../../../latex/murayama_loop_diagrams.tex), and their values are in [the estimate record](../../../couplings-and-decays/finite-mixing/murayama-radial-up-only-estimate-result.md).

Every diagram contributes $S_d(m_d^2)^{2-r_d}/(16\pi^2)$ times its vertex trace, as in [the estimate record](../../../couplings-and-decays/finite-mixing/murayama-radial-up-only-estimate-result.md). The residues of the external legs are in classes 7 and 8: class 7 is nonzero on the charged legs and zero on the $a$ and $\pi^0$ legs; class 8 is nonzero on the $a$ and $\pi^0$ legs and zero on the charged legs.

## Radial exchange in the tree

Tree-level $\sigma$ exchange has two $\sigma$ ends, one at the $\sigma aP$ vertex and one at the $\sigma PP$ vertex. Counting as above, it is of relative order $\epsilon^2$ with respect to the leading tree term. It is therefore of the same order as classes 2, 3, 5, 6, R2 and R3. It is kept exact in the tree denominator, as the crossing record requires, and is not counted as a loop class. This observation is relevant to the form of the crossing in task 4.1.

## Limits of the counting

- The powers are the leading end counting. They are not verified by evaluating any diagram, as for the ALP-only Result 4.
- Cancellations between classes of the same power are listed in [cancellation record](../../../../loop-corrections/finite-mixing/power-counting-and-validity/crossing-and-validity/alp-loop-estimate-cancellations.md) and are not reassessed here.
- The counting does not give the size of any contraction. The contractions are symbolic.
