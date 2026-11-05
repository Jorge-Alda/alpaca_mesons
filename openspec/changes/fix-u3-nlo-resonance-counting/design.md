# Design

## Context

See proposal.md (Why). The relevant existing records are:

- [lagrangian.md](../../../memory/nlo/lagrangian.md): the $\delta$ counting of Kaiser:1998ds, with $\mathcal L^{(1)}$ ($L_5,L_8,\Lambda_1,\Lambda_2$) at $\mathcal O(\delta)$ and $\mathcal L^{(2)}$ ($L_4,L_6,L_7,L_{18},L_{25}$) at $\mathcal O(\delta^2)$, in the reduced $V_4=0$ source convention.
- [one-loop-effects-and-renormalized-lecs.md](../../../memory/nlo/one-loop-effects-and-renormalized-lecs.md): the $U(3)$ coefficients $\Gamma_i$, the $Z_A$ behaviour of the singlet constants, and $\Gamma_{47}=0$ under the truncation of Herrera-Siklody:1996tqr. It describes loops qualitatively but does not say which graphs are included.
- [resonance-saturation-and-dynamical-masses.md](../../../memory/nlo/resonance-saturation-and-dynamical-masses.md): the $SU(3)$ tree-level saturation of Ecker:1988te and the coupled Goldstone–pseudoscalar-resonance pole system. It states that $\Lambda_{1,2}$, $L_{18}$ and $L_{25}$ are not fixed by standard saturation.
- [rcht/foundations/](../../../memory/rcht/foundations/INDEX.md): the vector nonet in the antisymmetric-tensor representation and the width prescriptions. Axial-vector, scalar and pseudoscalar resonance nonets have no LO record yet.
- Extracts available: Kaiser:2000gs (the counting of the $U(1)_A$ anomaly and the topological source), Herrera-Siklody:1996tqr, Rosell:2004mn and Pich:2010sm (resonance loops).

All conventions and notation follow Alda:2026nrd. Any new symbol needs the user's approval before it is used.

## Goals / Non-Goals

**Goals:**
- One set of records that alp-gii.2, alp-gii.3 and alp-gii.8 can cite for "the declared order".
- A script that encodes the order assignments and the scale running, so a later change to the counting fails the check visibly.

**Non-Goals:**
- Evaluating any loop integral or self-energy (alp-gii.2).
- Coefficient-by-coefficient resonance matching (alp-gii.3).
- Numerical values of constants or resonance couplings (alp-gii.5).
- Resonance loops. These are excluded by decision, and only their size enters the uncertainty.

## Decisions

1. **$\delta$ counting of Kaiser:1998ds, not plain $1/N_c$ RChT counting.** The existing NLO records already use it, and it is the counting in which $\eta_0$ and the topological source are treated consistently. *Alternative:* RChT counting at leading order in $1/N_c$ with chiral orders only. It does not place $\Lambda_{1,2}$ and the $\eta'$ mass consistently and would require rewriting the existing records.

2. **Resonances at tree level with unexpanded propagators.** Couplings scale as $\sqrt{N_c}$ and masses as $N_c^0$. At $\sqrt s\lesssim1\,$GeV, $s/M_R^2$ is not small, so exchange graphs keep the full propagator. Their low-energy limit is $\mathcal O(\delta)$ ($N_cp^4$), which is how it is identified with the corresponding part of $L_i$ under Decision 6. *Alternative:* integrate the resonances out. This would reduce the project to $U(3)$ ChPT, which is inadequate near $1\,$GeV.

3. **Spin-one fields in the antisymmetric-tensor representation (Ecker:1988te, Ecker:1989yg)** for both vector and axial-vector nonets. This is consistent with [rcht-basics.md](../../../memory/rcht/foundations/rcht-basics.md). It avoids bilinear $A$–$\partial P$ mixing and gives the $\mathcal O(p^4)$ saturation without extra local terms. *Alternative:* the Proca or hidden-local-symmetry representation. It needs extra local terms to agree with QCD short-distance behaviour, and it would conflict with the existing vector records.

4. **Loop content: Goldstone-only one loops with LO $U(3)$ vertices (user decision).** These graphs are $\mathcal O(\delta^2)$. They include LO ALP-source vertices at first order in $1/f_a$. Loops with an internal resonance line are $\mathcal O(1/N_c)$ without chiral suppression. They are excluded, and the size quoted from Rosell:2004mn and Pich:2010sm enters the uncertainty. *Correction after a check of the sources (task 3.1):* these loops already renormalize the $\mathcal O(p^2)$ constants $F^2$ and $F^2B_0$ (Rosell:2004mn, Rosell:2005ai). They therefore start at $\delta^1$, below the declared order. The exclusion holds below the $R\pi$ thresholds only if $F$, $B_0$, $\tau$, $\Lambda_{1,2}$ and $L_i^r(\mu)$ are taken from data, not from large-$N_c$ or saturation estimates. This condition applies to alp-gii.3 and alp-gii.5. Goldstone loops with one resonance vertex are excluded for the same reason. The executable classification in the script states this explicitly.

5. **Resonance widths.** The Goldstone-loop self-energy of a resonance is formally a $1/N_c$ correction. It is not part of the strict $\mathcal O(\delta^2)$ Goldstone/ALP result. Near a resonance pole it is included only through the existing [width prescription](../../../memory/rcht/foundations/vector-width-lineshape-prescriptions.md), labelled as a resummation outside the strict counting. Away from poles the narrow propagator is used.

6. **Scheme and scales.** Use $\lambda(\mu_\chi)$ of Gasser:1983yg with the $U(3)$ $\Gamma_i$ of Kaiser:1998ds: $\Gamma_{18}=-1/4$, $\Gamma_{25}=0$, and $\Gamma_{47}=0$ at the retained $1/N_c$ order. The default is $\mu_\chi=M_\rho=0.77\,$GeV, the usual scale for saturation estimates, with variation in $[0.5,1.0]\,$GeV. The QCD-scale dependence through $Z_A$ is kept separate from $\mu_\chi$. Singlet inputs are quoted together with the QCD scale given by their source. Resonance exchange is $\mu_\chi$ independent, so the residual constant carries all of the $\mu_\chi$ running. The saturation hypothesis "residual $\approx0$ at $\mu_\chi\approx M_R$" is labelled as an assumption. **Notation (approved by the user, task 1.2):** the residual constants are written $L_i^r(\mu)$, with $\mu$ the chiral renormalization scale (written $\mu_\chi$ in this change's artifacts when it must be distinguished from the QCD scale). In the framework with explicit resonances, $L_i^r(\mu)$ therefore denotes the renormalized local constant that remains after explicit resonance exchange is removed. It is not the full ChPT constant, which also contains the resonance-exchange part. Every record that uses $L_i^r(\mu)$ states which of the two meanings applies.

7. **Truncation uncertainty.** Estimate the omitted $\mathcal O(\delta^3)$ terms as the size of the $\mathcal O(\delta^2)$ correction to an observable times $\delta_{\rm est}(s)=\max(1/N_c,\ s/\Lambda_\chi^2,\ M_K^2/\Lambda_\chi^2)$, with $\Lambda_\chi=4\pi f_\pi$. Add in quadrature the size of the excluded resonance loops from the extracts. The prescription is implemented as a function in the script. *Alternative:* vary $\mu_\chi$ only. This underestimates the error, because the scale-independent $\mathcal O(\delta^3)$ terms are not probed.

8. **Record layout.** Create a topic folder `memory/nlo/counting/` with an `INDEX.md` and atomic records: field content; power counting; loop classes; subtraction scheme and scales; counterterm basis; separation of resonance exchange and residual constants; validity range and truncation. Update `memory/nlo/INDEX.md` with a topic-level description. Add a cross-link from the existing loop and saturation records instead of duplicating their content.

9. **Delegation.** The `coder` agent writes `scripts/derive_nlo_u3_resonance_counting.py`, which uses sympy for the $\Gamma_i$ running and exact rational orders. Any new INSPIRE record goes through `bibliography-manager`. No `latex/` edits.

## Risks / Trade-offs

- [Near $1\,$GeV, $\delta_{\rm est}\sim0.7$, so convergence is poor] → State this in the validity record. Quote the uncertainty as a function of $s$ rather than one number.
- [Kaiser:1998ds and Herrera-Siklody:1996tqr use different source conventions ($V_4=0$ versus $W_i$)] → The scheme record states the translation, and the script checks the $\Gamma_i$ in both.
- [The resonance couplings to the singlet and topological sector are unknown, so $\Lambda_{1,2}$, $L_{18}$ and $L_{25}$ have no saturation estimate] → Keep them as free residual constants and pass them to alp-gii.5.
- [Excluding resonance loops removes $\mu_\chi$-dependent pieces that a full RChT one-loop treatment would cancel] → No such cancellation is needed at the declared order, because the excluded loops have their own counterterms. State this explicitly.
- [The width prescription mixes orders near poles] → Label it as a resummation. Downstream issues report both the strict-order and the resummed results.

## Open Questions

- Whether $\mathcal O(\delta^2)$ ALP-source operators beyond $L_{47}$ are needed depends on the source specialization chosen in alp-gii.8. The counterterm table leaves a marked row for them.
