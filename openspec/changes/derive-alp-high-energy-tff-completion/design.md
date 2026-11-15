# Design

## Context

See `proposal.md` for motivation. Available: the tree TFF `alp-doubly-virtual-tff.md`; the free-quark leading local OPE `alp-hard-current-ope.md`, which gives $F_{a\gamma^*\gamma^*}\supset \tfrac{2i}{(q^2+i0)m_a^2}\sum_fQ_f^2p^\sigma\langle0|\bar q_f\gamma_\sigma\gamma_5q_f|a(p)\rangle$ for nearly symmetric hard photons; and the gate record `alp-isoscalar-tff-input-gates.md`, which lists scale, scheme, thresholds, and low-energy separation as open. The calculation plan requires a declared matching scale, finite-$m_a$ terms, and no use of the $\eta/\eta'$ coefficient. The low-energy representation is still being assembled (`assemble-low-energy-doubly-virtual-tff-spectra`) and the vector fit is nonconverged.

## Goals / Non-Goals

**Goals:**
- Obtain a source-linear asymptotic TFF in each of the declared virtuality regimes with scheme and threshold provenance.
- Obtain a matching prescription tested on explicit low-energy forms that are independent of the fit.

**Non-Goals:**
- Numerical calibration against the vector fit, which waits for `alp-6os.1`.
- Higher-order QCD corrections beyond the order fixed in task 2, and the lepton-loop integration (`alp-6os.4`).
- Choosing a numerical ALP-source benchmark.

## Decisions

1. **Express matrix elements in the no-mixing source basis.** Relate $\langle0|\bar q_f\gamma_\sigma\gamma_5q_f|a\rangle$ to $U_a,Y_a$ and the contact coefficients via the LO ChPT source construction, split into isovector, octet-isoscalar, and singlet parts. Alternative, using physical-meson decay constants, was rejected because it would import $\eta/\eta'$ information.
2. **Handle the singlet axial current and gluonic operator explicitly.** The singlet current has an anomalous dimension and the anomaly links it to $G\tilde G$; the ALP–gluon coupling enters through this operator. Derive what survives at leading order and flag the scheme dependence, rather than assuming nonsinglet scaling for all flavors.
3. **Use the local OPE for the symmetric limit and the light-cone expansion for the asymmetric limit.** The singly virtual limit needs the nonlocal axial-current matrix element (an ALP light-cone distribution amplitude). Derive the leading-twist structure and the asymptotic distribution-amplitude assumption separately, labelled as an assumption; do not claim it from the local OPE.
4. **Finite $m_a$ via expansion in $m_a^2/Q^2$ of the hard kernel**, keeping $p^2=m_a^2$ exact in the tensor projection. Alternative, setting $m_a\to0$, is rejected because the scalar projection divides by $m_a^2$.
5. **Heavy-quark thresholds by decoupling.** Include $c$ and $b$ as active flavors above their thresholds with matching of the axial-current matrix elements, and also record the heavy-quark mass-suppressed contribution below threshold. Check against the supported range $2m_e<m_a<2m_\tau$.
6. **Matching in the spectral representation.** Following Messerli:2025rnv, construct a subtracted asymptotic double-spectral contribution above a declared $s_{\rm m}$ with the independently derived ALP coefficient. Use effective-pole variants to supply residual normalization and singly virtual constraints when their parameters admit a solution. Preserve the assigned low-energy discontinuities where the added spectral support is disjoint; elsewhere require an explicit overlap subtraction or state the conditional model and its unverified physical overlap. Treat sharp switches and interpolations as Euclidean diagnostics unless an analytic continuation is established; they are not physical-cut prescriptions. Test the analytic construction on symbolic low-energy forms and quantify matching-scale and completion dependence analytically in terms of symbolic parameters, separately from hadronic and scheme inputs.
7. **Pair every derivation with a script and a memory record**, per `AGENTS.md`; scripts are written by the `coder` agent, SymPy-based, with a `--test` invocation as in `derive_alp_hard_current_ope.py`.

## Risks / Trade-offs

- [The singlet and gluonic sector is scheme dependent and may not be fixed by the low-energy sources] → State the scheme, derive only what the sources fix, and list the remainder as a named input.
- [The singly virtual limit depends on a distribution amplitude that is not derivable here] → Provide the structure with the shape as a declared assumption and vary it as a model uncertainty.
- [Literature results for $\eta/\eta'$ may be mistaken for ALP results] → Use them only as cross-checks of the method, with the source dependence made explicit.
- [Matching to a fit-calibrated form is blocked] → Test on benchmark symbolic forms and label the result conditional.

## Open Questions

- Which order in $\alpha_s$ is required for the Wilson coefficient? Default: leading order, with the NLO correction located in the literature and quoted as an uncertainty estimate if not derived.
- Do the extra derivative and twist-4 operators matter at the accuracy needed for the loop integral? Default: bound them as a model uncertainty.
