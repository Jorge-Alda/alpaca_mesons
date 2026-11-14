# Tasks

Scripts are written by the `coder` agent. Each memory record has the front matter required by AGENTS.md, links its inputs, and is added to its `INDEX.md` and every parent index. Track with `bd` (alp-sjp.3, LO part only). Do not close alp-sjp.3. Do not commit unless asked. All work is analytical and symbolic; no numerical phenomenological inputs are needed.

## 1. Common amplitude construction

- [x] 1.1 Script: joint expansion of $\mathcal L_{\rm LO}$ in $(P_i,a)$ to a chosen total order on the neutral basis (plus charged fields where needed), with the ALP entering through general symbolic $M(x)$, $C_{R,L}(x)$, $\Theta(x)$ (Taylor coefficients kept as symbols), producing symmetrized amputated vertex tensors. Verify that the linear-spurion special case reproduces `derive_alp_phi2_lo_chipt.py` and `derive_alp_phi3_lo_chipt.py`, and that the quadratic blocks agree with `solve_alp_meson_spectrum.py`; if the archived blocks assumed linear sources, derive the general blocks and report any difference to the LO spectrum.
- [x] 1.2 Script: pole-projected amplitude evaluator contracting vertices with $v^{(n)}$ for parent and daughters and with $D^{-1}(s)$ for exchanges, using the spectrum solver. Verify on-shell equality of the full propagator with the pole sum and invariance under a random constant $T$ of the on-shell amplitude.
- [x] 1.3 Script: transformed-vertex route (pole-aligned exact chiral map with all nonlinear terms) and comparison with 1.2 at random kinematic points for the two-meson and three-meson vertices; report the separate size of the nonlinear contact terms.
- [x] 1.4 Memory record: the construction (inputs linked, assumptions on the $a$-dependence of the sources, route 1 versus route 2, validity limits), with the minimal invocation of each script verified.

## 2. Channels

- [x] 2.1 One and two pseudoscalars ($a\to P$ mixing, generic $a\to P_iP_j$ including flavour-changing and $K^0$-related rules): generic amplitudes with daughter projections; reproduce `alp-charged-u3-vector-redefinition.md`; weak-mixing and Bose checks; record in `memory/alp/two-pseudoscalar/`.
- [x] 2.2 Three pseudoscalars: derive the generic $a\to P_iP_jP_k$ amplitude (indices open, outgoing ALP admixture included, symmetry factors for identical daughters), then check by substitution that it reproduces `alp-finite-mixing-three-pion-amplitudes.md` for $3\pi^0$ and $\pi^0\pi^+\pi^-$; no channel-by-channel derivations; record in `memory/alp/three-pseudoscalar/`.
- [x] 2.3 Photon channels at LO: $a\to\gamma\gamma$ summed over meson components of the pole, and $a\to\gamma PP$ from the WZW sector with mixing; Ward identity check with the photon momentum; reproduce the first-order coefficient of Alda:2026nrd in the weak-mixing limit; record in `memory/alp/photon/` and `memory/alp/photon-two-pseudoscalar/`.
- [ ] 2.4 Baryon-antibaryon: amplitude structure from the LO CCWZ vertex with the physical-pole projection, threshold test in terms of the physical parent mass, reproduction of the existing baryon couplings; record in `memory/alp/baryon/`.

## 3. Inventory and rules

- [ ] 3.1 Script and record: exact selection rules from the pole vectors (charge, CP, flavour, exact $K^0$ pole without ALP component, channels that open or close at finite mixing); verify on random CP-conserving and CP-violating sources with the solver and on the unmixed limit; record in `memory/alp/shared/`.
- [ ] 3.2 Record: process inventory with LO status per channel (derived here / existing link / deferred with reason), thresholds as analytical inequalities in the physical masses (no numerical inputs), and the stable-versus-broad daughter classification.
- [ ] 3.3 Update every affected `INDEX.md` and parent index; run all new and related scripts and their tests together; check relative links and front matter.

## 4. Closing

- [ ] 4.1 Check each alp-sjp.3 acceptance criterion against the evidence for the LO part. Add a note to alp-sjp.3 listing what is done and what remains (resonance exchange and vector-meson channels, $\ell^+\ell^-$, NLO terms and widths), create a follow-up bead for the beyond-LO part linked to alp-sjp.3 and alp-sjp.7 as appropriate, and keep alp-sjp.3 open. Report changed files and a proposed commit message.
