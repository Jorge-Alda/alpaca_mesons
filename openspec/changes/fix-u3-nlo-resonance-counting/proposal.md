# Proposal

## Why

The epic alp-gii (beads issue alp-gii.1, GitHub #31) needs one declared expansion before the loop self-energies (alp-gii.2), the resonance-saturation matching (alp-gii.3) and the loop corrections to ALP decay vertices (alp-gii.8) can be calculated. These three issues are blocked by alp-gii.1. The current `memory/nlo/` records contain the local $U(3)$ NLO basis of Kaiser:1998ds, a qualitative note on loops and the $U(3)$ counterterm coefficients, and a tree-level note on resonance saturation. None of them fixes, in one place, which fields, graphs and counterterms belong to the calculated order once explicit resonances are kept. Without that, the later issues can count a contribution twice or leave scale dependence uncancelled.

## What Changes

- Declare the field content: the pseudoscalar nonet including $\eta_0$, the ALP through the sources of Alda:2026nrd, and explicit nonets of scalar, pseudoscalar, vector and axial-vector resonances.
- Declare the combined counting of Kaiser:1998ds, $p^2\sim m_q\sim 1/N_c\sim\delta$, extended to the resonance sector. Resonances are tree-level fields whose propagators are not expanded in $p^2/M_R^2$. The calculated order is $\mathcal O(\delta^2)$ in the Goldstone/ALP sector.
- Fix the included loop classes. One-loop graphs with only nonet Goldstone internal lines and LO $U(3)$ vertices, including LO ALP-source vertices, are included. Loops with internal resonance lines are **excluded** and enter only the truncation uncertainty. This choice was confirmed by the user.
- Fix the subtraction scheme: the $\overline{\rm MS}-1$ convention of Gasser:1983yg with the $U(3)$ coefficients $\Gamma_i$ of Kaiser:1998ds and Herrera-Siklody:1996tqr. Declare the chiral renormalization scale $\mu_\chi$ and, separately, the QCD scale dependence of the singlet constants through $Z_A$.
- List the renormalized operator basis: for each local constant ($L_4$–$L_8$, $L_{18}$, $L_{25}$, $\Lambda_1$, $\Lambda_2$, $L_{47}$ and the ALP-source operators), give its $\delta$ order, its chiral counterterm coefficient, its QCD-scale behaviour, and whether resonance exchange contributes to it.
- Fix the rule against double counting: after explicit resonance exchange is removed, the remaining local constant $L_i^r(\mu_\chi)$ is a residual constant. The matching itself is left to alp-gii.3.
- Declare the validity range: external momenta and ALP masses up to about $1\,$GeV. State the truncation-uncertainty estimate and the treatment of resonance widths near poles.
- Add an executable script that tabulates the $\delta$ order of operators and graph classes and checks the scale independence of renormalized combinations. Add step-by-step memory records that explain the script.

## Capabilities

### New Capabilities

- `nlo-u3-resonance-counting`: The declared power counting, field content, included/excluded loop classes, subtraction scheme, renormalized counterterm basis, singlet/topological treatment, double-counting rule and validity range of the combined $U(3)$ NLO plus explicit-resonance theory used by all alp-gii calculations.

### Modified Capabilities

None. There are no existing specs under `openspec/specs/`.

## Impact

- New memory records under `memory/nlo/`, with `memory/nlo/INDEX.md` and `memory/README.md` updated. [one-loop-effects-and-renormalized-lecs.md](../../../memory/nlo/one-loop-effects-and-renormalized-lecs.md) and [resonance-saturation-and-dynamical-masses.md](../../../memory/nlo/resonance-saturation-and-dynamical-masses.md) are linked from the new records and corrected only if they conflict with the declared scheme.
- One new script under `scripts/`, written by the `coder` agent.
- Possible additions to `bibliography.bib` if a further INSPIRE record is consulted (for example on the resonance couplings in $1/N_c$ counting).
- Downstream: unblocks alp-gii.2, alp-gii.3 and alp-gii.8; the conventions fixed here bind alp-gii.4–alp-gii.7.
- No changes to `latex/`.
