# Proposal

## Why

Beads issue alp-sjp.2 (GitHub issue 21) asks for the physical poles, state compositions and residues of the combined neutral pseudoscalar-nonet plus ALP system, with derivative and mass mixing treated together and no expansion in $f_\pi/f_a$. The records [finite-alp-meson-mixing.md](../../../memory/alp/shared/finite-alp-meson-mixing.md) and [alp-inverse-propagator-congruence.md](../../../memory/alp/shared/alp-inverse-propagator-congruence.md) give the LO blocks $K$, $M^2$ and the covariance of $D(s)=sK-M^2$ under $\Phi=T\Phi'$, and alp-sjp.1 is closed. What is missing is a reproducible solver that returns the full spectrum, eigenvectors and residues, checks the weak-mixing limit, and reports where the LO problem is unstable or degenerate. alp-sjp.3, .5 and .6 depend on it.

## What Changes

- Build the LO inverse propagator for $\Phi_I=(P_i,a)$, $N_f=3$ neutral generators (including the singlet), from the blocks of `finite-alp-meson-mixing.md`: no NLO local terms, no loop self-energies, no resonance mixing.
- Solve the generalized eigenproblem $(M^2-m^2K)\,\Phi=0$ exactly: all physical $m_n^2$, eigenvectors normalized with $K$, and residues of $D^{-1}(s)$.
- Diagnostics: positivity of $K$ (at LO $\det K=1$), nonnegativity of $M^2$ (vacuum stability), near-degeneracy and avoided crossings, and pole tracking along a continuous path in $f_\pi/f_a$ or in $m_{a,0}^2$.
- Checks: weak-mixing limit (recovers the $O(1/f_a)$ results of Alda:2026nrd), invariance under allowed constant field-coordinate changes, agreement with the single-pole eigenvector and normalization of `finite-alp-meson-mixing.md`, and the completeness relation $\sum_nv^{(n)}v^{(n)\mathsf T}=K^{-1}$.
- State when the LO chiral truncation is not reliable even though the quadratic diagonalization is exact.
- Memory record(s) and Python script(s) per the repository protocol.
- **Out of scope (deferred to a future extension):** NLO local terms, loop self-energies, resonance mixing, momentum-dependent (non-polynomial) $D(s)$, and the full-width/complex-pole treatment. The solver interface SHALL accept general $K$, $M^2$ so the extension only supplies new blocks.

## Capabilities

### New Capabilities
- `lo-finite-mixing-spectrum`: exact LO physical poles, compositions and residues of the combined $(P_i,a)$ system, with stability, degeneracy and tracking diagnostics.

### Modified Capabilities
<!-- None: openspec/specs/ is empty. -->

## Impact

- New memory records under `memory/alp/shared/` with `INDEX.md` updated.
- New scripts under `scripts/` (written by the `coder` agent), reusing `derive_alp_lo_finite_mixing.py`, `derive_alp_finite_mixing.py`, `derive_alp_inverse_propagator_congruence.py`.
- The Beads issue alp-sjp.2 is only partly resolved: a note must record that NLO terms, loops and resonance mixing remain open. Do not close it; track the deferred part with a new bead.
- Provides the solver consumed by alp-sjp.3, .5, .6; alp-sjp.4 supplies the numerical inputs (here only illustrative or symbolic ones are used).
