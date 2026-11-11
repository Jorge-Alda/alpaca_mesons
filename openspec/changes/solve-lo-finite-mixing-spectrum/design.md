# Design

## Context

The LO quadratic theory is fixed in `finite-alp-meson-mixing.md`: $K_{ij}=\delta_{ij}$, $K_{aa}=1+K_{ia}K_{ia}$, and $M^2$ with the $\chi_0$, $m_0^2$, $c_G$ and $c_q^P$ terms. `alp-inverse-propagator-congruence.md` proves covariance under constant $T$. Existing scripts (`derive_alp_lo_finite_mixing.py`, `derive_alp_finite_mixing.py`) verify the blocks but do not return a full spectrum with residues and diagnostics. The NLO record `memory/nlo/neutral-meson-matrices.md` and the loop records are explicitly outside this change.

## Goals / Non-Goals

**Goals:**
- Exact LO poles, $K$-normalized eigenvectors and residues, from $(K,M^2)$ jointly.
- Stability, degeneracy, tracking diagnostics and a statement of the LO validity domain.
- A solver interface taking arbitrary real symmetric $(K,M^2)$, so NLO local terms can be added later without rewriting.

**Non-Goals:**
- NLO local terms, loop self-energies, resonance mixing, momentum-dependent $D(s)$, complex poles/widths.
- Decay amplitudes (alp-sjp.3), numerical phenomenological inputs (alp-sjp.4), performance optimization (alp-sjp.5), shifts of physical observables (alp-sjp.6).

## Decisions

1. **Generalized symmetric eigenproblem.** Since $K>0$ at LO, solve $M^2v=m^2Kv$ by Cholesky reduction $K=LL^{\mathsf T}$, $L^{-1}M^2L^{-\mathsf T}$ symmetric; eigenvectors are $K$-orthonormal and residues are $v^{(n)}v^{(n)\mathsf T}$. Alternative (root-finding of $\det D$) rejected: loses degeneracies and is not needed while $D$ is linear in $s$. The linear-in-$s$ form is what the NLO extension will break; the interface keeps a hook for a general $D(s)$ solver.
2. **Two implementations.** Exact rational/symbolic with SymPy for $N_f=3$ test data and the weak-mixing series; floating-point NumPy/mpmath for scans. Cross-check one against the other. Scripts are written by the `coder` agent; units: masses squared in GeV$^2$, $f_\pi,f_a$ in GeV, mixing parameters dimensionless.
3. **Kinetic metric analytically.** $K=\begin{pmatrix}1&K_{\cdot a}\\K_{\cdot a}^{\mathsf T}&1+K_{ia}K_{ia}\end{pmatrix}$ has $\det K=1$, so positivity of $K$ is automatic at LO for real $K_{ia}$; the diagnostic is still implemented for the general interface. Consequence recorded in the memory note (stability reduces to $M^2\ge0$ at LO).
4. **Weak-mixing check** against the first-order expansion of Alda:2026nrd (via `alp-dressed-mass-first-order-invariants.md` and the LO extracts), through symbolic series in $\epsilon=f_\pi/f_a$ with $c_q^{A,P},c_G$ fixed.
5. **Degeneracy handling.** Report gap and, for exact degeneracy, an arbitrary $K$-orthonormal basis plus the statement that the composition of individual states is not unique while the projector is.
6. **Pole tracking** via maximal overlap of $K$-weighted eigenvectors at consecutive parameter points with adaptive step refinement near small gaps.
7. **Validity domain** stated in terms of the eigenvector components, $|K_{ia}|$, and the chiral scale; heuristic and labelled as such (assumption), not a derived bound.

8. **Neutral basis.** The solver works in $(P_3,P_8,P_0,P_6,P_7,a)$ (6x6), not in all nine $U(3)$ generators. This is exact for $c_q^{A,P}$ and $\chi_0$ commuting with $Q$ (only d-s off-diagonal entries), shown symbolically by the coder scripts; the charged fields decouple. With diagonal $\chi_0$ (always reachable by a d-s rotation) the flavour-diagonal and $K^0$ blocks mix only through $a$. The block shrinks to 5x5 when $\operatorname{Im}(c^A_{ds}c^{P*}_{ds})=0$, and to 4x4 for flavour-diagonal sources. For any phases, $s=m_{K^0}^2$ is an exact pole with zero ALP component. Charge-violating sources are outside scope.

## Risks / Trade-offs

- LO is not a controlled description when $f_a\simeq f_\pi$ and $O(1)$ mixing arises: NLO terms of relative size $m^2/\Lambda_\chi^2$ are not suppressed by $f_\pi/f_a$. The result is a diagnostic of the LO model, not a prediction.
- Numerical inputs are illustrative until alp-sjp.4; no physical claim about a specific parameter point is made.
- The solver interface must remain general; risk of over-fitting to the LO closed forms is mitigated by tests with random $(K,M^2)$.
