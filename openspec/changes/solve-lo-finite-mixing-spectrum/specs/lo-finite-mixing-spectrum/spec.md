# Spec Delta

## Purpose

Defines the exact LO spectrum, state compositions and residues of the neutral pseudoscalar nonet plus ALP from the full combined inverse propagator, and the diagnostics that bound its validity.

## ADDED Requirements

### Requirement: Full LO inverse propagator
The solver SHALL build $D(s)=sK-M^2$ in the basis $\Phi=(P_i,a)$ from the LO blocks of `finite-alp-meson-mixing.md`, with kinetic and mass mixing present simultaneously, and SHALL NOT expand in $f_\pi/f_a$ or diagonalize the two off-diagonal blocks separately. The blocks SHALL be taken as inputs of a general-$(K,M^2)$ interface so that NLO or loop corrections can be supplied later.

#### Scenario: LO blocks
- **WHEN** $N_f=3$, $c_q^{A,P}$, $c_G$, $f_a$, $m_0^2$, $m_{a,0}^2$ and $\chi_0$ are given
- **THEN** the solver returns $K$ and $M^2$ equal to the closed forms of `finite-alp-meson-mixing.md`, with $K_{aa}=1+K_{ia}K_{ia}$ and $M^2_{aa}$ including the $c_G^2$ term

### Requirement: Physical poles, compositions and residues
The solver SHALL return all zeros of $\det D(s)$, the eigenvectors normalized by $\Phi^{\mathsf T}K\Phi=1$, the composition of each pole in $(P_i,a)$ and the residue matrix of $D^{-1}(s)$ at each pole.

#### Scenario: Residues and completeness
- **WHEN** all poles are simple
- **THEN** each residue equals $v^{(n)}v^{(n)\mathsf T}$, $\sum_nv^{(n)}v^{(n)\mathsf T}=K^{-1}$, and $D^{-1}(s)=\sum_nv^{(n)}v^{(n)\mathsf T}/(s-m_n^2)$ agrees with a direct matrix inverse at test values of $s$

#### Scenario: Agreement with the single-pole record
- **WHEN** the ALP-like pole is selected
- **THEN** its eigenvector in the normalization $(w_i^{(a)}(m_a^2),1)$ and its LSZ factor match `finite-alp-meson-mixing.md` to numerical or exact-rational precision

### Requirement: Weak-mixing limit
The solver SHALL reproduce the first-order mixing of Alda:2026nrd as $f_\pi/f_a\to0$ and SHALL report the size of the corrections beyond first order.

#### Scenario: Small $K_{ia}$, small $M^2_{ia}$
- **WHEN** $f_a/f_\pi$ is large
- **THEN** the exact eigenvectors and mass shifts agree with the $O(f_\pi/f_a)$ and $O(f_\pi^2/f_a^2)$ expansions, with the residual scaling as the next order

### Requirement: Invariance under field-coordinate changes
The spectrum and residues SHALL be invariant under constant invertible $T$ ($D'=T^{\mathsf T}DT$) and the poles under any $T$ that preserves positivity.

#### Scenario: Random $T$
- **WHEN** a random rational invertible $T$ is applied to $(K,M^2)$
- **THEN** $m_n^2$ are unchanged and eigenvectors and residues map by $\Phi=T\Phi'$ and $R'=T^{-1}RT^{-\mathsf T}$

### Requirement: Stability and domain diagnostics
The solver SHALL report: whether $K$ is positive definite (smallest eigenvalue, $\det K$); whether $M^2$ is positive semidefinite and every $m_n^2$ is real and nonnegative; near-degenerate pairs with the minimal gap and the mixing angle there; and the failure region of the exact formulas (singular $m^2K-M^2_{PP}$ block).

#### Scenario: Unstable vacuum
- **WHEN** parameters give a negative eigenvalue of $M^2$ or a nonpositive eigenvalue of $K$
- **THEN** the solver raises a diagnostic naming the violated condition and does not return the spectrum as physical

#### Scenario: Degeneracy
- **WHEN** two eigenvalues differ by less than a stated tolerance
- **THEN** the solver reports the pair as near-degenerate, gives the splitting $2|M^2_{P\hat a,\rm can}|$ form for the two-level case, and returns a basis of the degenerate subspace instead of a single vector

### Requirement: Pole tracking
The solver SHALL label poles by continuity along a parameter path and report avoided crossings.

#### Scenario: Scan in $m_{a,0}^2$
- **WHEN** $m_{a,0}^2$ crosses a meson mass
- **THEN** the ALP-like branch is identified by its largest $a$-component, the minimal gap is reported, and no crossing of labels occurs

### Requirement: Limits of the LO truncation
The records SHALL state conditions under which the LO chiral expansion or perturbative counting fails even though the quadratic diagonalization is exact, and SHALL state that NLO local terms, loop self-energies, resonance mixing and complex poles are not included.

#### Scenario: Large mixing
- **WHEN** $|K_{ia}|\gtrsim1$ or an eigenvector component of $a$ in a meson-like pole is $O(1)$, or $m_n^2$ approaches the chiral-expansion scale
- **THEN** the report flags the result as LO-only and outside the controlled domain
