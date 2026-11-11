# Spec Delta

## Purpose

Defines the all-orders transformation of ALP-dependent quark-level and CCWZ sources under general chiral field redefinitions, so that physical results can be shown to be independent of the chosen coordinates.

## ADDED Requirements

### Requirement: Restricted partonic Lagrangian is closed under redefinitions
The derivation SHALL work within the restricted Lagrangian of the general chiral-redefinition record (QCD, mass function $M(x)$, currents $C_{R,L}(x)$ with one $\partial_\mu a$, $\Phi(x)$, $\Theta(x)$) and SHALL state that omitted operators are not included.

#### Scenario: General $g_{R,L}(x)$ with $g(0)=\mathbf1$
- **WHEN** $q_{R,L}'=g_{R,L}(x)q_{R,L}$ is applied with noncommuting matrices
- **THEN** the transformed functions $M',C'_{R,L},\Theta'$ stay in the same operator set and are given by closed formulas at every Taylor order

### Requirement: Transport matrices give invariant coefficients
The derivation SHALL define $W_{R,L}$ by $\partial_xW_\eta=iC_\eta W_\eta$ and SHALL show that $\widehat M=W_R^\dagger MW_L$ is invariant, with a recursion for its Taylor coefficients.

#### Scenario: Invariance under two successive redefinitions
- **WHEN** two redefinitions $g^{(1)}$ then $g^{(2)}$ are composed
- **THEN** $\widehat M$ and $\widehat U=W_R^\dagger UW_L$ computed from the final sources equal those from the initial sources, verified symbolically through third order

### Requirement: Anomaly and topological terms transform with the redefinition
The derivation SHALL give the shift of $\Theta(x)$ and of the Wess-Zumino and $X_{TC}$ terms induced by a general, not only constant, singlet axial part of $g_{R,L}(x)$, in the normalization of Alda:2026nrd.

#### Scenario: Singlet axial rotation
- **WHEN** $g_R=g_L^\dagger=e^{i\theta(x)}\mathbf1$ with $\theta(x)$ a general function
- **THEN** the shift of $\Theta$ equals the Jacobian of the rotation and the meson-level topological source shifts consistently, so $X_{TC}$-dependent quantities are unchanged

### Requirement: Nonlinear coordinate map is complete
The derivation SHALL give the meson coordinate map induced by $g_{R,L}(x)$ to all orders in the meson fields, not only its tangent-space part, together with the transformed $u_\mu$, $\Gamma_\mu$, $\chi_\pm$.

#### Scenario: Off-shell two-meson contact term
- **WHEN** the off-shell amputated vertex with one ALP and two mesons is computed in original and transformed coordinates
- **THEN** the transformed vertex is reproduced by the pole-vector and nonlinear-map contact terms, and omitting the contact terms fails the off-shell comparison; a noncommuting $U(2)$ calculation may be used as a formal algebraic check

#### Scenario: Physical flavour-changing two-meson amplitude
- **WHEN** a physical on-shell $aPP$ amplitude from an off-diagonal $c_q^{VS}$ is checked
- **THEN** the calculation uses a charge-conserving $U(3)$ channel and compares the amplitude in the original and transformed coordinates; a formal $U(2)$ example is not presented as such a physical channel

### Requirement: Reproduction of known limits
The general results SHALL reduce to the first-order $c_q^{VS}=c_q^V-c_q^S$ and $c_q^{AP}=c_q^A-c_q^P$ combinations of Alda:2026nrd and to the linear pole-alignment shift of the finite-mixing record.

#### Scenario: First-order limit
- **WHEN** the all-orders formulas are expanded to first order in $1/f_a$
- **THEN** they agree with the invariants of Alda:2026nrd
