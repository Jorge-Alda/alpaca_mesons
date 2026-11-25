# Spec Delta

## Purpose

Defines the analytical, implicit crossing between the one-loop and the exact tree-level size of an ALP-sourced three-meson amplitude in the Murayama model and its DiLuzio generalization, with the radial mode included. It covers the general result in complex radial-phase invariants and its up-only specialization, with no expansion in $\epsilon=f_\pi/f_a$ in the main result.

## ADDED Requirements

### Requirement: General result in phi-valued invariants
The record SHALL derive the implicit crossing equation first in the invariants that are functions of the complex radial-phase field $\varphi$, before any parametrization by charges or by a specific quark-mass spurion. The DiLuzio $n_q$ parametrization SHALL appear only as a special case of this general result.

#### Scenario: Order of derivation
- **WHEN** the crossing equation is inspected
- **THEN** its general $\varphi$-valued form is stated before any specialization, and the $n_q$ and up-only forms are obtained from it

### Requirement: Loop-size estimate of the ALP and radial loops
The one-loop size in the crossing equation SHALL be a defensible estimate built from the loops that contain at least one internal ALP-like or radial-mode line, with their $\epsilon$ powers, reference scales and invariant contractions stated. Meson-only loops SHALL be excluded from the numerator, because they carry no power of $f_\pi/f_a$ beyond the ordinary chiral correction. The estimate SHALL NOT be truncated to the leading mixing-dependent order in $\epsilon$.

#### Scenario: Meson-only loop
- **WHEN** a loop contains no ALP-like or radial internal line
- **THEN** it is excluded from the numerator and is not evaluated

#### Scenario: Loop with a radial internal line
- **WHEN** a loop contains an internal radial-mode propagator
- **THEN** it is included in the numerator with its vertices, its propagator and its $\epsilon$ power stated, and is not evaluated

### Requirement: No loop computation
The change SHALL NOT compute any loop contribution. It SHALL NOT evaluate a loop integral, a scalar loop function, a loop coefficient, a derivative of a loop with respect to the renormalization scale, or a one-loop amplitude in any frame. Only the class classification, the $\epsilon$ powers, the invariant contractions and the size prescription SHALL be stated.

#### Scenario: Loop value requested
- **WHEN** a derivation would require a loop value to proceed
- **THEN** the record states the estimate in terms of the prescription and does not evaluate the loop value

### Requirement: Naturalness assumption for the neglected term
The record SHALL state that the neglected next-order term, including any renormalization-scale-dependent and counterterm part, is taken to have the size of the loop-size estimate, with no accidental cancellation. It SHALL list the cancellations for which the estimate overestimates the term. It SHALL NOT analyse what counterterms or $L_i$ absorb, and the estimate SHALL contain no renormalization scale $\mu$.

#### Scenario: Cancellation among classes
- **WHEN** two loop classes have the same power of $\epsilon$ and their invariant contractions can cancel
- **THEN** the record states the cancellation condition and that the estimate does not apply there

#### Scenario: Scale in the estimate
- **WHEN** the estimate is written
- **THEN** it contains no $\mu$ and no absorption analysis is given

### Requirement: Exact tree denominator with radial exchange
The tree-level reference SHALL be the exact finite-mixing amplitude including radial-mode exchange, with daughter poles and exchange propagators kept exact, for $\widehat\Theta_1=0$, the case of the Murayama model and its generalizations.

#### Scenario: Radial exchange in the tree
- **WHEN** the tree amplitude is evaluated for a three-meson final state
- **THEN** it includes the three radial-mode exchange channels with their exact propagators

### Requirement: Symbolic threshold and fixed kinematics
The criterion SHALL keep the threshold $R_{\max}$ symbolic and SHALL use the fixed independent Dalitz invariants of the existing ALP-only derivation.

#### Scenario: Threshold
- **WHEN** the crossing equation is stated
- **THEN** $R_{\max}$ appears as a symbol and no numerical value is assumed

### Requirement: Up-only specialization with analytical parameters
The general result SHALL be specialized to the up-only model with $f_a=\sqrt2\,m_u/\kappa_u$, keeping $\kappa_u$ and the radial quartic coupling $\lambda$ as analytical parameters, and SHALL give the radial mass $m_\sigma$ and the ALP-like mass $m_a$ in that model in analytical form.

#### Scenario: Specialization
- **WHEN** the up-only case is derived from the general result
- **THEN** the specialized equation depends on $\kappa_u$ and $\lambda$ symbolically and no value is substituted

### Requirement: Frame invariance
The loop-to-tree ratio of the general result SHALL be shown, with an executable symbolic check, to be unchanged under a constant axial redefinition of the chiral frame.

#### Scenario: Constant axial redefinition
- **WHEN** the sources are transformed by a constant axial redefinition
- **THEN** the loop-to-tree ratio agrees in both frames

### Requirement: Stated caveats not computed
The record SHALL state, without computing, the near-degeneracy of the ALP-like state with $\eta'$ and the validity of leading-order chiral perturbation theory for $m_a$ of about $1$ to $2.3$ GeV.

#### Scenario: Caveat statement
- **WHEN** the validity of the crossing is inspected
- **THEN** both caveats are stated and no value for either region is derived
