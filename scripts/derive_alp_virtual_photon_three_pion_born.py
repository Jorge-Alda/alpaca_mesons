"""Derive the LO charged-pion Born current for a gamma* -> pi+ pi- pi0.

Run with ``.venv/bin/python scripts/derive_alp_virtual_photon_three_pion_born.py``.
The input is the no-mixing LO a Phi^3 interaction of Alda:2026nrd,
recorded in memory/alp/three-pseudoscalar/alp-phi3-lo-chipt-derivation.md.
Natural units, metric (+---), S=1+iT, incoming ALP p and photon q,
and outgoing pion momenta p_+, p_-, p_0 are used. The electromagnetic
coupling is included in the photon amplitude; divide by e for the current.
The neutral generator is left diagonal and symbolic, preserving physical
neutral mixing. This is a gauged LO contribution, not a complete resonance
amplitude or a dispersive solution. Source matrices are diagonal; no new
source relation or isospin approximation is imposed in the general result.
The printed source vertex uses flavour traces of c_q^A, {c_q^P,m_q},
and the pion generators. Explicit diagonal entries are used internally
only to retain the existing independent coefficient calculation.
"""

import sympy as sp


def derive() -> dict[str, sp.Expr]:
    """Expand flavour traces and derive pole and contact Ward contractions."""
    axial_up, axial_down, axial_strange = sp.symbols("c_Au c_Ad c_As")
    spurion_up, spurion_down, spurion_strange = sp.symbols("P_u P_d P_s")
    neutral_up, neutral_down, neutral_strange = sp.symbols("lambda_0u lambda_0d lambda_0s")
    axial = sp.diag(axial_up, axial_down, axial_strange)
    pseudoscalar_mass = sp.diag(spurion_up, spurion_down, spurion_strange)
    positive = sp.zeros(3)
    positive[0, 1] = sp.sqrt(2)
    negative = positive.T
    neutral = sp.diag(neutral_up, neutral_down, neutral_strange)
    generators = (positive, negative, neutral)
    fields = sp.symbols("pi_plus pi_minus pi_zero")
    derivatives = sp.symbols("dpi_plus dpi_minus dpi_zero")
    phi = sum((field * generator for field, generator in zip(fields, generators)), sp.zeros(3))
    dphi = sum((field * generator for field, generator in zip(derivatives, generators)), sp.zeros(3))
    first_commutator = phi * dphi - dphi * phi
    kinetic = sp.expand(sp.trace(axial * (phi * first_commutator - first_commutator * phi)))
    differentiated_weights = []
    for differentiated in range(3):
        monomial = derivatives[differentiated]
        for other in range(3):
            if other != differentiated:
                monomial *= fields[other]
        differentiated_weights.append(sp.expand(kinetic).coeff(monomial))
    mass_weight = sp.expand(sp.trace(phi**3 * pseudoscalar_mass))
    for field in fields:
        mass_weight = mass_weight.coeff(field, 1)
    mass_weight = sp.factor(mass_weight)
    decay_scale, pion_decay_constant, b_zero = sp.symbols("f_a f_pi B_0", nonzero=True)
    alp_mass_squared, alp_dot_neutral = sp.symbols("m_a2 p_dot_p0")
    alp_dot_positive, alp_dot_negative, alp_dot_photon = sp.symbols(
        "p_dot_pplus p_dot_pminus p_dot_q"
    )
    # At either charged-pion attachment, the momentum entering the source
    # vertex is p_+-q or p_--q; their sum with p_0 equals p.
    kinetic_factor = 1 / (6 * decay_scale * pion_decay_constant)
    mass_factor = -b_zero / (6 * decay_scale * pion_decay_constant)
    source_vertex = kinetic_factor * (
        differentiated_weights[0] * alp_dot_positive
        + differentiated_weights[1] * alp_dot_negative
        + differentiated_weights[2] * alp_dot_neutral
    ) + mass_factor * mass_weight
    # Keep the two off-shell source vertices distinct until their extracted
    # flavour coefficients and momentum conservation establish their equality.
    positive_source_vertex = source_vertex.subs(
        alp_dot_positive, alp_dot_positive - alp_dot_photon
    )
    negative_source_vertex = source_vertex.subs(
        alp_dot_negative, alp_dot_negative - alp_dot_photon
    )
    # q is INCOMING: p+q=p_++p_-+p_0. Taking the scalar product
    # with p gives p.p_- = m_a^2 + p.q - p.p_+ - p.p_0.
    # The opposite sign of p.q belongs to an outgoing-photon convention
    # and would also require shifting the pion vertices by +q, not -q.
    momentum_conservation = {
        alp_dot_negative: alp_mass_squared + alp_dot_photon
        - alp_dot_positive - alp_dot_neutral
    }
    positive_source_vertex = sp.factor(positive_source_vertex.subs(momentum_conservation))
    negative_source_vertex = sp.factor(negative_source_vertex.subs(momentum_conservation))
    source_amplitude = positive_source_vertex
    # Gauging d pi+ and d pi- gives opposite electric charges. Extract
    # the source coefficient; its photon-current contribution is -e K times
    # this coefficient times p^mu in the pole convention used below.
    contact_weight = sp.factor(differentiated_weights[0] - differentiated_weights[1])
    photon_squared, photon_dot_positive, photon_dot_negative = sp.symbols("q2 q_dot_pplus q_dot_pminus")
    positive_denominator = photon_squared - 2 * photon_dot_positive
    negative_denominator = photon_squared - 2 * photon_dot_negative
    electric_charge = sp.Symbol("e")
    # q.(2p_r-q)=-( (p_r-q)^2-m_pi^2 ) on the external mass shell.
    contact_current_coefficient = -electric_charge * kinetic_factor * contact_weight
    ward = sp.factor(electric_charge * (
        positive_source_vertex * (2 * photon_dot_positive - photon_squared) / positive_denominator
        - negative_source_vertex * (2 * photon_dot_negative - photon_squared) / negative_denominator
    ) + contact_current_coefficient * alp_dot_photon)
    isospin_substitutions = {neutral_up: 1, neutral_down: -1, neutral_strange: 0}
    charge_up, charge_down, charge_strange = sp.symbols("Q_u Q_d Q_s")
    charge = sp.diag(charge_up, charge_down, charge_strange)
    return {
        "positive_derivative_weight": sp.factor(differentiated_weights[0]),
        "negative_derivative_weight": sp.factor(differentiated_weights[1]),
        "neutral_derivative_weight": sp.factor(differentiated_weights[2]),
        "pseudoscalar_mass_weight": mass_weight,
        "contact_weight": contact_weight,
        "charged_derivative_weight_difference": sp.factor(differentiated_weights[0] - differentiated_weights[1]),
        "positive_source_vertex": positive_source_vertex,
        "negative_source_vertex": negative_source_vertex,
        "source_vertex_difference": sp.factor(positive_source_vertex - negative_source_vertex),
        "contact_current_coefficient": contact_current_coefficient,
        "source_amplitude": source_amplitude,
        "ward_contraction": ward,
        "isospin_source_amplitude": sp.factor(source_amplitude.subs(isospin_substitutions)),
        "isoscalar_source_amplitude": sp.factor(source_amplitude.subs(isospin_substitutions).subs({axial_down: axial_up, spurion_down: spurion_up})),
        "positive_charge_commutator": charge * positive - positive * charge,
        "negative_charge_commutator": charge * negative - negative * charge,
        "isoscalar_charge_commutator": (charge * positive - positive * charge).subs(charge_down, charge_up),
        "positive_denominator": positive_denominator,
        "negative_denominator": negative_denominator,
    }


def main() -> None:
    """Print the source amplitude as traces and the invariant zero checks."""
    print("LO no-mixing diagonal sources; natural units; p+q=p_++p_-+p_0.")
    print("Photon factor e retained; no resonance or rescattering terms.")
    print(r"Flavour trace: \langle M\rangle=\operatorname{Tr}_{f}M.")
    print(r"\Phi=\pi^+\lambda_++\pi^-\lambda_-+\pi^0\lambda_{\pi^0}; "
          r"\lambda_\pm=(\lambda_1\pm i\lambda_2)/\sqrt{2},\quad"
          r"\operatorname{Tr}_{f}(\lambda_i\lambda_j)=2\delta_{ij}.")
    print(r"\lambda_{\pi^0} is diagonal, with the retained physical neutral mixing.")
    print(r"\mathcal V_+=\mathcal V_-=\mathcal V_{a\pi^+\pi^-\pi^0}="
          r"\frac{3p\cdot p_0-m_a^2}{12f_af_\pi}"
          r"\left\langle c_q^A\left("
          r"[\lambda_+,[\lambda_-,\lambda_{\pi^0}]]"
          r"+[\lambda_-,[\lambda_+,\lambda_{\pi^0}]]"
          r"\right)\right\rangle"
          r"-\frac{B_0}{6f_af_\pi}\left\langle\{c_q^P,m_q\}"
          r"\left(\lambda_+\lambda_-\lambda_{\pi^0}"
          r"+\lambda_+\lambda_{\pi^0}\lambda_-"
          r"+\lambda_-\lambda_+\lambda_{\pi^0}"
          r"+\lambda_-\lambda_{\pi^0}\lambda_+"
          r"+\lambda_{\pi^0}\lambda_+\lambda_-"
          r"+\lambda_{\pi^0}\lambda_-\lambda_+\right)\right\rangle.")
    expressions = derive()
    for label in (
        "charged_derivative_weight_difference",
        "source_vertex_difference",
        "contact_current_coefficient",
        "ward_contraction",
        "isoscalar_source_amplitude",
        "isoscalar_charge_commutator",
    ):
        expression = expressions[label]
        print(f"{label}:\n{sp.latex(expression)}\n")
    print(r"J^\mu_{a\to3\pi,\rm Born}=e\mathcal V_{a\pi^+\pi^-\pi^0}"
          r"\left[\frac{(2p_+-q)^\mu}{D_++i0}"
          r"-\frac{(2p_--q)^\mu}{D_-+i0}\right],"
          r"\qquad D_\pm=q^2-2q\cdot p_\pm.")
    print(r"The isoscalar source check uses \lambda_{\pi^0}=\lambda_3 and isoscalar light-quark sources.")
    print(r"The isoscalar-current check uses [Q^{I=0},\lambda_\pm]=0.")
    print("The pole prescription is D_+ + i0 and D_- + i0; Ward identities are algebraic away from poles.")


if __name__ == "__main__":
    main()
