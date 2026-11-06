#!/usr/bin/env python3
"""Derive conditional rho-to-two-pion widths and a subtracted response.

Run: .venv/bin/python scripts/derive_alp_rho_pipi_self_energy.py
Verify: .venv/bin/python scripts/derive_alp_rho_pipi_self_energy.py --test

Inputs: memory/rcht/foundations/rcht-vpp-vpa-photon-pa-expansion.md,
memory/rcht/foundations/vector-gv-fit.md, and
memory/alp/lepton-pair/alp-isoscalar-vector-self-energy.md.
Use natural units, metric (+---), relativistic states, S=1+iT, and the
upper boundary s+i0. Masses and widths are in GeV; s and self-energies are
in GeV². Only the exact-isospin charged two-pion cut and the leading G_V
vertex are included. Off-shell tensor normalization is a stated convention,
not an assertion of field-redefinition-independent off-shell physics.
"""

import sys
import unittest

import sympy as sp


def show(label, expression):
    """Print an intermediate expression with SymPy's LaTeX printer."""
    print(f"{label}:\n{sp.latex(expression)}\n")


def derive():
    """Evaluate the rho flavor trace, two-body cut and response identities."""
    invariant = sp.Symbol("s", positive=True)
    spectral_invariant = sp.Symbol("x", positive=True)
    pion_mass = sp.Symbol(r"M_\pi", positive=True)
    rho_mass = sp.Symbol(r"M_\rho", positive=True)
    tensor_mass = sp.Symbol(r"M_{\mathrm{tensor}}", positive=True)
    bare_rho_mass = sp.Symbol(r"M_{0\rho}", positive=True)
    decay_constant = sp.Symbol("F", positive=True)
    vector_coupling = sp.Symbol("G_V", real=True, nonzero=True)
    pion_momentum = sp.Symbol(r"q_{\pi\pi}", positive=True)
    pole_pion_momentum = sp.Symbol(r"q_{\pi\pi}(M_\rho^2)", positive=True)
    angle_cosine = sp.Symbol("z", real=True)
    azimuth = sp.Symbol(r"\varphi", real=True)

    # U=exp(i sqrt(2) Phi/F), Tr(T_a T_b)=delta_ab.  The two
    # charged pion fields are distinct, so the vertex is counted once.
    rho_generator = sp.diag(1, -1, 0) / sp.sqrt(2)
    pion_plus_generator = sp.zeros(3)
    pion_plus_generator[0, 1] = 1
    pion_minus_generator = pion_plus_generator.T
    flavor_trace = sp.trace(rho_generator * (pion_plus_generator * pion_minus_generator - pion_minus_generator * pion_plus_generator))
    vertex_coefficient = sp.sqrt(2) * vector_coupling * flavor_trace / decay_constant**2

    # The common overall i and the sign fixed by charged-field phases do not
    # affect the width.  The tensor polarization is i(q_mu eps_nu-q_nu
    # eps_mu)/sqrt(s) for a unit spin-one projector at timelike q²=s.
    metric = sp.diag(1, -1, -1, -1)
    rho_momentum = sp.Matrix([sp.sqrt(invariant), 0, 0, 0])
    pion_plus_momentum = sp.Matrix([sp.sqrt(invariant) / 2, 0, 0, pion_momentum])
    pion_minus_momentum = sp.Matrix([sp.sqrt(invariant) / 2, 0, 0, -pion_momentum])
    polarization = sp.Matrix([0, 0, 0, 1])
    scalar_product = lambda left, right: (left.T * metric * right)[0]
    tensor_contraction = sp.expand(
        scalar_product(rho_momentum, pion_plus_momentum) * scalar_product(polarization, pion_minus_momentum)
        - scalar_product(rho_momentum, pion_minus_momentum) * scalar_product(polarization, pion_plus_momentum)
    )
    pion_difference = pion_plus_momentum - pion_minus_momentum
    unit_tensor_amplitude = sp.simplify(vertex_coefficient * tensor_contraction / sp.sqrt(invariant))
    unit_tensor_current_coefficient = sp.simplify(unit_tensor_amplitude / scalar_product(polarization, pion_difference))
    fixed_tensor_amplitude = sp.simplify(vertex_coefficient * tensor_contraction / tensor_mass)
    fixed_tensor_current_coefficient = sp.simplify(fixed_tensor_amplitude / scalar_product(polarization, pion_difference))
    transverse_spin_sum = sp.simplify(
        -scalar_product(pion_difference, pion_difference)
        + scalar_product(rho_momentum, pion_difference) ** 2 / invariant
    )
    averaged_amplitude_squared = sp.factor(unit_tensor_current_coefficient**2 * transverse_spin_sum / 3)
    fixed_averaged_amplitude_squared = sp.factor(fixed_tensor_current_coefficient**2 * transverse_spin_sum / 3)
    two_body_angular_density = pion_momentum / (16 * sp.pi**2 * sp.sqrt(invariant))
    integrated_two_body_phase_space = sp.integrate(two_body_angular_density, (azimuth, 0, 2 * sp.pi), (angle_cosine, -1, 1))
    unit_tensor_width_in_momentum = sp.factor(averaged_amplitude_squared * integrated_two_body_phase_space / (2 * sp.sqrt(invariant)))
    fixed_tensor_width_in_momentum = sp.factor(fixed_averaged_amplitude_squared * integrated_two_body_phase_space / (2 * sp.sqrt(invariant)))
    pion_momentum_squared = sp.solve(sp.Eq(pion_plus_momentum[0]**2 - pion_momentum**2, pion_mass**2), pion_momentum**2)[0]
    physical_pion_momentum = sp.sqrt(pion_momentum_squared)
    threshold = 4 * pion_mass**2
    unit_tensor_width = sp.factor(unit_tensor_width_in_momentum.subs(pion_momentum, physical_pion_momentum))
    fixed_tensor_width = sp.factor(fixed_tensor_width_in_momentum.subs(pion_momentum, physical_pion_momentum))
    width_with_threshold = sp.Piecewise((0, invariant <= threshold), (unit_tensor_width, True))

    # The fixed-M_tensor continuation reproduces the project's on-shell
    # tensor-to-Proca formula; it is a separate off-shell choice.
    fixed_width_at_pole = fixed_tensor_width_in_momentum.subs({invariant: rho_mass**2, pion_momentum: pole_pion_momentum})
    recorded_on_shell_width = flavor_trace**2 * vector_coupling**2 * rho_mass**2 * pole_pion_momentum**3 / (12 * sp.pi * decay_constant**4 * tensor_mass**2)
    unit_width_at_pole = unit_tensor_width_in_momentum.subs({invariant: rho_mass**2, pion_momentum: pole_pion_momentum})
    unit_normalized_shape = sp.factor(unit_tensor_width_in_momentum / unit_width_at_pole)
    fixed_normalized_shape = sp.factor(fixed_tensor_width_in_momentum / fixed_width_at_pole)
    phenomenological_shape = rho_mass**2 / invariant * (pion_momentum / pole_pion_momentum)**3

    imaginary_self_energy = sp.Symbol(r"\operatorname{Im}\Delta_\rho^{\mathrm{self}}", nonnegative=True)
    real_self_energy = sp.Symbol(r"\operatorname{Re}\Delta_\rho^{\mathrm{self}}", real=True)
    decay_integral = 2 * sp.sqrt(invariant) * unit_tensor_width_in_momentum
    optical_equation = sp.Eq(2 * imaginary_self_energy, decay_integral)
    imaginary_from_unit_width = sp.solve(optical_equation, imaginary_self_energy)[0]
    imaginary_above_threshold = sp.factor(imaginary_from_unit_width.subs(pion_momentum, physical_pion_momentum))
    inverse_response = bare_rho_mass**2 - invariant - real_self_energy - sp.I * imaginary_self_energy
    spectral_density = sp.factor(sp.im(sp.expand_complex(1 / inverse_response)) / sp.pi)

    # Cauchy subtraction at zero, below the two-pion threshold.
    subtraction_count = sp.Symbol("N", integer=True, positive=True)
    summation_index = sp.Symbol("k", integer=True, nonnegative=True)
    cauchy_terms = sp.summation(invariant**summation_index / spectral_invariant**(summation_index + 1), (summation_index, 0, subtraction_count - 1))
    # SymPy nests the geometric sum's x=s branch inside a product; x=s is
    # excluded from the off-cut identity used to construct the dispersion.
    cauchy_terms = cauchy_terms.replace(sp.Piecewise, lambda *branches: branches[-1][0])
    cauchy_remainder = sp.simplify(1 / (spectral_invariant - invariant) - cauchy_terms)
    subtraction_coefficients = sp.IndexedBase(r"c_{\rho}")
    absorptive_function = sp.Function(r"\operatorname{Im}\Delta_\rho^{\mathrm{self}}")(spectral_invariant)
    subtraction_polynomial = sp.Sum(subtraction_coefficients[summation_index] * invariant**summation_index, (summation_index, 0, subtraction_count - 1))
    general_dispersion_integrand = absorptive_function * cauchy_remainder / sp.pi
    principal_value = sp.Function(r"\mathrm{PV}")
    general_real_dispersion = subtraction_polynomial + principal_value(sp.Integral(general_dispersion_integrand, (spectral_invariant, threshold, sp.oo)))
    actual_absorptive_part = imaginary_above_threshold.subs(invariant, spectral_invariant)
    twice_subtracted_integrand = actual_absorptive_part * invariant**2 / (sp.pi * spectral_invariant**2 * (spectral_invariant - invariant))
    thrice_subtracted_integrand = actual_absorptive_part * invariant**3 / (sp.pi * spectral_invariant**3 * (spectral_invariant - invariant))
    twice_subtracted_uv = sp.simplify(sp.limit(spectral_invariant * twice_subtracted_integrand, spectral_invariant, sp.oo))
    thrice_subtracted_uv = sp.simplify(sp.limit(spectral_invariant**2 * thrice_subtracted_integrand, spectral_invariant, sp.oo))
    thrice_subtracted_real = sum(subtraction_coefficients[index] * invariant**index for index in range(3)) + principal_value(sp.Integral(thrice_subtracted_integrand, (spectral_invariant, threshold, sp.oo)))

    # Real-axis conditions are a convention for an unstable resonance.  The
    # derivative is taken after forming the real boundary-value function.
    real_self_energy_function = sp.Function(r"\operatorname{Re}\Delta_\rho^{\mathrm{self}}")
    real_at_pole = real_self_energy_function(rho_mass**2)
    real_derivative_at_pole = sp.Subs(sp.Derivative(real_self_energy_function(invariant), invariant), invariant, rho_mass**2)
    subtracted_real_self_energy = real_self_energy_function(invariant) - real_at_pole - (invariant - rho_mass**2) * real_derivative_at_pole
    mass_condition = sp.Eq(bare_rho_mass**2 - rho_mass**2 - real_at_pole, 0)
    bare_mass_solution = sp.solve(mass_condition, bare_rho_mass**2)[0]
    real_inverse_slope = -1 - real_derivative_at_pole
    second_sheet_pole = sp.Symbol("s_p")
    second_sheet_self_energy = sp.Function(r"\Delta_\rho^{\mathrm{self,II}}")
    second_sheet_inverse = bare_rho_mass**2 - invariant - second_sheet_self_energy(invariant)
    complex_pole_residue = sp.simplify(1 / sp.diff(second_sheet_inverse, invariant).subs(invariant, second_sheet_pole))

    regulator = sp.Symbol(r"\varepsilon", positive=True)
    stable_density = regulator / (sp.pi * ((rho_mass**2 - invariant)**2 + regulator**2))
    stable_weight = sp.integrate(stable_density, (invariant, 0, sp.oo))
    stable_limit = sp.limit(stable_weight, regulator, 0, dir="+")
    return locals()


def main():
    """Print the derivation, including its off-shell and subtraction limits."""
    result = derive()
    print("Exact isospin; charged pions are distinct; only the LO G_V rho-pipi vertex is retained.")
    print("Unit spin-one tensor projector: i(q_mu eps_nu-q_nu eps_mu)/sqrt(s), for timelike s>0.")
    print("Natural units: F, G_V, masses and widths in GeV; s and self-energies in GeV².")
    show("1. Rho charged-pion flavor trace", result["flavor_trace"])
    show("2. Tensor-vertex coefficient apart from its overall phase", result["vertex_coefficient"])
    show("3. Rest-frame tensor contraction", result["tensor_contraction"])
    show("4. Unit-projector coefficient of epsilon dot (p_plus-p_minus)", result["unit_tensor_current_coefficient"])
    show("5. Transverse polarization sum", result["transverse_spin_sum"])
    show("6. Spin-averaged squared amplitude", result["averaged_amplitude_squared"])
    show("7. Integrated two-body phase space", result["integrated_two_body_phase_space"])
    show("8. Width from phase space and flux", result["unit_tensor_width_in_momentum"])
    show("9. On-shell pion momentum squared", result["pion_momentum_squared"])
    show("10. Unit-projector width with physical threshold", result["width_with_threshold"])
    show("11. Threshold coefficient above 4 M_pi²", sp.simplify(result["unit_tensor_width"] / (result["invariant"] - result["threshold"])**sp.Rational(3, 2)))
    print("The tensor normalization above is a virtual spin-one projector convention, not an external unstable-particle LSZ state.")
    show("12. Fixed-tensor-mass continuation of the current coefficient", result["fixed_tensor_current_coefficient"])
    show("13. Fixed-tensor-mass width above threshold", result["fixed_tensor_width"])
    show("14. Difference from the recorded on-shell tensor-to-Proca width", sp.simplify(result["fixed_width_at_pole"] - result["recorded_on_shell_width"]))
    show("15. Unit-projector width divided by its value at M_rho²", result["unit_normalized_shape"])
    show("16. Fixed-tensor-mass width divided by its value at M_rho²", result["fixed_normalized_shape"])
    show("17. Previously recorded phenomenological P-wave shape", result["phenomenological_shape"])
    print("The three shapes agree at s=M_rho² but differ away from it; no off-shell equivalence is assumed.")
    show("18. Optical theorem in the chosen self-energy convention", result["optical_equation"])
    show("19. Im Delta_rho from the two-pion cut, above threshold", result["imaginary_above_threshold"])
    show("20. Positive response spectral density in terms of Im Delta_rho", result["spectral_density"])
    show("21. Cauchy subtraction remainder for x != s", result["cauchy_remainder"])
    show("22. General N-subtracted real part; PV only on the cut", result["general_real_dispersion"])
    show("23. UV coefficient of x times the twice-subtracted integrand", result["twice_subtracted_uv"])
    show("24. UV coefficient of x² times the thrice-subtracted integrand", result["thrice_subtracted_uv"])
    show("25. Three-subtraction real part for the unregulated unit-projector cut", result["thrice_subtracted_real"])
    print("For this unregulated cut N=2 is logarithmically divergent and N=3 converges; a form factor or UV completion can change this requirement.")
    show("26. Real-axis mass condition solved for M_0rho²", sp.Eq(result["bare_rho_mass"]**2, result["bare_mass_solution"]))
    show("27. Real-axis inverse-response slope at M_rho²", result["real_inverse_slope"])
    show("28. Real self-energy with value and slope subtracted at M_rho²", result["subtracted_real_self_energy"])
    show("29. Subtracted real value at M_rho²", sp.simplify(result["subtracted_real_self_energy"].subs(result["invariant"], result["rho_mass"]**2)))
    show("30. Subtracted real slope at M_rho²", sp.simplify(sp.diff(result["subtracted_real_self_energy"], result["invariant"]).subs(result["invariant"], result["rho_mass"]**2)))
    show("31. Residue of the second-sheet response in the (s-s_p)^(-1) convention", result["complex_pole_residue"])
    print("A unit real-axis slope is a field convention and must be reflected in the numerator; it does not locate the complex pole.")
    print("The third subtraction remains an independent input for the unregulated cut.")
    print("On the cut, differentiate the real boundary-value function, not a naive divergent squared-pole integral.")
    show("32. Regulated stable spectral density", result["stable_density"])
    show("33. Its weight on s>=0", result["stable_weight"])
    show("34. Zero-regulator weight for M_rho²>0", result["stable_limit"])
    print("The delta-function limit assumes residual local self-energy terms vanish or are absorbed consistently as G_V tends to zero.")
    print("These response identities do not supply a complete pion form factor, ALP production wave, or mixed-term matching.")


class RhoSelfEnergyTests(unittest.TestCase):
    """Check independent phase-space, threshold, subtraction, and pole identities."""

    @classmethod
    def setUpClass(cls):
        cls.result = derive()

    def test_flavor_spin_and_on_shell_normalization(self):
        result = self.result
        self.assertEqual(result["flavor_trace"], sp.sqrt(2))
        self.assertEqual(result["transverse_spin_sum"], 4 * result["pion_momentum"]**2)
        self.assertEqual(sp.simplify(result["fixed_width_at_pole"] - result["recorded_on_shell_width"]), 0)

    def test_threshold_and_positive_cut(self):
        result = self.result
        expected = result["vector_coupling"]**2 / (48 * sp.pi * result["decay_constant"]**4)
        self.assertEqual(sp.simplify(result["unit_tensor_width"] / (result["invariant"] - result["threshold"])**sp.Rational(3, 2) - expected), 0)
        self.assertEqual(result["width_with_threshold"].subs(result["invariant"], 2 * result["pion_mass"]**2), 0)
        self.assertTrue(sp.ask(sp.Q.positive(result["unit_tensor_width_in_momentum"])))
        self.assertTrue(sp.ask(sp.Q.positive(result["imaginary_from_unit_width"])))

    def test_optical_sign_and_spectral_positivity(self):
        result = self.result
        imaginary = result["imaginary_self_energy"]
        real_denominator = result["bare_rho_mass"]**2 - result["invariant"] - result["real_self_energy"]
        expected = imaginary / (sp.pi * (real_denominator**2 + imaginary**2))
        self.assertEqual(sp.simplify(result["spectral_density"] - expected), 0)
        self.assertEqual(sp.simplify(result["imaginary_from_unit_width"] - sp.sqrt(result["invariant"]) * result["unit_tensor_width_in_momentum"]), 0)

    def test_subtractions_and_real_axis_conditions(self):
        result = self.result
        expected = result["invariant"]**result["subtraction_count"] / (result["spectral_invariant"]**result["subtraction_count"] * (result["spectral_invariant"] - result["invariant"]))
        self.assertEqual(sp.simplify(result["cauchy_remainder"] - expected), 0)
        self.assertNotEqual(result["twice_subtracted_uv"], 0)
        self.assertNotEqual(result["thrice_subtracted_uv"], 0)
        self.assertEqual(sp.simplify(result["subtracted_real_self_energy"].subs(result["invariant"], result["rho_mass"]**2)), 0)
        self.assertEqual(sp.simplify(sp.diff(result["subtracted_real_self_energy"], result["invariant"]).subs(result["invariant"], result["rho_mass"]**2)), 0)

    def test_stable_pole_and_normalized_shapes(self):
        result = self.result
        self.assertEqual(result["stable_limit"], 1)
        self.assertEqual(sp.simplify(result["fixed_normalized_shape"] / result["unit_normalized_shape"] - result["invariant"] / result["rho_mass"]**2), 0)
        self.assertEqual(sp.simplify(result["phenomenological_shape"] / result["unit_normalized_shape"] - result["rho_mass"]**2 / result["invariant"]), 0)


if __name__ == "__main__":
    if sys.argv[1:] == ["--test"]:
        unittest.main(argv=[__file__])
    elif not sys.argv[1:]:
        main()
    else:
        raise SystemExit("Use no arguments for the derivation or --test for its checks.")
