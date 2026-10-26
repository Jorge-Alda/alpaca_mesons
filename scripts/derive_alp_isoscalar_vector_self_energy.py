#!/usr/bin/env python
"""Derive isoscalar vector dressing, pole replacement and on-shell 3π width.

Run with ``.venv/bin/python scripts/derive_alp_isoscalar_vector_self_energy.py``
or append ``--test`` for the minimal unit tests. Inputs and conventions:
memory/alp/lepton-pair/alp-isoscalar-tff-source-decomposition.md and
memory/alp/lepton-pair/alp-three-pion-discontinuity.md and
memory/alp/lepton-pair/alp-three-pion-vector-partial-wave.md. No numerical or
channel-truncation input is selected. Natural units are used; s has GeV²,
the self-energy has GeV², and the width has GeV. The displayed self-energy
notation is supplied by the user through the coordinating agent.
"""

import argparse
import unittest

import sympy as sp


def show(label, expression):
    """Print a step using SymPy's standard LaTeX renderer."""
    print(f"{label}:\n{sp.latex(expression)}\n")


def derive():
    """Construct the algebraic steps, leaving physical input functions free."""
    invariant = sp.Symbol("s", positive=True)
    spectral_invariant = sp.Symbol("x", positive=True)
    bare_mass = sp.Symbol("M_{0i}", positive=True)
    real_self_energy = sp.Symbol(r"\operatorname{Re}\Delta_i^{\mathrm{self}}", real=True)
    imaginary_self_energy = sp.Symbol(r"\operatorname{Im}\Delta_i^{\mathrm{self}}", real=True)
    width = sp.Function(r"\Gamma_i")(invariant)
    decay_integral = sp.Symbol(r"\sum_n\int d\Phi_n\overline{|\mathcal M_{i\to n}|^2}", nonnegative=True)
    width_expression = decay_integral / (2 * sp.sqrt(invariant))
    inverse_propagator = bare_mass**2 - invariant - real_self_energy - sp.I * imaginary_self_energy
    propagator = 1 / inverse_propagator
    propagator_imaginary = sp.simplify(sp.im(sp.expand_complex(propagator)))
    optical_equation = sp.Eq(2 * imaginary_self_energy, decay_integral)
    imaginary_solution = sp.solve(optical_equation, imaginary_self_energy)[0]
    imaginary_width = imaginary_solution.subs(decay_integral, 2 * sp.sqrt(invariant) * width)
    spectral_density = sp.factor(propagator_imaginary.subs(imaginary_self_energy, imaginary_width) / sp.pi)

    subtraction_count = sp.Symbol("N", integer=True, positive=True)
    summation_index = sp.Symbol("k", integer=True, nonnegative=True)
    geometric_sum = sp.Sum(invariant**summation_index / spectral_invariant**(summation_index + 1), (summation_index, 0, subtraction_count - 1)).doit()
    # Cauchy's integrand is derived off its pole, x != s. SymPy's finite
    # geometric sum also returns an x=s branch, which is excluded here.
    geometric_sum = geometric_sum.replace(sp.Piecewise, lambda *branches: branches[-1][0])
    remainder = sp.simplify(1 / (spectral_invariant - invariant) - geometric_sum)
    threshold = sp.Symbol("s_{\\mathrm{th}}", nonnegative=True)
    subtraction_coefficient = sp.IndexedBase("c_i")
    subtraction_polynomial = sp.Sum(subtraction_coefficient[summation_index] * invariant**summation_index, (summation_index, 0, subtraction_count - 1))
    absorptive_part = sp.Function(r"\operatorname{Im}\Delta_i^{\mathrm{self}}")(spectral_invariant)
    dispersion = subtraction_polynomial + sp.Integral(absorptive_part * remainder / sp.pi, (spectral_invariant, threshold, sp.oo))

    regulator = sp.Symbol(r"\varepsilon", positive=True)
    mass = sp.Symbol("M_i", positive=True)
    stable_density = sp.simplify(sp.im(1 / (mass**2 - invariant - sp.I * regulator)) / sp.pi)
    centered_invariant = sp.Symbol("u", real=True)
    density_normalization = sp.integrate(stable_density.subs(invariant, centered_invariant + mass**2), (centered_invariant, -sp.oo, sp.oo))
    residue = sp.Symbol("R_i", real=True)
    pole = residue / (mass**2 - invariant)
    subtracted_pole = residue / mass**2 + invariant * residue / (mass**2 * (mass**2 - invariant))

    # The charged channel has three distinct particles: no 1/3! is present.
    pair_invariant = sp.Symbol("u", positive=True)
    angle_cosine = sp.Symbol("z", real=True)
    plus_mass, minus_mass, neutral_mass = sp.symbols("m_+ m_- m_0", positive=True)
    pair_kallen = sp.expand(pair_invariant**2 + plus_mass**4 + minus_mass**4 - 2 * pair_invariant * plus_mass**2 - 2 * pair_invariant * minus_mass**2 - 2 * plus_mass**2 * minus_mass**2)
    spectator_kallen = sp.expand(invariant**2 + pair_invariant**2 + neutral_mass**4 - 2 * invariant * pair_invariant - 2 * invariant * neutral_mass**2 - 2 * pair_invariant * neutral_mass**2)
    phase_space_density = sp.sqrt(spectator_kallen * pair_kallen) / (256 * sp.pi**3 * invariant * pair_invariant)
    current_norm = (1 - angle_cosine**2) * spectator_kallen * pair_kallen / (16 * pair_invariant)
    spin_average = current_norm / 3
    width_density = sp.factor(phase_space_density * spin_average / (2 * sp.sqrt(invariant)))
    width_prefactor = sp.simplify(width_density / ((1 - angle_cosine**2) * (spectator_kallen * pair_kallen)**sp.Rational(3, 2) / pair_invariant**2))
    pair_energy = (pair_invariant + plus_mass**2 - minus_mass**2) / (2 * sp.sqrt(pair_invariant))
    neutral_energy = (invariant - pair_invariant - neutral_mass**2) / (2 * sp.sqrt(pair_invariant))
    other_pair_invariant = plus_mass**2 + neutral_mass**2 + 2 * pair_energy * neutral_energy - angle_cosine * sp.sqrt(pair_kallen * spectator_kallen) / (2 * pair_invariant)
    on_shell_decay_scalar = sp.Function(r"F_i^{\mathrm{on}}")(pair_invariant, other_pair_invariant.subs(invariant, mass**2))
    on_shell_width = sp.Integral(width_density.subs(invariant, mass**2) * sp.Abs(on_shell_decay_scalar)**2,
                                 (angle_cosine, -1, 1),
                                 (pair_invariant, (plus_mass + minus_mass)**2, (mass - neutral_mass)**2))
    return locals()


def main():
    """Print the conditional derivation and its assumptions."""
    result = derive()
    print("Assumptions: one diagonal transverse response; relativistic states; S=1+iT;")
    print("upper boundary s+i0; D_i=M_0i²-s-Δ_i^self; Schwarz reflection;")
    print("complete open-channel sum, no chosen channel truncation; spin average is in the decay integral.")
    print("The off-shell width is a definition conditional on the selected transverse self-energy.")
    print("No off-shell vertices, mixing, channel phases, subtraction constants or high-energy behavior are fixed.")
    show("1. Channel sum defining the width", sp.Eq(result["width"], result["width_expression"]))
    show("2. Chosen inverse propagator on the upper boundary", result["inverse_propagator"])
    show("3. Optical-theorem absorptive equation in this convention", result["optical_equation"])
    show("4. Solve the absorptive equation", sp.Eq(result["imaginary_self_energy"], result["imaginary_solution"]))
    show("5. Substitute the width definition: the sign is positive", sp.Eq(result["imaginary_self_energy"], result["imaginary_width"]))
    show("6. Imaginary part of the dressed propagator", result["propagator_imaginary"])
    show("7. Dressed propagator spectral density", result["spectral_density"])
    show("8. First N Cauchy terms", result["geometric_sum"])
    show("9. Cauchy remainder after subtraction", result["remainder"])
    show("10. Conditional N-subtracted self-energy representation", result["dispersion"])
    print("The integral in step 10 is evaluated off the cut, then continued as s+i0;")
    print("its subtraction point 0 must lie below the first threshold and N must ensure convergence.")
    print("Each coefficient c_i[k] is free and real under Schwarz reflection.")
    show("11. Regulated stable propagator spectral density", result["stable_density"])
    show("12. Unit weight of the centered Cauchy density on the full real line", result["density_normalization"])
    show("13. Distributional zero-width limit", sp.Eq(sp.Symbol(r"\lim_{\varepsilon\to0^+}\rho_i(s)"), sp.DiracDelta(result["invariant"] - result["mass"]**2)))
    show("14. Stable pole reconstructed by its subtracted spectral delta", sp.Eq(result["pole"], result["subtracted_pole"], evaluate=False))
    show("15. Difference of those two representations", sp.simplify(result["pole"] - result["subtracted_pole"]))
    print("The dressed spectral density replaces the same stable-pole strength; adding both counts it twice.")
    print("For a dressed pole, mass and residue renormalization also depend on the real self-energy and its derivative.")
    print("Dressing this response alone does not derive complete ALP production vertices or multiparticle unitarity.")
    print("Three-pion channel: use the established decay scalar only on shell, s=M_i².")
    print("The scalar argument order below is (u,t); this is a relabeling of the stored (t,u) Dalitz arguments.")
    show("16. Three-pion phase-space density dΦ3/(du dz)", result["phase_space_density"])
    show("17. Spin-averaged squared geometric decay current", result["spin_average"])
    show("18. Assemble flux, phase-space density and spin average", result["width_density"])
    show("19. Extract coefficient multiplying (1-z²)(λ(s,u,m0²)λ(u,m+²,m-²))^(3/2)/u²", result["width_prefactor"])
    show("20. Other pair invariant t(u,z;s)", result["other_pair_invariant"])
    show("21. Physical on-shell charged-three-pion partial width", result["on_shell_width"])
    print("For a running width Γ_i→3π(s), the same integral requires a specified off-shell scalar")
    print("F_i(s;u,t), not merely F_i^on(u,t). Its continuation is not selected here.")
    print("The physical on-shell decay normalization is part of F_i^on; the ALP relative phase is not fixed by this width.")


class SelfEnergyTests(unittest.TestCase):
    """Check sign, subtracted Cauchy identity and stable-pole normalization."""

    def test_optical_sign(self):
        result = derive()
        expected = result["imaginary_self_energy"] / (sp.pi * ((result["bare_mass"]**2 - result["invariant"] - result["real_self_energy"])**2 + result["imaginary_self_energy"]**2))
        self.assertEqual(sp.simplify(result["propagator_imaginary"] / sp.pi - expected), 0)
        self.assertEqual(result["imaginary_width"], sp.sqrt(result["invariant"]) * result["width"])

    def test_subtracted_cauchy(self):
        result = derive()
        expected = result["invariant"]**result["subtraction_count"] / (result["spectral_invariant"]**result["subtraction_count"] * (result["spectral_invariant"] - result["invariant"]))
        self.assertEqual(sp.simplify(result["remainder"] - expected), 0)

    def test_stable_pole(self):
        result = derive()
        self.assertEqual(result["density_normalization"], 1)
        self.assertEqual(sp.simplify(result["pole"] - result["subtracted_pole"]), 0)

    def test_three_pion_prefactor(self):
        result = derive()
        # Independent integer factors: flux 2, spin average 3,
        # phase-space normalization 256, geometric Gram coefficient 16.
        expected = 1 / (sp.Integer(2) * 3 * 256 * 16 * sp.pi**3 * result["invariant"]**sp.Rational(3, 2))
        self.assertEqual(sp.simplify(result["width_prefactor"] - expected), 0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--test", action="store_true")
    arguments = parser.parse_args()
    if arguments.test:
        unittest.main(argv=[__file__])
    else:
        main()
