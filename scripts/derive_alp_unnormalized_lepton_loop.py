"""Derive the unnormalized ALP dilepton loop and its pointlike UV behavior.

Run .venv/bin/python scripts/derive_alp_unnormalized_lepton_loop.py
or append --test. Inputs: indexed Messerli:2025rnv functional and existing
ALP normalization. Natural units (GeV), (+---), epsilon^0123=+1, leading QED.
No hadronic TFF, subtraction constant, or local Wilson coefficient is chosen.
"""

from functools import lru_cache
import sys
import unittest

import sympy as sp

from derive_alp_lepton_pair_normalization import derive_normalization


@lru_cache(maxsize=1)
def derive() -> dict[str, sp.Expr]:
    """Remove normalization, derive the cut/width, and integrate the UV log."""
    invariant_mass_squared, lepton_mass, alpha = sp.symbols(
        "s m_l alpha_em", positive=True)
    photon_one_squared, photon_two_squared, lepton_denominator = sp.symbols(
        "k_squared (q-k)_squared D_l", real=True)
    scalar_product = sp.Symbol("q_dot_k", real=True)
    real_photon_form_factor = sp.Symbol("F_agg")
    full_form_factor = sp.Function("F_a_gamma_gamma")(photon_one_squared,
                                                       photon_two_squared)
    numerator = invariant_mass_squared * photon_one_squared - scalar_product**2
    denominator = photon_one_squared * photon_two_squared * lepton_denominator
    reduced_integrand = 2 * sp.I / (sp.pi**2 * invariant_mass_squared) * (
        numerator / denominator) * full_form_factor / real_photon_form_factor
    # The signed loop-coefficient convention is declared here; matching a
    # local Lagrangian coefficient to it is a separate calculation.
    unnormalized_integrand = sp.cancel(
        2 * alpha**2 * lepton_mass * real_photon_form_factor * reduced_integrand)

    source_one, source_two = sp.symbols("c_1 c_2")
    first_response = sp.Function("F_1")(photon_one_squared, photon_two_squared)
    second_response = sp.Function("F_2")(photon_one_squared, photon_two_squared)
    source_resolved_integrand = sp.expand(unnormalized_integrand.subs(
        full_form_factor, source_one * first_response + source_two * second_response))
    source_linearity_difference = sp.simplify(source_resolved_integrand
        - source_one * unnormalized_integrand.subs(full_form_factor, first_response)
        - source_two * unnormalized_integrand.subs(full_form_factor, second_response))
    # A formal counterexample, not a selected hadronic model: equal values
    # at zero do not imply equal functions away from the real-photon point.
    remainder = first_response - second_response
    cancellation_integrand = sp.simplify(unnormalized_integrand.subs(
        full_form_factor, remainder))
    cancellation_at_origin = remainder.subs({photon_one_squared: 0,
        photon_two_squared: 0}).subs(sp.Function("F_2")(0, 0), sp.Function("F_1")(0, 0))

    normalization = derive_normalization()
    velocity = sp.Symbol("sigma", positive=True)
    cosine = sp.Symbol("z", real=True)
    cut_integrand = normalization["Two-photon cut angular integrand"]
    # Reintegrate the angular function with physical real logarithms. The
    # derivative check below relates this primitive to the actual integrand.
    cut_primitive = sp.integrate(cut_integrand, cosine)
    real_primitive = cut_primitive.xreplace({
        sp.log(2 * velocity * cosine - 2): sp.log(1 - velocity * cosine)})
    reduced_cut = sp.simplify(real_primitive.subs(cosine, 1)
                              - real_primitive.subs(cosine, -1))
    coefficient_cut = sp.simplify(
        2 * alpha**2 * lepton_mass * real_photon_form_factor * reduced_cut)

    loop_coefficient, local_coefficient = sp.symbols("g_l_TFF g_l_local")
    spin_sum = normalization["Lepton spin sum divided by |g_all|^2"]
    width = sp.simplify(spin_sum * velocity / (16 * sp.pi * sp.sqrt(
        invariant_mass_squared)) * sp.Abs(loop_coefficient + local_coefficient)**2)
    interference = sp.expand_complex((loop_coefficient + local_coefficient)
        * sp.conjugate(loop_coefficient + local_coefficient)).expand()

    dimension = sp.Symbol("d", positive=True)
    ultraviolet_squared = sp.Symbol("K_squared", positive=True)
    ultraviolet_numerator_average = sp.simplify(numerator.subs(
        scalar_product**2, invariant_mass_squared * photon_one_squared / dimension))
    ultraviolet_integrand_minkowski = sp.simplify(unnormalized_integrand.subs({
        scalar_product**2: invariant_mass_squared * photon_one_squared / dimension,
        photon_two_squared: photon_one_squared,
        lepton_denominator: photon_one_squared,
        full_form_factor: real_photon_form_factor}).subs(dimension, 4))
    # Wick rotation d^4k_M=i d^4K_E and k_M^2=-K_E^2. At leading
    # logarithmic order external momenta/masses are negligible in denominators.
    ultraviolet_integrand_euclidean = sp.simplify(sp.I
        * ultraviolet_integrand_minkowski.subs(photon_one_squared, -ultraviolet_squared))
    radius, lower_scale, cutoff = sp.symbols("K mu Lambda", positive=True)
    sphere_area = sp.simplify(2 * sp.pi**(sp.Rational(4, 2)) / sp.gamma(sp.Rational(4, 2)))
    radial_integrand = sp.simplify(sphere_area * radius**3
        * ultraviolet_integrand_euclidean.subs(ultraviolet_squared, radius**2))
    ultraviolet_shell = sp.integrate(radial_integrand, (radius, lower_scale, cutoff))
    local_cutoff_slope = sp.simplify(-sp.diff(ultraviolet_shell, cutoff) * cutoff)
    falloff_power = sp.Symbol("n", positive=True)
    convergent_tail = sp.integrate(radius**(-1-falloff_power), (radius, lower_scale, sp.oo))

    return {
        "Reduced normalized loop integrand (integrate d^4k with Feynman i0)": reduced_integrand,
        "Unnormalized g_l TFF integrand": unnormalized_integrand,
        "Source-resolved loop integrand": source_resolved_integrand,
        "Source-linearity difference": source_linearity_difference,
        "Finite integrand for responses cancelling only at zero": cancellation_integrand,
        "Value at origin when F_1(0,0)=F_2(0,0)": cancellation_at_origin,
        "Two-real-photon angular integrand": cut_integrand,
        "Real angular primitive on 0<sigma<1": real_primitive,
        "Primitive derivative difference": sp.simplify(sp.diff(real_primitive, cosine) - cut_integrand),
        "Unnormalized two-real-photon cut coefficient (Disc divided by 2i)": coefficient_cut,
        "Threshold reduced cut": sp.limit(reduced_cut, velocity, 0, dir="+"),
        "Spin sum": spin_sum,
        "Coherent total width": width,
        "Coherent modulus squared expanded": interference,
        "Lorentz-symmetric hard numerator average": ultraviolet_numerator_average,
        "Leading pointlike Minkowski g_l integrand": ultraviolet_integrand_minkowski,
        "Leading pointlike Euclidean g_l integrand": ultraviolet_integrand_euclidean,
        "Four-dimensional Euclidean sphere area": sphere_area,
        "Pointlike radial integrand": radial_integrand,
        "Pointlike UV shell between mu and Lambda": ultraviolet_shell,
        "Required local derivative with respect to log Lambda": local_cutoff_slope,
        "Radial tail for falloff K^(-n), n>0": convergent_tail,
        "Two-real-photon massless helicity limit": sp.limit(
            lepton_mass * sp.log(lepton_mass**2 / invariant_mass_squared), lepton_mass, 0, dir="+"),
    }


class UnnormalizedLoopTests(unittest.TestCase):
    """Independent checks of cancellation, physical cuts and UV scaling."""

    @classmethod
    def setUpClass(cls):
        cls.results = derive()

    def test_no_real_photon_normalization_division(self):
        expression = self.results["Unnormalized g_l TFF integrand"]
        self.assertNotIn(sp.Symbol("F_agg"), expression.free_symbols)
        self.assertEqual(self.results["Source-linearity difference"], 0)
        self.assertEqual(self.results["Value at origin when F_1(0,0)=F_2(0,0)"], 0)
        self.assertNotEqual(self.results["Finite integrand for responses cancelling only at zero"], 0)

    def test_cut_integral_and_threshold(self):
        self.assertEqual(self.results["Primitive derivative difference"], 0)
        self.assertEqual(self.results["Threshold reduced cut"], -sp.pi)
        sigma = sp.Symbol("sigma", positive=True)
        alpha, mass = sp.symbols("alpha_em m_l", positive=True)
        form_factor = sp.Symbol("F_agg")
        expected = sp.pi * alpha**2 * mass * form_factor / sigma * sp.log((1-sigma)/(1+sigma))
        self.assertEqual(sp.simplify(self.results[
            "Unnormalized two-real-photon cut coefficient (Disc divided by 2i)"] - expected), 0)

    def test_width_interference_and_helicity(self):
        s, velocity = sp.symbols("s sigma", positive=True)
        loop, local = sp.symbols("g_l_TFF g_l_local")
        expected = sp.sqrt(s) * velocity / (8*sp.pi) * sp.Abs(loop+local)**2
        self.assertEqual(sp.simplify(self.results["Coherent total width"] - expected), 0)
        self.assertEqual(self.results["Two-real-photon massless helicity limit"], 0)

    def test_pointlike_log_and_local_cancellation(self):
        alpha, mass, cutoff, scale = sp.symbols("alpha_em m_l Lambda mu", positive=True)
        form_factor = sp.Symbol("F_agg")
        expected = -6*alpha**2*mass*form_factor*sp.log(cutoff/scale)
        self.assertEqual(sp.simplify(sp.expand_log(self.results[
            "Pointlike UV shell between mu and Lambda"] - expected, force=True)), 0)
        self.assertEqual(sp.simplify(self.results[
            "Required local derivative with respect to log Lambda"] - 6*alpha**2*mass*form_factor), 0)
        falloff, scale = sp.symbols("n mu", positive=True)
        self.assertEqual(sp.simplify(self.results[
            "Radial tail for falloff K^(-n), n>0"] - scale**(-falloff)/falloff), 0)


def main():
    """Print reproducible intermediate expressions and final relations."""
    print("Natural units GeV; (+---); epsilon^0123=+1; leading QED.")
    print("Loop-coefficient convention: g_l|TFF=2 alpha_em^2 m_l F_agg A_l.")
    print("Its unnormalized integral is defined directly when F_agg=0.")
    print("Pointlike UV shell uses a Euclidean radial cutoff, not a finite full amplitude.")
    for label, expression in derive().items():
        print(f"{label}:\n{sp.latex(expression)}\n")


if __name__ == "__main__":
    if "--test" in sys.argv:
        unittest.main(argv=[sys.argv[0]])
    else:
        main()
