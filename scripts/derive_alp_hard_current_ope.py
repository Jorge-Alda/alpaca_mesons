"""Derive the leading parity-odd local electromagnetic-current OPE.

Run .venv/bin/python scripts/derive_alp_hard_current_ope.py [--test].
Natural units; metric (+---), epsilon^{0123}=+1, gamma5=i gamma0 gamma1
gamma2 gamma3. Fourier convention i integral exp(i q.x) T{Jmu(x/2) Jnu(-x/2)}.
The approximation is a massless free hard quark at leading alpha_s, fixed
soft external momentum, and nearly opposite deep-Euclidean photons.
"""

from functools import lru_cache
import sys
import unittest

import sympy as sp
from sympy.physics.matrices import mgamma


def show(label, expression):
    """Print expressions through the standard SymPy LaTeX renderer."""
    print(f"{label}:\n{sp.latex(expression)}\n")


@lru_cache(maxsize=1)
def derive():
    """Project the two Wick contractions onto the axial Clifford basis.

    SymPy supplies its explicit standard Dirac matrices; these are a fixed
    representation of the declared Clifford algebra, not generic commuting
    replacements for abstract flavor or source matrices. The trace projection
    is representation independent. No custom symbolic algebra is used.
    """
    metric = sp.diag(1, -1, -1, -1)
    gamma = tuple(mgamma(index) for index in range(4))
    gamma5 = mgamma(5)
    hard_momentum = sp.Matrix(sp.symbols("q_0 q_1 q_2 q_3", real=True))
    hard_lower = metric * hard_momentum
    hard_squared = (hard_momentum.T * metric * hard_momentum)[0]
    slash = sum((hard_lower[index] * gamma[index] for index in range(4)), sp.zeros(4))

    # The two single contractions have S_F(q)=i slash(q)/(q^2+i0)
    # and S_F(-q)=-S_F(q). Multiplication by the Fourier prefactor i
    # gives minus the ordered gamma difference divided by q^2.
    ordered_difference = tuple(tuple(
        sp.expand(gamma[first] * slash * gamma[second]
                  - gamma[second] * slash * gamma[first])
        for second in range(4)) for first in range(4))
    axial_coefficients = tuple(tuple(tuple(sp.simplify(
        -sp.trace(gamma[axial] * gamma5 * ordered_difference[first][second]) / 4)
        for axial in range(4)) for second in range(4)) for first in range(4))
    fourier_coefficients = tuple(tuple(tuple(sp.simplify(-coefficient / hard_squared)
        for coefficient in axial_coefficients[first][second])
        for second in range(4)) for first in range(4))
    expected_coefficients = tuple(tuple(tuple(sp.simplify(
        2 * sp.I * sum(sp.LeviCivita(first, second, momentum, axial)
                      * hard_lower[momentum] for momentum in range(4)) / hard_squared)
        for axial in range(4)) for second in range(4)) for first in range(4))
    reconstructed = tuple(tuple(sum((axial_coefficients[first][second][axial]
        * metric[axial, axial] * gamma[axial] * gamma5 for axial in range(4)), sp.zeros(4))
        for second in range(4)) for first in range(4))
    hard_ward = tuple(tuple(sp.simplify(sum(hard_lower[first]
        * fourier_coefficients[first][second][axial] for first in range(4)))
        for axial in range(4)) for second in range(4))

    # Derive the scalar coefficient in a nondegenerate ALP rest frame.
    # Lorentz covariance forces the axial matrix element to be parallel to p.
    # The symbol below denotes its time component, not a newly assigned
    # source-dependent decay constant or a numerical matching input.
    alp_mass = sp.Symbol("m_a", positive=True)
    axial_time = sp.Symbol(r"{\langle 0|\bar{q}\gamma^{0}\gamma_{5}q|a\rangle}")
    soft_momentum = sp.Matrix([alp_mass, 0, 0, 0])
    frame_substitution = {hard_momentum[1]: 0, hard_momentum[2]: 0}
    frame_hard = hard_momentum.subs(frame_substitution)
    first_photon = frame_hard + soft_momentum / 2
    second_photon = -frame_hard + soft_momentum / 2
    first_lower, second_lower = metric * first_photon, metric * second_photon
    epsilon_tensor = sp.Matrix(4, 4, lambda first, second: sp.expand(sum(
        sp.LeviCivita(first, second, third, fourth) * first_lower[third] * second_lower[fourth]
        for third in range(4) for fourth in range(4))))
    current_tensor = sp.Matrix(4, 4, lambda first, second:
        sp.simplify(fourier_coefficients[first][second][0].subs(frame_substitution) * axial_time))
    epsilon_lower = metric * epsilon_tensor * metric
    norm = sp.simplify(sum(epsilon_lower[first, second] * epsilon_tensor[first, second]
        for first in range(4) for second in range(4)))
    projection = sp.simplify(sum(epsilon_lower[first, second] * current_tensor[first, second]
        for first in range(4) for second in range(4)) / norm)
    scalar_reconstruction = (current_tensor - epsilon_tensor * projection).applyfunc(sp.simplify)
    first_ward = (first_lower.T * current_tensor).applyfunc(sp.simplify)
    second_ward = (current_tensor * second_lower).applyfunc(sp.simplify)
    charges = sp.diag(sp.Rational(2, 3), -sp.Rational(1, 3), -sp.Rational(1, 3))
    flavor_weights = charges * charges
    return locals()


class HardCurrentTests(unittest.TestCase):
    """Independent Clifford, tensor, and scalar checks of the derivation."""

    @classmethod
    def setUpClass(cls):
        cls.result = derive()

    def test_clifford_and_axial_orientation(self):
        result = self.result
        gamma, metric = result["gamma"], result["metric"]
        for first in range(4):
            for second in range(4):
                self.assertEqual(gamma[first] * gamma[second] + gamma[second] * gamma[first],
                                 2 * metric[first, second] * sp.eye(4))
        self.assertEqual(sp.trace(gamma[0] * gamma[1] * gamma[2] * gamma[3] * result["gamma5"]),
                         -4 * sp.I)

    def test_trace_projection_and_wick_sign(self):
        result = self.result
        for first in range(4):
            for second in range(4):
                self.assertEqual((result["ordered_difference"][first][second]
                    - result["reconstructed"][first][second]).applyfunc(sp.simplify), sp.zeros(4))
                for axial in range(4):
                    self.assertEqual(sp.simplify(result["fourier_coefficients"][first][second][axial]
                        - result["expected_coefficients"][first][second][axial]), 0)

    def test_hard_ward_and_bose(self):
        result = self.result
        for row in result["hard_ward"]:
            self.assertEqual(row, (0, 0, 0, 0))
        reverse = {component: -component for component in result["hard_momentum"]}
        for first in range(4):
            for second in range(4):
                for axial in range(4):
                    self.assertEqual(sp.simplify(result["fourier_coefficients"][second][first][axial]
                        .subs(reverse, simultaneous=True) - result["fourier_coefficients"][first][second][axial]), 0)

    def test_scalar_projection_and_light_flavor_weights(self):
        result = self.result
        self.assertEqual(result["scalar_reconstruction"], sp.zeros(4))
        self.assertEqual(result["first_ward"], sp.zeros(1, 4))
        self.assertEqual(result["second_ward"], sp.zeros(4, 1))
        self.assertEqual(result["flavor_weights"], sp.diag(sp.Rational(4, 9), sp.Rational(1, 9), sp.Rational(1, 9)))
        q0, _, _, q3 = result["hard_momentum"]
        expected = 2 * sp.I * result["axial_time"] / (result["alp_mass"] * (q0**2 - q3**2))
        self.assertEqual(sp.simplify(result["projection"] - expected), 0)


def main():
    result = derive()
    print("Leading alpha_s=0, massless hard propagator, fixed soft p; metric (+---), epsilon0123=+1.")
    print("Current excludes e. Fourier i integral exp(i q.x), insertions at x/2,-x/2;")
    print("S_F(q)=i slash(q)/(q^2+i0), q=(q1-q2)/2, p=q1+q2 fixed and soft.")
    show("hard momentum squared", result["hard_squared"])
    show("gamma5 trace orientation", sp.trace(result["gamma"][0] * result["gamma"][1]
         * result["gamma"][2] * result["gamma"][3] * result["gamma5"]))
    show("ordered Wick gamma difference, mu=1 nu=2", result["ordered_difference"][1][2])
    show("axial Clifford coefficients before Fourier prefactor, mu=1 nu=2",
         sp.Matrix(result["axial_coefficients"][1][2]))
    show("Fourier OPE axial coefficients, mu=1 nu=2", sp.Matrix(result["fourier_coefficients"][1][2]))
    show("hard-current Ward contractions", sp.Matrix(result["hard_ward"]))
    show("light quark charge-squared weights", result["flavor_weights"])
    show("two-photon epsilon tensor in ALP rest frame", result["epsilon_tensor"])
    show("epsilon scalar norm", result["norm"])
    show("one-flavor scalar coefficient from contraction", result["projection"])
    show("scalar reconstruction difference", result["scalar_reconstruction"])
    show("first physical photon Ward contraction after axial matrix element", result["first_ward"])
    show("second physical photon Ward contraction after axial matrix element", result["second_ward"])
    print("Remaining inputs: matched ALP axial-current matrix elements, active-flavor/scale scheme,")
    print("higher orders and operators, unequal-ratio/light-cone limits, local photon term.")


if __name__ == "__main__":
    if "--test" in sys.argv:
        unittest.main(argv=[sys.argv[0]])
    else:
        main()
