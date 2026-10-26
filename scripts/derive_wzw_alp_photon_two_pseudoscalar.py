"""Derive the leading WZW ALP photon two-pseudoscalar vertex.

Invocation: .venv/bin/python scripts/derive_wzw_alp_photon_two_pseudoscalar.py
Tests: .venv/bin/python scripts/derive_wzw_alp_photon_two_pseudoscalar.py --test

Natural units, U=exp(i Phi/f_pi), v=-e Q A, epsilon^{0123}=+1.
Retain O(e/f_a), tree level and leading anomalous chiral order. Sources
commute with Q. The WZW orientation follows the existing gamma3pi record
and remains explicit. The source-only anomalous variation under a local
axial rotation has no one-photon/two-meson term. Consequently the pure
axial-gradient source can be removed by Phi -> Phi-2 f_pi c_q^A a/f_a;
only the linear part of that field change enters gamma Phi^3 at this order.
This anomalous Ward identity is an input; the script derives its flavor
coefficients, LSZ combination and redefinition cancellation, not the full
gauged WZW functional from its five-dimensional definition.

Inputs: memory/rcht/odd-parity/wzw-gamma3pi-contribution.md;
memory/extracts/Alda:2026nrd/appendix-a-expansion.md;
memory/rcht/odd-parity/odd-alp-source-lsz-assembly.md.
"""

import argparse
import itertools
import unittest

import sympy as sp
from sympy.combinatorics import Permutation


def show(label, expression):
    """Print a derivation expression using SymPy's LaTeX renderer."""
    print(f"{label}:\n{sp.latex(expression)}\n")


def antisymmetric_trace(charge_matrix, first_matrix, second_matrix, third_matrix):
    """Extract a matrix-valued three-form coefficient by six permutations.

    SymPy has no matrix-valued exterior-product coefficient operation. Its
    matrix multiplication, Trace and permutation signature supply this
    minimal finite expansion without assuming flavor matrices commute.
    """
    matrices = (first_matrix, second_matrix, third_matrix)
    terms = [
        Permutation(order).signature()
        * sp.trace(sp.expand(charge_matrix * matrices[order[0]]
                   * matrices[order[1]] * matrices[order[2]])).doit()
        for order in itertools.permutations(range(3))
    ]
    return sp.Add(*terms).expand()


number_colors, electric_charge, pion_decay_constant, alp_decay_constant = sp.symbols(
    "N_c e f_pi f_a", positive=True, real=True
)
orientation = sp.Symbol("orientation", real=True)
charge_matrix = sp.diag(sp.Rational(2, 3), -sp.Rational(1, 3), -sp.Rational(1, 3))
axial_source = sp.MatrixSymbol("c_q^A", 3, 3)
mass_solution = sp.MatrixSymbol("S", 3, 3)
axial_rotation = sp.MatrixSymbol("kappa", 3, 3)
first_generator = sp.MatrixSymbol("lambda_1^dagger", 3, 3)
second_generator = sp.MatrixSymbol("lambda_2^dagger", 3, 3)
abstract_charge = sp.MatrixSymbol("Q", 3, 3)
pseudoscalar_source = sp.MatrixSymbol("c_q^P", 3, 3)
quark_mass = sp.MatrixSymbol("m_q", 3, 3)
gluon_source, chiral_condensate, singlet_mass_squared = sp.symbols(
    "c_G B_0 m_0^2", real=True
)


def mass_operator_input(axial, pseudoscalar, gluon):
    """Construct Alda's invariant C entering S=M_a^{-1}(C)."""
    difference = axial - pseudoscalar
    return sp.expand(-chiral_condensate * (quark_mass*difference + difference*quark_mass)
        - singlet_mass_squared/3 * (gluon + sp.trace(axial).doit()) * sp.Identity(3))

# The source-linear current product follows directly from dU=i dPhi/f_pi.
# Its two current placements are equal at cubic order.
current_coefficient = 2 * (sp.I / pion_decay_constant)**3
photon_three_meson_coefficient = sp.simplify(
    orientation * number_colors / (48 * sp.pi**2)
    * (-electric_charge) * current_coefficient
)
field_overlap_coefficient = 2 * pion_decay_constant / alp_decay_constant
photon_two_meson_coefficient = sp.simplify(
    photon_three_meson_coefficient * field_overlap_coefficient
)
# Fourier factors: outgoing P1,P2 give +i p1,+i p2; incoming a gives
# -i p_a.  With iM=i integral L, divide the vertex coefficient by i.
# Lowering the four contracted vectors converts epsilon^upper to epsilon_lower.
fourier_derivative_coefficient = sp.I**2 * (-sp.I)
amplitude_coefficient = sp.simplify(
    sp.I * photon_two_meson_coefficient * fourier_derivative_coefficient / sp.I
)


def direct_flavor(source):
    """Return the pure axial-gradient direct flavor coefficient."""
    return -antisymmetric_trace(abstract_charge, first_generator, second_generator, source)


def mixing_flavor(source, solution):
    """Insert the LO ALP-meson overlap in gamma Phi^3."""
    return antisymmetric_trace(abstract_charge, first_generator, second_generator, source + solution)


def combined_flavor(source, solution):
    """Assemble direct and meson-overlap coefficients without commuting them."""
    return sp.expand(direct_flavor(source) + mixing_flavor(source, solution))


# Explicit entries are needed only for the requested physical channel table.
source_up, source_down, source_strange, source_down_strange, source_strange_down = sp.symbols(
    "S_u S_d S_s S_ds S_sd"
)
solution_components = sp.Matrix([
    [source_up, 0, 0], [0, source_down, source_down_strange],
    [0, source_strange_down, source_strange],
])


def elementary_matrix(row, column):
    """Construct a charged-field generator with Alda normalization."""
    return sp.sqrt(2) * sp.eye(3)[:, row] * sp.eye(3)[column, :]


# Final-state outgoing fields require lambda_P^dagger.
outgoing_generators = {
    "pi+": elementary_matrix(1, 0),
    "pi-": elementary_matrix(0, 1),
    "K+": elementary_matrix(2, 0),
    "K-": elementary_matrix(0, 2),
    "K0": elementary_matrix(2, 1),
    "K0bar": elementary_matrix(1, 2),
    "pi0": sp.diag(1, -1, 0),
    "eta8": sp.diag(1, 1, -2) / sp.sqrt(3),
    "eta0": sp.sqrt(sp.Rational(2, 3)) * sp.eye(3),
}


def channel_coefficients():
    """Compute every unordered two-pseudoscalar channel in a U(3) basis."""
    return {
        (first_name, second_name): sp.simplify(antisymmetric_trace(
            charge_matrix, first_matrix, second_matrix, solution_components
        ))
        for (first_name, first_matrix), (second_name, second_matrix)
        in itertools.combinations(outgoing_generators.items(), 2)
    }


class WzwAlpPhotonTests(unittest.TestCase):
    """Independent anomaly normalization and exact basis identities."""

    def test_gamma_three_pion_normalization(self):
        """Recover the existing charged-field anomaly normalization."""
        pion_trace = antisymmetric_trace(
            charge_matrix, outgoing_generators["pi-"],
            outgoing_generators["pi+"], outgoing_generators["pi0"]
        )
        self.assertEqual(pion_trace, 2)
        self.assertEqual(sp.simplify(
            photon_three_meson_coefficient * pion_trace
            - sp.I * orientation * number_colors * electric_charge
            / (12 * sp.pi**2 * pion_decay_constant**3)
        ), 0)

    def test_abstract_basis_cancellation(self):
        """Cancel arbitrary axial sources before choosing flavor entries."""
        expected = antisymmetric_trace(
            abstract_charge, first_generator, second_generator, mass_solution
        )
        self.assertEqual(sp.expand(combined_flavor(axial_source, mass_solution) - expected), 0)

    def test_axial_redefinition(self):
        """An arbitrary axial rotation shifts direct and mixing oppositely."""
        difference = combined_flavor(axial_source + axial_rotation, mass_solution) \
            - combined_flavor(axial_source, mass_solution)
        self.assertEqual(sp.expand(difference), 0)

    def test_no_mixing_basis(self):
        """The hatted axial source -S has a vanishing overlap."""
        self.assertEqual(sp.expand(mixing_flavor(-mass_solution, mass_solution)), 0)
        self.assertEqual(sp.expand(
            direct_flavor(-mass_solution) - combined_flavor(axial_source, mass_solution)
        ), 0)

    def test_mass_input_invariance(self):
        """The same source rotation leaves C, and therefore its resolvent, fixed."""
        original = mass_operator_input(axial_source, pseudoscalar_source, gluon_source)
        transformed = mass_operator_input(
            axial_source + axial_rotation, pseudoscalar_source + axial_rotation,
            gluon_source - sp.trace(axial_rotation)
        )
        self.assertEqual(sp.expand(transformed-original), sp.ZeroMatrix(3, 3))

    def test_flavor_and_bose_checks(self):
        """Antisymmetry, singlet anomaly, and neutral-kaon charge checks."""
        coefficients = channel_coefficients()
        self.assertEqual(coefficients[("pi+", "pi-")], -2*source_up)
        self.assertEqual(coefficients[("K+", "K-")], -2*source_up)
        self.assertEqual(coefficients[("K0", "K0bar")], 2*(source_down-source_strange))
        for generator in outgoing_generators.values():
            self.assertEqual(sp.simplify(antisymmetric_trace(
                charge_matrix, generator, generator, solution_components
            )), 0)
        self.assertEqual(sp.simplify(coefficients[("pi+", "pi-")].subs({
            source_up:1, source_down:1, source_strange:1
        })), -2)

    def test_fourier_sign_and_ward_identity(self):
        """Check the canonical tensor routing and the photon Ward identity."""
        self.assertEqual(sp.simplify(amplitude_coefficient + orientation
            * number_colors * electric_charge
            / (12 * sp.pi**2 * alp_decay_constant * pion_decay_constant**2)), 0)
        photon = sp.Matrix(1, 4, sp.symbols("k_0:4"))
        first_momentum = sp.Matrix(1, 4, sp.symbols("p1_0:4"))
        second_momentum = sp.Matrix(1, 4, sp.symbols("p2_0:4"))
        alp_momentum = photon + first_momentum + second_momentum
        self.assertEqual(sp.simplify(sp.Matrix.vstack(
            photon, first_momentum, second_momentum, alp_momentum
        ).det()), 0)


def derive():
    """Print the executable derivation in field and channel order."""
    print("Natural units; O(e/f_a), tree level, leading WZW order; [Q,c_q^A]=[Q,S]=0.")
    print("The WZW orientation is explicit; displayed coefficients multiply epsilon A dP1 dP2 da.")
    show("Cubic current coefficient 2(i/f_pi)^3", current_coefficient)
    show("gamma Phi^3 coefficient", photon_three_meson_coefficient)
    show("LO meson overlap Phi/a", field_overlap_coefficient * (axial_source + mass_solution))
    show("Common a gamma PP coefficient", photon_two_meson_coefficient)
    show("Fourier derivative coefficient (+i)^2(-i)", fourier_derivative_coefficient)
    show("Common amplitude coefficient in epsilon_lower epsilon_gamma p1 p2 k", amplitude_coefficient)
    show("Direct flavor coefficient", direct_flavor(axial_source))
    show("Mixing flavor coefficient", mixing_flavor(axial_source, mass_solution))
    show("Combined flavor coefficient", combined_flavor(axial_source, mass_solution))
    show("Chiral-redefinition difference", sp.expand(
        combined_flavor(axial_source + axial_rotation, mass_solution)
        - combined_flavor(axial_source, mass_solution)
    ))
    show("Invariant mass-operator input C", mass_operator_input(
        axial_source, pseudoscalar_source, gluon_source
    ))
    show("Mass-operator input C redefinition difference", sp.expand(
        mass_operator_input(axial_source + axial_rotation,
            pseudoscalar_source + axial_rotation, gluon_source - sp.trace(axial_rotation))
        - mass_operator_input(axial_source, pseudoscalar_source, gluon_source)
    ))
    show("No-mixing direct coefficient", direct_flavor(-mass_solution))
    print("Nonzero outgoing-generator coefficients (multiply the common coefficient):")
    for channel, coefficient in channel_coefficients().items():
        if coefficient != 0:
            show(" ".join(channel), coefficient)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--test", action="store_true")
    arguments = parser.parse_args()
    if arguments.test:
        unittest.main(argv=[__file__])
    else:
        derive()
