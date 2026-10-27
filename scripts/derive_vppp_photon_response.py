"""Derive the tree photon-source response of the local isoscalar VPPP vertex.

Run: .venv/bin/python scripts/derive_vppp_photon_response.py
Tests: .venv/bin/python -m unittest discover -s scripts -p derive_vppp_photon_response.py

Natural units (GeV), exact isospin, canonical tensor kinetic term, leading
cubic pion order, on-shell equal-mass pions, and no vector self energies.
The source is v_mu=-e Q A_mu. Scalars divide out the same common
i epsilon(p_plus,p_minus,p_zero,epsilon_gamma) as the projection script.
The e-stripped photon scalar is defined by division by e; no independent
external-current phase is introduced. This is a local resonance contribution,
excluding WZW, vector exchange into two pions, rescattering and source contacts.

Inputs: memory/rcht/odd-parity/isoscalar-vppp-three-pion-projection.md;
memory/rcht/odd-parity/vector-stage4-amplitudes.md;
memory/rcht/odd-parity/levi-civita-convention.md.
"""

import unittest

import sympy as sp

if __package__:
    from . import derive_vppp_isoscalar_three_pion as projection
else:
    import derive_vppp_isoscalar_three_pion as projection


def derive_photon_response():
    """Vary the quadratic tensor expression and insert its sourced solution.

    The quadratic expression has F_gamma squared / 4 removed. Its solution
    multiplies q^mu epsilon_gamma^nu-q^nu epsilon_gamma^mu. The previously
    derived external vector projection instead includes 1/M_i; multiplication
    by M_i removes that external-state factor before the source is inserted.
    All field and coupling dependence is kept symbolic.
    """
    electric_charge, source_residue = sp.symbols('e R_i', real=True)
    response_coefficient = sp.Symbol('r_i', real=True)
    quadratic_expression = (
        (projection.vector_mass**2-projection.total_invariant)*response_coefficient**2
        -2*sp.sqrt(2)*electric_charge*source_residue*response_coefficient)
    response = sp.solve(sp.diff(quadratic_expression, response_coefficient),
                        response_coefficient)[0]
    _, _, external_vector_scalar, _ = projection.derive_results()
    tensor_vertex_scalar = sp.simplify(projection.vector_mass*external_vector_scalar)
    photon_scalar = sp.factor(tensor_vertex_scalar*response)
    return quadratic_expression, response, tensor_vertex_scalar, photon_scalar


def derive_isoscalar_sum():
    """Assemble both physical isoscalar poles without changing their residues."""
    _, _, _, photon_scalar = derive_photon_response()
    omega_mass, phi_mass = sp.symbols('M_omega M_phi', positive=True)
    omega_residue, phi_residue = sp.symbols('R_omega R_phi', real=True)
    mixing_angle = sp.Symbol('theta_V', real=True)
    source_residue = sp.Symbol('R_i', real=True)
    omega_scalar = photon_scalar.subs({projection.vector_mass: omega_mass,
        source_residue: omega_residue,
        projection.non_strange_entry: sp.cos(mixing_angle)/sp.sqrt(2)})
    phi_scalar = photon_scalar.subs({projection.vector_mass: phi_mass,
        source_residue: phi_residue,
        projection.non_strange_entry: sp.sin(mixing_angle)/sp.sqrt(2)})
    return omega_scalar+phi_scalar


def derive_contact_degeneracy():
    """Separate one bare pole and test its dependence on a free source contact.

    Differentiating with respect to literal g_1 and g_5 varies exactly the
    two combinations of interest, while the other g_i remain fixed. The rows
    are the constant and pole-residue coefficients, which are independent
    functions of energy for a single bare pole. A free constant is a declared
    extra source term, not an additional interaction inferred from VPPP.
    """
    _, _, external_vector, _ = projection.derive_results()
    polynomial = sp.expand(external_vector*projection.vector_mass
        *projection.pion_decay_constant**3*projection.basis_mass
        /(16*projection.non_strange_entry))
    rational_expression = polynomial/(projection.vector_mass**2-projection.total_invariant)
    separated_expression = sp.apart(rational_expression, projection.total_invariant)
    pole_residue = polynomial.subs(projection.total_invariant, projection.vector_mass**2)
    constant = sp.simplify(separated_expression-pole_residue/(
        projection.vector_mass**2-projection.total_invariant))
    free_contact = sp.Symbol('C_contact', real=True)
    coefficient_matrix = sp.Matrix([constant, pole_residue]).jacobian(
        [projection.couplings[0], projection.couplings[4]])
    matrix_with_contact = sp.Matrix([constant+free_contact, pole_residue]).jacobian(
        [projection.couplings[0], projection.couplings[4], free_contact])
    # The identifiable coefficient dimension after allowing the contact to vary.
    contact_column = matrix_with_contact[:, 2:3]
    projected_rank = matrix_with_contact.rank()-contact_column.rank()
    return (rational_expression, separated_expression, constant, pole_residue,
            coefficient_matrix, matrix_with_contact, projected_rank)


class PhotonResponseTests(unittest.TestCase):
    """Check pole normalization, charge counting and coefficient sensitivity."""

    def test_stationarity_and_external_vector_residue(self):
        quadratic, response, tensor_vertex, photon = derive_photon_response()
        response_variable = sp.Symbol('r_i', real=True)
        self.assertEqual(sp.simplify(sp.diff(quadratic, response_variable).subs(
            response_variable, response)), 0)
        electric_charge, residue = sp.symbols('e R_i', real=True)
        _, _, external_vector, _ = projection.derive_results()
        pole_residue = sp.limit((projection.vector_mass**2-projection.total_invariant)
                               *photon, projection.total_invariant,
                               projection.vector_mass**2)
        expected = sp.sqrt(2)*electric_charge*residue*projection.vector_mass*external_vector.subs(
            projection.total_invariant, projection.vector_mass**2)
        self.assertEqual(sp.simplify(pole_residue-expected), 0)
        self.assertEqual(sp.diff(photon, electric_charge, 2), 0)

    def test_ideal_omega_universal_source(self):
        electric_charge, residue, vector_decay_constant = sp.symbols('e R_i F_V', real=True)
        _, _, _, photon = derive_photon_response()
        ideal = photon.subs({projection.non_strange_entry: 1/sp.sqrt(2),
                            residue: vector_decay_constant/(3*sp.sqrt(2))})
        _, _, external_vector, _ = projection.derive_results()
        expected = electric_charge*vector_decay_constant*projection.vector_mass/(
            3*(projection.vector_mass**2-projection.total_invariant))*external_vector.subs(
                projection.non_strange_entry, 1/sp.sqrt(2))
        self.assertEqual(sp.simplify(ideal-expected), 0)

    def test_rank_and_ward_identity(self):
        _, _, _, photon = derive_photon_response()
        coefficient_row = sp.Matrix([[sp.diff(photon, coupling)
                                     for coupling in projection.couplings]])
        for direction in (sp.Matrix([1,1,0,0,0]), sp.Matrix([1,0,1,0,0]),
                          sp.Matrix([0,0,0,1,-2])):
            self.assertEqual(sp.simplify((coefficient_row*direction)[0]), 0)
        # Momentum conservation makes the anomalous three-pion current transverse.
        self.assertEqual(sp.det(sp.Matrix.hstack(*projection.momentum_basis[:3],
                                                projection.total_momentum)), 0)

    def test_free_contact_removes_one_coefficient_constraint(self):
        original, separated, constant, residue, matrix, extended, profiled_rank = derive_contact_degeneracy()
        self.assertEqual(sp.cancel(original-separated), 0)
        self.assertEqual(matrix.rank(), 2)
        self.assertEqual(extended.rank(), 2)
        self.assertEqual(profiled_rank, 1)
        self.assertEqual(matrix.subs(projection.pion_mass_squared, 0).rank(), 1)
        for direction in extended.nullspace():
            self.assertEqual((extended*direction).applyfunc(sp.simplify), sp.zeros(2, 1))


def main():
    """Print the reproducible source insertion and its explicit limitations."""
    print('Tree level; GeV; exact isospin; bare canonical tensor poles.\n')
    quadratic, response, vertex, photon = derive_photon_response()
    electric_charge = sp.Symbol('e', real=True)
    for label, expression in (
        ('Quadratic tensor expression, F_gamma squared / 4 removed', quadratic),
        ('Stationary V_i / F_gamma response', response),
        ('Local tensor vertex without external-state 1/M_i', vertex),
        ('Local photon scalar', photon),
        ('Local photon scalar divided by e', sp.cancel(photon/electric_charge)),
        ('Physical omega plus phi local photon scalar', derive_isoscalar_sum())):
        print(f'{label}:\n{sp.latex(expression)}\n')
    original, separated, constant, residue, matrix, extended, profiled_rank = derive_contact_degeneracy()
    for label, expression in (
        ('One-pole local polynomial divided by M_i squared minus s', original),
        ('One-pole decomposition into constant plus pole', separated),
        ('Non-pole constant', constant),
        ('Pole residue with denominator M_i squared minus s', residue),
        ('Two-combination coefficient matrix (vary literal g_1 and g_5)', matrix),
        ('Two-combination rank with fixed source contact', sp.Integer(matrix.rank())),
        ('Coefficient matrix with independent C_contact', extended),
        ('Total rank including independent source contact', sp.Integer(extended.rank())),
        ('Two-combination rank after allowing source contact to vary', sp.Integer(profiled_rank)),
        ('Null direction in g_1, g_5, C_contact', sp.Matrix.hstack(*extended.nullspace()))):
        print(f'{label}:\n{sp.latex(expression)}\n')
    print('A timelike fit additionally requires a declared vector self-energy, '
          'local-source contact and rescattering prescription. No width is inserted here.')


if __name__ == '__main__':
    main()
