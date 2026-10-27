"""Conditional VPPP matching to the Hoferichter subtraction normalization.

Run: .venv/bin/python scripts/derive_vppp_subtraction_matching.py
Tests: .venv/bin/python -m unittest discover -s scripts -p derive_vppp_subtraction_matching.py

Inputs are the local projection and photon response scripts and the indexed
Hoferichter:2023bjm/Hoferichter:2025lcz extracts. Natural units, GeV. This
diagnostic assumes exact isospin and bare vector poles. It chooses the common
WZW orientation so the Hoferichter current scalar has positive WZW F_3pi,
and matches a Dalitz-independent local full photon scalar to Delta a equal
to that scalar divided by 3e BEFORE KT rescattering. This matching choice is
an assumption, not a determination of the full amplitude from tree RChT.
The spectral pole is zero-width; no conformal completion or fit is supplied.
"""

import unittest

import sympy as sp

if __package__:
    from . import derive_vppp_photon_response as source
else:
    import derive_vppp_photon_response as source

projection = source.projection


def derive_single_pole_matching():
    """Extract the zero-virtuality subtraction and spectral residue by limits."""
    electric_charge = sp.Symbol('e', real=True)
    _, _, _, photon_scalar = source.derive_photon_response()
    subtraction_shift = sp.cancel(photon_scalar/(3*electric_charge))
    anomaly_shift = sp.simplify(subtraction_shift.subs(projection.total_invariant, 0))
    spectral_residue_shift = sp.simplify(sp.limit(
        (projection.vector_mass**2-projection.total_invariant)*subtraction_shift,
        projection.total_invariant, projection.vector_mass**2))
    subtracted_pole = (spectral_residue_shift*projection.total_invariant/(
        projection.vector_mass**2*(projection.vector_mass**2-projection.total_invariant)))
    sum_rule_mismatch = sp.simplify(anomaly_shift-spectral_residue_shift/projection.vector_mass**2)
    coefficient_matrix = sp.Matrix([anomaly_shift, spectral_residue_shift]).jacobian(
        [projection.couplings[0], projection.couplings[4]])
    return (subtraction_shift, anomaly_shift, spectral_residue_shift,
            subtracted_pole, sum_rule_mismatch, coefficient_matrix)


def derive_two_pole_matching():
    """Sum physical omega and phi sources with their established generators."""
    _, anomaly, residue, _, mismatch, _ = derive_single_pole_matching()
    omega_mass, phi_mass = sp.symbols('M_omega M_phi', positive=True)
    omega_residue, phi_residue = sp.symbols('R_omega R_phi', real=True)
    angle = sp.Symbol('theta_V', real=True)
    omega_substitution = {projection.vector_mass: omega_mass,
        sp.Symbol('R_i', real=True): omega_residue,
        projection.non_strange_entry: sp.cos(angle)/sp.sqrt(2)}
    phi_substitution = {projection.vector_mass: phi_mass,
        sp.Symbol('R_i', real=True): phi_residue,
        projection.non_strange_entry: sp.sin(angle)/sp.sqrt(2)}
    observables = sp.Matrix([anomaly.subs(omega_substitution)+anomaly.subs(phi_substitution),
        residue.subs(omega_substitution), residue.subs(phi_substitution)])
    matrix = observables.jacobian([projection.couplings[0], projection.couplings[4]])
    return (observables, matrix, sp.simplify(mismatch.subs(omega_substitution)
        +mismatch.subs(phi_substitution)), sp.factor(matrix[1:3, :].det()))


class SubtractionMatchingTests(unittest.TestCase):
    """Check exact dispersion identities and physically relevant rank limits."""

    def test_subtracted_dispersion_and_sum_rule(self):
        full, anomaly, _, pole, mismatch, _ = derive_single_pole_matching()
        self.assertEqual(sp.cancel(full-anomaly-pole), 0)
        self.assertEqual(sp.simplify(sp.limit(full, projection.total_invariant, sp.oo)-mismatch), 0)
        self.assertEqual(sp.diff(mismatch, projection.couplings[4]), 0)

    def test_ideal_omega_rank(self):
        *_, matrix = derive_single_pole_matching()
        ideal_matrix = matrix.subs(projection.non_strange_entry, 1/sp.sqrt(2))
        self.assertEqual(ideal_matrix.rank(), 2)
        self.assertEqual(ideal_matrix[1:2, :].rank(), 1)
        self.assertEqual(ideal_matrix.subs(projection.pion_mass_squared, 0).rank(), 1)

    def test_nonideal_pole_residue_rank(self):
        _, matrix, _, determinant = derive_two_pole_matching()
        self.assertNotEqual(determinant, 0)
        self.assertEqual(matrix.rank(), 2)
        self.assertEqual(matrix[1:3, :].rank(), 2)
        omega_mass, phi_mass = sp.symbols('M_omega M_phi', positive=True)
        self.assertEqual(matrix[1:3, :].subs(phi_mass, omega_mass).rank(), 1)
        self.assertEqual(matrix[1:3, :].subs(sp.Symbol('theta_V', real=True), 0).rank(), 1)


def main():
    """Print the matching assumptions, spectral map and rank diagnostics."""
    print('Conditional subtraction-seed matching BEFORE KT rescattering; GeV.\n'
          'Positive WZW F_3pi in the Hoferichter current convention; '
          'Delta a = local full photon scalar / (3e).\n'
          'Exact isospin; bare zero-width vector poles; subtraction at q^2=0.\n')
    full, anomaly, residue, pole, mismatch, matrix = derive_single_pole_matching()
    observables, two_pole_matrix, two_pole_mismatch, determinant = derive_two_pole_matching()
    for label, expression in (
        ('Single-vector Delta a(q squared)', full),
        ('Delta alpha_A from Delta a(0)', anomaly),
        ('Delta c_i: residue for denominator M_i squared minus q squared', residue),
        ('Once-subtracted zero-width spectral pole', pole),
        ('Delta alpha_A minus Delta c_i / M_i squared', mismatch),
        ('Jacobian of Delta alpha_A, Delta c_i with literal g_1 and g_5', matrix),
        ('Ideal omega rank including subtraction', sp.Integer(matrix.rank())),
        ('Ideal omega rank using pole residue alone', sp.Integer(matrix[1:2, :].rank())),
        ('Nonideal Delta alpha_A, Delta c_omega, Delta c_phi', observables),
        ('Nonideal coefficient matrix', two_pole_matrix),
        ('Nonideal rank including subtraction', sp.Integer(two_pole_matrix.rank())),
        ('Nonideal pole-residue determinant', determinant),
        ('Nonideal pole-residue rank for distinct masses and nonzero overlaps',
         sp.Integer(two_pole_matrix[1:3, :].rank())),
        ('Two-pole spectral sum-rule mismatch', two_pole_mismatch)):
        print(f'{label}:\n{sp.latex(expression)}\n')
    print('Pole-only saturation requires the displayed mismatch to vanish. '
          'For a nonzero total source overlap this sets g_1-g_2-g_3=0. '
          'Otherwise additional spectral strength or a changed asymptotic '
          'prescription is needed; no inelastic completion is inferred here. '
          'Finite widths require a new dispersive calculation of the pole terms.')


if __name__ == '__main__':
    main()
