#!/usr/bin/env python3
"""Derive mixed UV coefficients under local field transformations.

Run: .venv/bin/python scripts/derive_mixed_uv_field_renormalization.py
Conventions: Q=s K-M^2+Sigma; GeV, s in GeV^2; first order in 1/f_a.
"""
import unittest
import sympy as sp
import derive_nlo_loop_self_energies as loop


def meson_transformation():
    """Return the established meson field redefinition coefficient of lambda(mu)."""
    coefficients = {
        "{chi_0,Phi}": -sp.Rational(1, 2),
        "<chi_0>Phi": -sp.Rational(1, 3),
        "<chi_0 Phi>1": sp.Rational(1, 3),
        "<Phi>chi_0": sp.Rational(1, 3),
    }
    maps = loop.field_redefinition_maps("nonet")
    return sum((value * maps[name] for name, value in coefficients.items()), sp.zeros(9))


def mixed_poles_and_transformation():
    """Compute arbitrary source coefficients before and after the meson redefinition.

    The matrices within traces multiply c_q^A and c_q^P. The separate ALP
    shift is fixed by the coefficient linear in s, then the residual is tested.
    """
    chi = loop.chi_zero()
    mass = loop.leading_mass_matrix(loop.NONET_LABELS)
    pole_a, pole_p = loop.mixed_tadpole_coefficients(2 * mass)
    local_a, local_p = loop.mixed_local_pole_coefficients()
    redefinition = meson_transformation()
    # From u_mu^(1) and chi_+^(a Phi) in quadratic-lagrangian.md.
    leading_a = [-loop.F_PI * loop.MANDELSTAM_S / loop.F_A * loop.GENERATORS[i]
                 for i in loop.NONET_LABELS]
    leading_p = [loop.F_PI / (2 * loop.F_A)
                 * (loop.GENERATORS[i] * chi + chi * loop.GENERATORS[i])
                 for i in loop.NONET_LABELS]
    transformed = {}
    for name, pole, local, leading in (("axial", pole_a, local_a, leading_a),
                                      ("pseudoscalar", pole_p, local_p, leading_p)):
        initial = [(loop.MANDELSTAM_S if name == "axial" else 1) * (p + l)
                   for p, l in zip(pole, local)]
        induced = [sum((redefinition[j, i] * leading[j] for j in range(9)), sp.zeros(3))
                   for i in range(9)]
        after_meson = [(initial[i] + induced[i]).applyfunc(sp.expand) for i in range(9)]
        shift = [(-entry.diff(loop.MANDELSTAM_S)).applyfunc(sp.expand) for entry in after_meson]
        after_shift = [(after_meson[i] + loop.MANDELSTAM_S * shift[i]
                        - sum((mass[i, j] * shift[j] for j in range(9)), sp.zeros(3)))
                       .applyfunc(sp.expand) for i in range(9)]
        transformed[name] = (initial, induced, after_meson, shift, after_shift)
    return transformed


class Checks(unittest.TestCase):
    def test_meson_transformation(self):
        coefficients = loop.gamma_table("nonet") | {
            "field:{chi_0,Phi}": -sp.Rational(1, 2),
            "field:<chi_0>Phi": -sp.Rational(1, 3),
            "field:<chi_0 Phi>1": sp.Rational(1, 3),
            "field:<Phi>chi_0": sp.Rational(1, 3),
        }
        self.assertTrue(loop.total_pole("nonet", coefficients).is_zero_matrix)

    def test_slope_removed(self):
        for data in mixed_poles_and_transformation().values():
            for entry in data[-1]:
                self.assertTrue(entry.diff(loop.MANDELSTAM_S).is_zero_matrix)

    def test_all_mixed_poles_cancel(self):
        for data in mixed_poles_and_transformation().values():
            self.assertTrue(all(entry.is_zero_matrix for entry in data[2]))

    def test_gluonic_mixing_unchanged(self):
        self.assertTrue(meson_transformation().row(0).is_zero_matrix)

    def test_shift_rule(self):
        s, mass_squared, mixed, correction, shift = sp.symbols("s m2 b r h")
        change = sp.Matrix([[1 + correction, shift], [0, 1]])
        inverse_propagator = sp.Matrix([[s - mass_squared, mixed], [mixed, 0]])
        transformed = change.T * inverse_propagator * change
        self.assertEqual(sp.expand(sp.expand(transformed[0, 1]).coeff(correction) - mixed - shift * (s - mass_squared)), 0)
        self.assertEqual(sp.expand(sp.expand(transformed[0, 1]).subs(correction, 0) - mixed - shift * (s - mass_squared)), 0)


def main():
    print("Q=s K-M^2+Sigma; f_pi=93 MeV; lambda coefficient; first order in 1/f_a")
    result = mixed_poles_and_transformation()
    for name, data in result.items():
        print(f"Source {name}: initial pole, induced term, sum, ALP-dependent shift, residual")
        for index in (3, 8, 0):
            print(sp.latex(loop.MESON_FIELDS[index]), *(sp.latex(data[j][index]) for j in (0, 1, 2, 3, 4)), sep="\n")
        print("All nonet residual coefficients vanish:", all(matrix.is_zero_matrix for matrix in data[4]))
    unittest.main(argv=[__file__], verbosity=2)


if __name__ == "__main__":
    main()
