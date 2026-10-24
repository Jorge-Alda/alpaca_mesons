#!/usr/bin/env python
"""Derive charged-pion P-wave angular and two-body-cut normalizations.

Natural units, metric (+---), epsilon^{0123}=+1.  States have relativistic
normalization and pi+ pi- are distinct particles, so no symmetry factor occurs.
Run with .venv/bin/python scripts/derive_alp_pion_partial_wave_unitarity.py
or append --test for the independent projector and phase-space unit tests.
The scalar used in the cut is the charge-stripped current scalar, with its
phase specified separately when converting a decay amplitude.
Messerli:2025rnv Eqs. (3.3)--(3.5) use this current P wave and a symmetric
double-spectral density followed by the explicit exchanged ordering.  Its
1/(192 pi) coefficient is half the single-cut 1/(96 pi) coefficient; it
does not rescale the production P wave.
"""
from __future__ import annotations

import argparse
import unittest

import sympy as sp


pion_angle = sp.Symbol('z', real=True)
azimuth = sp.Symbol('phi', real=True)
dipion_mass_squared, pion_velocity = sp.symbols('s sigma_pi', positive=True)
current_momentum = sp.Symbol('k_z', nonzero=True, real=True)
electric_charge = sp.Symbol('e', positive=True)
pion_momentum = sp.sqrt(dipion_mass_squared) * pion_velocity / 2


def derivation() -> dict[str, sp.Expr]:
    """Integrate the transverse tensor component in the dipion rest frame."""
    # Choose mu=1 for the vector current and nu=2 for the second current;
    # stripping the common epsilon(1,2,q1,k) sign leaves the following ratio.
    pion_direction_x = sp.sqrt(1 - pion_angle**2) * sp.cos(azimuth)
    pion_vector_current_x = 2 * pion_momentum * pion_direction_x
    production_tensor_y = (
        sp.sqrt(dipion_mass_squared) * pion_momentum
        * current_momentum * pion_direction_x
    )
    tff_tensor_xy = sp.sqrt(dipion_mass_squared) * current_momentum
    phase_space_density = pion_velocity / (32 * sp.pi**2)
    angular_cut = sp.simplify(sp.integrate(
        phase_space_density * pion_vector_current_x
        * production_tensor_y / tff_tensor_xy,
        (azimuth, 0, 2 * sp.pi),
    ))
    projector_coefficient = 1 / sp.integrate(
        1 - pion_angle**2, (pion_angle, -1, 1)
    )
    # 2 Im T = sum_n int dPhi_n T* T.  The factor 1/2 is unitarity,
    # not an identical-pion symmetry factor.
    single_cut_coefficient = sp.simplify(
        angular_cut / (1 - pion_angle**2) / projector_coefficient / 2
    )
    quoted_density_coefficient = (
        dipion_mass_squared * pion_velocity**3 / (192 * sp.pi)
    )
    return {
        'dPhi2/(dphi dz)': phase_space_density,
        'azimuth-integrated transverse contraction': angular_cut,
        'P-wave projector normalization': projector_coefficient,
        'Im_s F coefficient multiplying Fpi* f_current': single_cut_coefficient,
        'quoted double-spectral coefficient': quoted_density_coefficient,
        'double-spectral per-ordering allocation': sp.simplify(
            quoted_density_coefficient / single_cut_coefficient
        ),
        'Messerli current-wave identification calF/f_current': sp.simplify(
            single_cut_coefficient / 2 / quoted_density_coefficient
        ),
    }


class PionPartialWaveTests(unittest.TestCase):
    """Check angular orthogonality and independently integrated phase space."""

    def test_derivative_legendre_orthogonality(self) -> None:
        for angular_momentum in range(1, 7):
            projected = sp.Rational(3, 4) * sp.integrate(
                (1 - pion_angle**2)
                * sp.diff(sp.legendre(angular_momentum, pion_angle), pion_angle),
                (pion_angle, -1, 1),
            )
            self.assertEqual(projected, int(angular_momentum == 1))

    def test_phase_space_total(self) -> None:
        integrated = sp.integrate(
            pion_velocity / (32 * sp.pi**2),
            (azimuth, 0, 2 * sp.pi), (pion_angle, -1, 1),
        )
        self.assertEqual(integrated, pion_velocity / (8 * sp.pi))

    def test_tensor_average(self) -> None:
        # Isotropic tensor integral independently supplies 1/3 of |p|^2.
        isotropic_result = sp.simplify(
            pion_velocity / (8 * sp.pi) * 2 * pion_momentum**2 / 3 / 2
        )
        self.assertEqual(
            sp.simplify(derivation()[
                'Im_s F coefficient multiplying Fpi* f_current'
            ] - isotropic_result), 0,
        )

    def test_constant_and_p3_scalar(self) -> None:
        constant, third_wave = sp.symbols('constant third_wave')
        amplitude = constant + third_wave * sp.diff(
            sp.legendre(3, pion_angle), pion_angle
        )
        projected = sp.Rational(3, 4) * sp.integrate(
            (1 - pion_angle**2) * amplitude, (pion_angle, -1, 1)
        )
        self.assertEqual(sp.expand(projected), constant)

    def test_symmetric_double_spectral_allocation(self) -> None:
        results = derivation()
        self.assertEqual(
            2 * results['quoted double-spectral coefficient'],
            results['Im_s F coefficient multiplying Fpi* f_current'],
        )
        self.assertEqual(
            results['Messerli current-wave identification calF/f_current'], 1,
        )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--test', action='store_true')
    arguments = parser.parse_args()
    if arguments.test:
        unittest.main(argv=['derive_alp_pion_partial_wave_unitarity.py'])
    else:
        for label, expression in derivation().items():
            print(f'{label}:\n{sp.latex(expression)}\n')
