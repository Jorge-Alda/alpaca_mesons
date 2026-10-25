#!/usr/bin/env python
"""Derive the conditional three-pion cut and subtracted dispersion factors.

Natural units GeV; metric (+---), epsilon^{0123}=+1; relativistic states,
S=1+iT. The amplitudes, crossing phases, and subtraction data remain inputs.
Run .venv/bin/python scripts/derive_alp_three_pion_discontinuity.py
or append --test for phase-space, tensor, and dispersion unit tests.
"""
from __future__ import annotations

import argparse
from functools import lru_cache
import unittest

import sympy as sp


def kallen(first: sp.Expr, second: sp.Expr, third: sp.Expr) -> sp.Expr:
    """Return the standard Kallen polynomial using SymPy arithmetic."""
    return sp.expand(first**2 + second**2 + third**2
                     - 2 * (first * second + first * third + second * third))


def epsilon_tensor(first: sp.Matrix, second: sp.Matrix) -> sp.Matrix:
    """Form epsilon^{mu nu rho sigma} first_rho second_sigma.

    SymPy LeviCivita supplies components. Explicit finite Lorentz contraction
    is used because it preserves this repository's component conventions.
    """
    metric = sp.diag(1, -1, -1, -1)
    first_lower, second_lower = metric * first, metric * second
    return sp.Matrix(4, 4, lambda mu, nu: sum(
        sp.LeviCivita(mu, nu, rho, sigma)
        * first_lower[rho] * second_lower[sigma]
        for rho in range(4) for sigma in range(4)
    ))


@lru_cache(maxsize=1)
def derive() -> dict[str, sp.Expr]:
    """Build normalization results from recursive phase space and tensors."""
    transition, transition_adjoint = sp.symbols('T T_dagger', commutative=False)
    unitarity_expansion = sp.expand(
        (1-sp.I*transition_adjoint)*(1+sp.I*transition)-1
    )
    x, u = sp.symbols('x u', positive=True)
    mass_plus_squared, mass_minus_squared, mass_zero_squared = sp.symbols(
        'm_plus_squared m_minus_squared m_zero_squared', nonnegative=True
    )
    z, phi = sp.symbols('z phi', real=True)
    two_body_first = sp.sqrt(kallen(x, u, mass_zero_squared)) / (
        32 * sp.pi**2 * x
    )
    two_body_second = sp.sqrt(kallen(
        u, mass_plus_squared, mass_minus_squared
    )) / (32 * sp.pi**2 * u)
    # Integrate total-system orientation and the pair's internal azimuth.
    recursive_density = sp.simplify(
        two_body_first * two_body_second * 4 * sp.pi * 2 * sp.pi
        / (2 * sp.pi)
    )
    pair_energy_plus = (u + mass_plus_squared - mass_minus_squared) / (
        2 * sp.sqrt(u)
    )
    pair_energy_zero = (x - u - mass_zero_squared) / (2 * sp.sqrt(u))
    pair_momentum_plus = sp.sqrt(kallen(
        u, mass_plus_squared, mass_minus_squared
    )) / (2 * sp.sqrt(u))
    pair_momentum_zero = sp.sqrt(kallen(
        x, u, mass_zero_squared
    )) / (2 * sp.sqrt(u))
    pair_t = (mass_plus_squared + mass_zero_squared
              + 2 * pair_energy_plus * pair_energy_zero
              - 2 * pair_momentum_plus * pair_momentum_zero * z)
    jacobian = -sp.diff(pair_t, z)
    dalitz_density = sp.simplify(recursive_density / jacobian)

    # The geometric current norm is obtained from its three-vector Gram
    # determinant. This finite component identity is verified independently.
    dot_plus_minus, dot_plus_zero, dot_minus_zero = sp.symbols(
        'p_plus_dot_p_minus p_plus_dot_p_zero p_minus_dot_p_zero', real=True
    )
    gram = sp.Matrix([
        [mass_plus_squared, dot_plus_minus, dot_plus_zero],
        [dot_plus_minus, mass_minus_squared, dot_minus_zero],
        [dot_plus_zero, dot_minus_zero, mass_zero_squared],
    ])
    gram_dalitz = sp.factor(gram.det().subs({
        dot_plus_minus: (u - mass_plus_squared - mass_minus_squared) / 2,
        dot_plus_zero: (pair_t - mass_plus_squared - mass_zero_squared) / 2,
        dot_minus_zero: (x + mass_plus_squared - pair_t - u) / 2,
    }))
    direction = sp.Matrix([
        sp.sqrt(1 - z**2) * sp.cos(phi),
        sp.sqrt(1 - z**2) * sp.sin(phi), z,
    ])
    angular_average = sp.Matrix(3, 3, lambda row, col: sp.simplify(
        sp.integrate(direction[row] * direction[col],
                     (phi, 0, 2 * sp.pi), (z, -1, 1)) / (4 * sp.pi)
    ))
    orientation_denominator = sp.symbols('orientation_denominator')
    spatial_trace = sp.trace(sp.eye(3) / orientation_denominator)
    isotropic_coefficient = sp.solve(
        sp.Eq(spatial_trace, 1), orientation_denominator
    )[0]

    metric = sp.diag(1, -1, -1, -1)
    q1_energy, q1_z, q2_energy, q2_z = sp.symbols(
        'q1_0 q1_z q2_0 q2_z', real=True
    )
    q1 = sp.Matrix([q1_energy, 0, 0, q1_z])
    q2 = sp.Matrix([q2_energy, 0, 0, q2_z])
    tensor = epsilon_tensor(q1, q2)
    tensor_lower = metric * tensor * metric
    tensor_norm = sp.expand(sum(
        tensor[row, col] * tensor_lower[row, col]
        for row in range(4) for col in range(4)
    ))
    q1_squared = (q1.T * metric * q1)[0]
    q2_squared = (q2.T * metric * q2)[0]
    alp_mass_squared = ((q1 + q2).T * metric * (q1 + q2))[0]
    norm_from_invariants = sp.expand(
        kallen(alp_mass_squared, q1_squared, q2_squared) / 2
    )
    # Solve a scalar coefficient equation after contraction, never divide by
    # an antisymmetric tensor. Factor i is from S=1+iT unitarity.
    tensor_contracted_current_product, disc_form_factor = sp.symbols(
        'tensor_contracted_current_product disc_F'
    )
    invariant_lambda = sp.Symbol('lambda_ma2_s1_s2', nonzero=True)
    scalar_cut = sp.solve(sp.Eq(
        invariant_lambda * disc_form_factor / 2,
        sp.I * tensor_contracted_current_product,
    ), disc_form_factor)[0]

    s1, s2, spectral_x, spectral_y = sp.symbols('s1 s2 x_spectral y_spectral')
    # Sequential subtraction: (1/(x-s)-1/x) is derived before multiplying
    # the two spectral-variable factors.
    subtracted_first = sp.factor(1 / (spectral_x - s1) - 1 / spectral_x)
    subtracted_second = sp.factor(1 / (spectral_y - s2) - 1 / spectral_y)
    double_denominator = sp.factor(subtracted_first * subtracted_second)
    first_cauchy_residue = sp.residue(subtracted_first, spectral_x, s1)
    crossed_double = double_denominator.xreplace({
        s1: s2, s2: s1, spectral_x: spectral_y, spectral_y: spectral_x,
    })
    pole_mass_squared, residue = sp.symbols('M_i_squared pole_residue')
    pole = residue / (pole_mass_squared - s1)
    pole_subtraction = pole.subs(s1, 0)
    pole_from_cut = sp.factor(residue * subtracted_first.subs(
        spectral_x, pole_mass_squared
    ))
    massless_density = dalitz_density
    massless_total = sp.integrate(massless_density, (sp.Symbol('t'), 0, x-u),
                                  (u, 0, x))
    return {
        'S dagger S minus one for S=1+iT': unitarity_expansion,
        'recursive dPhi3/(du dz), all spatial orientations integrated': recursive_density,
        't in charged-pair rest frame': pair_t,
        'absolute dt/dz': sp.factor(jacobian),
        'dPhi3/(du dt), orientation-independent integrand only': dalitz_density,
        'Gram determinant after Dalitz substitution': gram_dalitz,
        'normalized spatial orientation average': angular_average,
        'spin-one trace denominator': isotropic_coefficient,
        'tensor norm': tensor_norm,
        'tensor norm minus lambda/2': sp.expand(tensor_norm-norm_from_invariants),
        'first photon Ward contraction': (metric*q1).T*tensor,
        'second photon Ward contraction': tensor*(metric*q2),
        'Disc F from contracted current product': scalar_cut,
        'once-subtracted first-variable denominator': subtracted_first,
        'Cauchy pole residue of subtracted denominator': first_cauchy_residue,
        'sequential double-spectral denominator': double_denominator,
        'Bose-exchange denominator residual': sp.simplify(double_denominator-crossed_double),
        'pole subtraction value': pole_subtraction,
        'pole contribution reconstructed from its discontinuity': pole_from_cut,
        'pole reconstruction residual': sp.simplify(pole-pole_subtraction-pole_from_cut),
        'integrated massless three-body phase space': sp.simplify(massless_total),
    }


class ThreePionDiscontinuityTests(unittest.TestCase):
    """Independent normalization, Lorentz, and analytic-continuation tests."""

    def test_massless_phase_space_recursion(self) -> None:
        x, u = sp.symbols('x u', positive=True)
        # Independently integrate two total two-body phase spaces, including
        # the invariant-mass recursion factor 1/(2pi).
        total = sp.integrate((x-u)/(8*sp.pi*x)/(8*sp.pi)/(2*sp.pi), (u,0,x))
        self.assertEqual(total, x/(256*sp.pi**3))
        self.assertEqual(derive()['integrated massless three-body phase space'], total)

    def test_general_tensor_norm_and_wards(self) -> None:
        first = sp.Matrix(sp.symbols('a0:4', real=True))
        second = sp.Matrix(sp.symbols('b0:4', real=True))
        metric = sp.diag(1,-1,-1,-1)
        tensor = epsilon_tensor(first, second)
        lowered = metric*tensor*metric
        norm = sum(tensor[i,j]*lowered[i,j] for i in range(4) for j in range(4))
        scalar = 2*((first.T*metric*second)[0]**2
                    -(first.T*metric*first)[0]*(second.T*metric*second)[0])
        self.assertEqual(sp.expand(norm-scalar), 0)
        self.assertEqual(((metric*first).T*tensor).applyfunc(sp.expand), sp.zeros(1,4))
        self.assertEqual((tensor*(metric*second)).applyfunc(sp.expand), sp.zeros(4,1))

    def test_three_pion_gram_identity(self) -> None:
        momenta = [sp.Matrix(v) for v in ((4,1,0,1),(3,0,1,0),(5,-1,-1,-1))]
        metric = sp.diag(1,-1,-1,-1)
        lowered = [metric*p for p in momenta]
        current = sp.Matrix([sum(
            sp.LeviCivita(mu, a,b,c)*lowered[0][a]*lowered[1][b]*lowered[2][c]
            for a in range(4) for b in range(4) for c in range(4)
        ) for mu in range(4)])
        gram = sp.Matrix(3,3,lambda i,j:(momenta[i].T*metric*momenta[j])[0])
        self.assertEqual(-(current.T*metric*current)[0], gram.det())
        self.assertEqual(((sum(momenta,sp.zeros(4,1))).T*metric*current)[0], 0)

    def test_spin_one_average(self) -> None:
        self.assertEqual(derive()['normalized spatial orientation average'], sp.eye(3)/3)
        self.assertEqual(derive()['spin-one trace denominator'], 3)

    def test_subtracted_pole_and_bose_symmetry(self) -> None:
        self.assertEqual(derive()['pole reconstruction residual'], 0)
        self.assertEqual(derive()['Bose-exchange denominator residual'], 0)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--test', action='store_true')
    arguments = parser.parse_args()
    if arguments.test:
        unittest.main(argv=['derive_alp_three_pion_discontinuity.py'])
    else:
        for label, expression in derive().items():
            print(f'{label}:\n{sp.latex(expression)}\n')
