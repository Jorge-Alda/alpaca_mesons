#!/usr/bin/env python3
"""Derive the running mixed inverse propagator; GeV, s in GeV^2.

Run: .venv/bin/python scripts/derive_mixed_running_block.py
One loop, first order in 1/f_a, arbitrary flavour matrices.
"""
import unittest
import sympy as sp
import derive_nlo_loop_self_energies as loop
from derive_nlo_quadratic_chipt import derive
from derive_mixed_uv_field_renormalization import meson_transformation


def derive_block():
    derivation = derive()
    symbols = derivation['symbols']
    axial = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f'cA_{i}{j}'))
    pseudoscalar = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f'cP_{i}{j}'))
    gluonic = symbols['c_G']
    constants = {name: sp.Symbol(name + '_r') for name in loop.LOCAL_CONSTANTS}
    chi = loop.chi_zero()
    local = []
    leading = []
    for label in loop.NONET_LABELS:
        generator = loop.GENERATORS[label]
        matrices = {'Phi': generator, 'partial_mu_Phi': generator, 'chi_0': chi,
                    'c_q^A': axial, 'c_q^P': pseudoscalar,
                    'c_q^S': sp.zeros(3), 'c_q^V': sp.zeros(3)}
        expression = 0
        for name, result in derivation['results'].items():
            if name not in constants:
                continue
            term = result['aPhi']
            term = term.xreplace({m: matrices[str(m)] for m in term.atoms(sp.MatrixSymbol)}).doit()
            expression += constants[name] * term
        expression = sp.expand(expression.subs({symbols['f_pi']: loop.F_PI, symbols['f_a']: loop.F_A}))
        local.append(sp.expand(loop.MANDELSTAM_S * expression.coeff(symbols['derivative_a']) + expression.coeff(symbols['a'])))
        leading.append(sp.expand(-loop.F_PI*loop.MANDELSTAM_S/loop.F_A*(axial*generator).trace()
            + loop.F_PI/(2*loop.F_A)*(pseudoscalar*(chi*generator+generator*chi)).trace()
            - loop.F_PI*loop.SINGLET_MASS**2/(3*loop.F_A)*gluonic*generator.trace()))
    return symbols, axial, pseudoscalar, constants, sp.Matrix(leading), sp.Matrix(local)


def scale_terms():
    symbols, axial, pseudoscalar, constants, leading, local = derive_block()
    gamma = loop.gamma_table('nonet')
    local_running = local.diff(next(iter(constants.values()))) * 0
    for name, constant in constants.items():
        factor = gamma[name]
        if name.startswith('Lambda'):
            factor *= loop.SINGLET_MASS**2 / loop.F_PI**2
        local_running -= local.diff(constant) * factor
    b_factor = loop.B_POLE_FROM_RECORD*loop.SINGLET_MASS**2/loop.F_PI**2
    b_running = -b_factor * loop.B_ZERO * leading.diff(loop.B_ZERO)
    pole_axial, pole_pseudoscalar = loop.mixed_tadpole_coefficients(2*loop.leading_mass_matrix(loop.NONET_LABELS))
    loop_running = sp.Matrix([-loop.MANDELSTAM_S*(axial*a).trace()-(pseudoscalar*p).trace()
                              for a,p in zip(pole_axial,pole_pseudoscalar)])
    induced = meson_transformation().T * leading
    return leading, local, local_running, b_running, loop_running, induced


class Checks(unittest.TestCase):
    def test_scale_cancellation(self):
        leading, local, running, b_running, tadpole, induced = scale_terms()
        for entry in running+b_running+tadpole-induced:
            self.assertEqual(sp.expand(entry), 0)


def main():
    leading, local, running, b_running, tadpole, induced = scale_terms()
    print('Q=s K-M^2+Sigma; GeV, s in GeV^2; first order in 1/f_a.')
    print('Each entry of the local mixed block, ordered 0,...,8:')
    for entry in local:
        print(sp.latex(entry))
    print('Scale variations per log(mu/mu0)/(16 pi^2):')
    for label, block in [('local',running),('B',b_running),('tadpole',tadpole),('finite field transformation',-induced)]:
        print(label, sp.latex(block))
    unittest.main(argv=[__file__], verbosity=2)

if __name__ == '__main__':
    main()
