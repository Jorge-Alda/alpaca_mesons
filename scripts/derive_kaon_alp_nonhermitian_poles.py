"""Derive constant-width neutral-kaon–ALP poles; energies are in GeV.

Run: .venv/bin/python scripts/derive_kaon_alp_nonhermitian_poles.py
"""
import sympy as sp
import unittest


def derive():
    """Return symbolic steps, using independent directed mixing entries."""
    lambda_S, lambda_L, lambda_a = sp.symbols('lambda_S lambda_L lambda_a')
    mu_Sa, mu_La, mu_aS, mu_aL = sp.symbols('mu_Sa mu_La mu_aS mu_aL')
    pole, mixing_order = sp.symbols('z t')
    p, q = sp.symbols('p q', nonzero=True)
    epsilon_K = sp.symbols('epsilon_K')
    epsilon_K_bar = sp.symbols('epsilon_K_bar')
    normalization = sp.sqrt(2 * (1 + epsilon_K * epsilon_K_bar))
    flavor_0a, flavor_bara, flavor_a0, flavor_abar = sp.symbols(
        'H_0a H_bara H_a0 H_abar')
    flavor_to_poles = sp.Matrix([[p, p, 0], [q, -q, 0], [0, 0, 1]])
    flavor_mixing = sp.Matrix([[0, 0, flavor_0a], [0, 0, flavor_bara],
                               [flavor_a0, flavor_abar, 0]])
    directed_mixing = sp.simplify(flavor_to_poles.inv() * flavor_mixing * flavor_to_poles)
    # Exact complex epsilon_K convention.  epsilon_K_bar denotes the complex
    # conjugate in the normalization and is kept as a separate symbol.
    cp_odd_flavor_mixing = flavor_mixing.subs({
        flavor_bara: -flavor_0a, flavor_abar: -flavor_a0})
    epsilon_substitutions = {p: (1 + epsilon_K) / normalization,
                             q: (1 - epsilon_K) / normalization}
    cp_odd_pole_mixing = sp.simplify(
        (flavor_to_poles.inv() * cp_odd_flavor_mixing * flavor_to_poles)
        .subs(epsilon_substitutions))
    general_epsilon_mixing = sp.simplify(directed_mixing.subs(epsilon_substitutions))
    hamiltonian = sp.Matrix([[lambda_S, 0, mu_Sa], [0, lambda_L, mu_La],
                             [mu_aS, mu_aL, lambda_a]])
    polynomial = sp.factor((pole * sp.eye(3) - hamiltonian).det())
    kaon_block = pole * sp.eye(2) - hamiltonian[:2, :2]
    secular = (pole - lambda_a
               - (hamiltonian[2, :2] * kaon_block.inv() * hamiltonian[:2, 2])[0])
    two_level = hamiltonian.extract([0, 2], [0, 2])
    two_polynomial = sp.factor((pole * sp.eye(2) - two_level).det())
    roots = sp.solve(two_polynomial, pole)
    splitting_squared = sp.factor((roots[1]-roots[0])**2)
    # A column of the adjugate is a null vector at any root of the determinant.
    two_level_vector = (pole*sp.eye(2)-two_level).adjugate()[:, 1]
    mass, gamma_i, gamma_a, mixing = sp.symbols('m Gamma_i Gamma_a mu', real=True)
    equal_mass_splitting = sp.simplify(splitting_squared.subs({
        lambda_S: mass-sp.I*gamma_i/2, lambda_a: mass-sp.I*gamma_a/2,
        mu_Sa: mixing, mu_aS: mixing}))
    # Introduce t for each off-diagonal element and solve the coefficient of t^2.
    shifted = hamiltonian.subs({mu_Sa: mixing_order*mu_Sa, mu_La: mixing_order*mu_La,
                                 mu_aS: mixing_order*mu_aS, mu_aL: mixing_order*mu_aL})
    perturbation_polynomial = (pole * sp.eye(3) - shifted).det()
    second_order_shift = sp.symbols('delta_lambda')
    shifts = {}
    for name, diagonal in [('S', lambda_S), ('L', lambda_L), ('a', lambda_a)]:
        coefficient = sp.expand(perturbation_polynomial.subs(
            pole, diagonal + mixing_order**2 * second_order_shift)).coeff(mixing_order, 2)
        shifts[name] = sp.factor(sp.solve(coefficient, second_order_shift)[0])
    # Solve (H-lambda)r=0 at first order, fixing the dominant coordinate to one.
    eigenvectors = {}
    first_S, first_L, first_a = sp.symbols('r_S r_L r_a')
    for name, diagonal, vector, unknowns in [
        ('S', lambda_S, sp.Matrix([1, mixing_order*first_L, mixing_order*first_a]), [first_L, first_a]),
        ('L', lambda_L, sp.Matrix([mixing_order*first_S, 1, mixing_order*first_a]), [first_S, first_a]),
        ('a', lambda_a, sp.Matrix([mixing_order*first_S, mixing_order*first_L, 1]), [first_S, first_L])]:
        residual = (shifted - diagonal*sp.eye(3))*vector
        equations = [sp.expand(entry).coeff(mixing_order, 1) for entry in residual]
        solution = sp.solve(equations, unknowns, dict=True)[0]
        eigenvectors[name] = sp.simplify(vector.subs(solution).subs(mixing_order, 1))
    return locals()


class DirectedMixingChecks(unittest.TestCase):
    """Verify the general transformation and the CP-odd specialization."""

    def test_general_flavor_entries(self):
        results = self.results
        p, q = results['p'], results['q']
        flavor_0a, flavor_bara = results['flavor_0a'], results['flavor_bara']
        flavor_a0, flavor_abar = results['flavor_a0'], results['flavor_abar']
        expected = sp.Matrix([[0, 0, flavor_0a/(2*p) + flavor_bara/(2*q)],
                              [0, 0, flavor_0a/(2*p) - flavor_bara/(2*q)],
                              [p*flavor_a0 + q*flavor_abar,
                               p*flavor_a0 - q*flavor_abar, 0]])
        for difference in results['directed_mixing'] - expected:
            self.assertEqual(sp.simplify(difference), 0)
        for difference in results['general_epsilon_mixing'] - expected.subs(results['epsilon_substitutions']):
            self.assertEqual(sp.simplify(difference), 0)

    def test_cp_odd_finite_epsilon_entries(self):
        results = self.results
        epsilon_K, normalization = results['epsilon_K'], results['normalization']
        flavor_0a, flavor_a0 = results['flavor_0a'], results['flavor_a0']
        expected = sp.Matrix([[0, 0, -normalization*epsilon_K*flavor_0a/(1-epsilon_K**2)],
                              [0, 0, normalization*flavor_0a/(1-epsilon_K**2)],
                              [2*epsilon_K*flavor_a0/normalization,
                               2*flavor_a0/normalization, 0]])
        for difference in results['cp_odd_pole_mixing'] - expected:
            self.assertEqual(sp.simplify(difference), 0)


def latex_output(expression):
    """Render flavor symbols with the neutral-kaon labels used in memory."""
    return (sp.latex(expression)
            .replace(r'H_{bara}', r'H_{\bar 0a}')
            .replace(r'H_{\bar{a}}', r'H_{a\bar 0}'))


def main():
    results = derive()
    for test_name in ('test_general_flavor_entries', 'test_cp_odd_finite_epsilon_entries'):
        check = DirectedMixingChecks(test_name)
        check.results = results
        check.debug()
    print('Constant widths; hbar = c = 1; pole energies and mixing entries in GeV.')
    for label, key in [('Mixing after flavor-to-pole similarity transformation', 'directed_mixing'),
                       ('General flavor mixing at finite epsilon_K', 'general_epsilon_mixing'),
                       ('CP-odd mixing at finite epsilon_K', 'cp_odd_pole_mixing'),
                       ('Effective Hamiltonian', 'hamiltonian'),
                       ('Three-pole characteristic polynomial', 'polynomial'),
                       ('Schur-complement equation', 'secular'),
                       ('Matching-CP two-level characteristic polynomial', 'two_polynomial'),
                       ('Exact two-level poles', 'roots'),
                       ('Squared complex splitting', 'splitting_squared'),
                       ('Two-level right eigenvector at a pole', 'two_level_vector'),
                       ('Squared splitting at equal masses and real symmetric mixing', 'equal_mass_splitting')]:
        print(f'{label}:\n{latex_output(results[key])}\n')
    for name in ('S', 'L', 'a'):
        print(f'{name} second-order pole shift:\n{latex_output(results["shifts"][name])}\n')
        print(f'{name} first-order right eigenvector:\n{latex_output(results["eigenvectors"][name])}\n')


if __name__ == '__main__':
    main()
