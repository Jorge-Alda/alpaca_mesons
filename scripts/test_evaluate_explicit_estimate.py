"""Unit tests for evaluate_explicit_estimate.py.

No loop is computed: the tests check tree-level vertex formulas, the completeness identity,
pole data and the explicit diagram list.

Run: PYTHONPATH=scripts .venv/bin/python -m unittest test_evaluate_explicit_estimate
"""
import itertools
import math
import unittest

import numpy as np

import evaluate_explicit_estimate as explicit

ZERO = np.zeros((3, 3), complex)


def gell_mann_generators():
    """Nine Hermitian U(3) generators with <lambda_i lambda_j> = 2 delta_ij."""
    generators = [np.sqrt(2 / 3) * np.eye(3, dtype=complex)]
    for b in range(3):
        for c in range(b + 1, 3):
            symmetric, antisymmetric = np.zeros((3, 3), complex), np.zeros((3, 3), complex)
            symmetric[b, c] = symmetric[c, b] = 1
            antisymmetric[b, c], antisymmetric[c, b] = -1j, 1j
            generators += [symmetric, antisymmetric]
    generators += [np.diag([1, -1, 0]).astype(complex), np.diag([1, 1, -2]).astype(complex) / np.sqrt(3)]
    return generators


class SetupMixin:
    @classmethod
    def setUpClass(cls):
        cls.estimate = explicit.ExplicitEstimate(0.3)


class VertexFormulaTest(SetupMixin, unittest.TestCase):
    def test_four_meson_kinetic_closed_form(self):
        """(1/48) sum_pi (-k_pi2.k_pi4) <[L_pi1, L_pi2][L_pi3, L_pi4]> equals the general T/R formula."""
        generator = np.random.default_rng(3)
        matrices = []
        for _ in range(4):
            raw = generator.normal(size=(3, 3)) + 1j * generator.normal(size=(3, 3))
            matrices.append(raw + raw.conj().T)
        names = 'abcd'
        products = generator.normal(size=(4, 4))
        products = products + products.T

        def dot(first, second):
            return products[names.index(first), names.index(second)]

        legs = [explicit.Leg(matrices[i], 0j, 0j, names[i]) for i in range(4)]
        general = self.estimate.meson_kinetic_vertex(legs, dot)
        closed = 0j
        for order in itertools.permutations(range(4)):
            first = matrices[order[0]] @ matrices[order[1]] - matrices[order[1]] @ matrices[order[0]]
            second = matrices[order[2]] @ matrices[order[3]] - matrices[order[3]] @ matrices[order[2]]
            closed += -products[order[1], order[3]] * np.trace(first @ second)
        self.assertAlmostEqual(abs(general - closed / 48), 0., places=12)

    def test_up_only_potential_closed_form(self):
        """Up-only potential: B0 m_u sum eps^{r+s} prod rho prod xi (-1)^{(s+m)/2}/m! sum_pi (Lambda_pi)_11
        over r <= 1, s + m even, plus the (0, 0) term with diag(0, m_d, m_s)."""
        e = self.estimate
        generator = np.random.default_rng(5)
        legs = []
        for index in range(5):
            raw = generator.normal(size=(3, 3)) + 1j * generator.normal(size=(3, 3))
            legs.append(explicit.Leg(raw + raw.conj().T, complex(generator.normal()),
                                     complex(generator.normal()), str(index)))
        lighter = np.diag([0., e.down_mass, e.strange_mass])
        closed = 0j
        for assignment in itertools.product('RAP', repeat=len(legs)):
            r, s = assignment.count('R'), assignment.count('A')
            mesons = [leg.meson for leg, role in zip(legs, assignment) if role == 'P']
            m = len(mesons)
            if r > 1 or (s + m) % 2:
                continue
            factor = e.epsilon ** (r + s) / math.factorial(m)
            for leg, role in zip(legs, assignment):
                factor *= leg.radial if role == 'R' else leg.alp if role == 'A' else 1
            ordered = e.ordered_sum(mesons)
            closed += e.b_zero * e.up_mass * factor * (-1) ** ((s + m) // 2) * ordered[0, 0]
            if r == 0 and s == 0:
                closed += e.b_zero * (-1) ** (m // 2) / math.factorial(m) * np.trace(lighter @ ordered)
        self.assertLess(abs(e.potential_vertex(legs) - closed), 1e-9 * abs(closed))

    def test_sigma_cubed_vertex(self):
        """V(sigma, sigma, sigma) = -6 lambda f_a (f_pi = 1: f_a = 1/eps), multilinear in rho."""
        e = self.estimate
        sigma = explicit.Leg(ZERO, 0j, 1. + 0j, explicit.LOOP)
        expected = -6 * e.quartic_coupling / e.epsilon
        self.assertAlmostEqual(e.vertex([sigma, sigma, sigma], e.dot_rule(1.)).real, expected, places=14)
        value = e.vertex([e.legs['beta'], sigma, sigma], e.dot_rule(e.radial_mass_squared))
        self.assertLess(abs(value - expected * e.legs['beta'].radial), 1e-12 * abs(value))

    def test_tree_has_no_sigma_cubed_term(self):
        """The sigma^3 vertex needs three radial legs, so the tree vertices are unchanged."""
        e = self.estimate
        sigma = explicit.Leg(ZERO, 0j, 1. + 0j, 'sigma')
        dot = e.dot_rule(0.)
        for first, second in (('A', '0'), ('+', '-')):
            with self.subTest(pair=(first, second)):
                self.assertEqual(e.radial_cubic_vertex([sigma, e.legs[first], e.legs[second]]), 0j)
                self.assertNotEqual(e.vertex([sigma, e.legs[first], e.legs[second]], dot), 0j)

    def test_closed_alp_line_slope_commutator_form(self):
        """d/dp^2 (1/2) V(l_e(p), l_e'(-p), l_A, l_A) = (1/24) <[L_A, L_e][L_A, L_e']>; zero for the
        diagonal a and pi0 pole vectors of the up-only model, so the class-7 residues 7.1 (a) and
        7.4 (pi0) vanish."""
        e = self.estimate
        alp = e.legs['A'].meson
        for first, second in (('A', 'A'), ('0', '0'), ('+', '-')):
            first_matrix, second_matrix = e.legs[first].meson, e.legs[second].meson
            closed = np.trace((alp @ first_matrix - first_matrix @ alp) @ (alp @ second_matrix - second_matrix @ alp)) / 24
            with self.subTest(legs=(first, second)):
                self.assertLess(abs(e.closed_alp_line_slope(first, second) - closed),
                                1e-9 * max(abs(closed), 1e-6))
        for name in ('A', '0'):
            self.assertLess(abs(e.closed_alp_line_slope(name, name)), 1e-12)

    def test_charged_leg_residue_closed_form(self):
        """pi+ and pi- slopes d/dp^2 [(1/2) V(l_+(p), l_-(-p), l_A, l_A)] from the vertex equal
        -(eps^2/(3N)) [(M^-1)_11 - (M^-1)_22]^2, M^-1 = M^{-1}_{m_a^2}(C), N = N(m_a^2)."""
        for kappa_u in (0.3, 1.0):
            e = explicit.ExplicitEstimate(kappa_u)
            inverse = e.pole_operator_inverse(e.alp_mass_squared)
            closed = -(e.epsilon ** 2 / (3 * e.normalization(e.alp_mass_squared))) * (inverse[0, 0] - inverse[1, 1]) ** 2
            for first, second in (('+', '-'), ('-', '+')):
                with self.subTest(kappa_u=kappa_u, legs=(first, second)):
                    self.assertLess(abs(e.closed_alp_line_slope(first, second) - closed) / abs(closed), 1e-9)

    def test_class_seven_residue_includes_charged_legs(self):
        """Class-7 residues 7.1-7.4 = (1/2) sum_e slope_e M0 over e = a, pi+, pi-, pi0; the a and pi0
        slopes vanish, so the sum equals (1/2)(2 x pi+ slope) M0 (decision of 2026-11-25)."""
        e = self.estimate
        tree = e.tree()[0]
        self.assertEqual(explicit.RESIDUE_DIAGRAM_LABELS['7'], ('7.1', '7.2', '7.3', '7.4'))
        lines, value = e.external_residues(tree)['7']
        self.assertEqual(lines, ('A',))
        expected = e.charged_leg_alp_slope() * tree
        self.assertLess(abs(value - expected), 1e-12 * abs(expected))
        self.assertGreater(abs(value), 0.)

    def test_charged_radial_vertex_vanishes(self):
        """V(pi+, pi-, sigma, sigma) = 0 (sources linear in sigma), so class 8 has no charged-leg part."""
        e = self.estimate
        sigma = explicit.Leg(ZERO, 0j, 1. + 0j, explicit.LOOP)
        legs = [explicit.Leg(e.legs['+'].meson, 0j, 0j, 'p'), explicit.Leg(e.legs['-'].meson, 0j, 0j, 'q'), sigma, sigma]
        for momentum_squared in (0., 1.):
            dot = e.dot_rule(e.radial_mass_squared, {('p', 'q'): -momentum_squared, ('q', 'p'): -momentum_squared})
            self.assertEqual(e.vertex(legs, dot), 0j)

    def test_quadratic_terms_give_unit_kinetic_matrix(self):
        """Two-leg vertex: d V/d p^2 is the unit kinetic matrix on the generators and the ALP."""
        legs = [explicit.Leg(g, 0j, 0j, 'p') for g in gell_mann_generators()] + [explicit.Leg(ZERO, 1. + 0j, 0j, 'p')]
        slopes = np.zeros((10, 10))
        for i, first in enumerate(legs):
            for j, second in enumerate(legs):
                values = []
                for momentum_squared in (0., 1.):
                    dot = self.estimate.dot_rule(0., {('p', 'q'): -momentum_squared, ('q', 'p'): -momentum_squared})
                    values.append(self.estimate.vertex(
                        [first, explicit.Leg(second.meson, second.alp, 0j, 'q')], dot).real)
                slopes[i, j] = values[1] - values[0]
        np.testing.assert_allclose(slopes, np.eye(10), atol=1e-12)

    def test_radial_coupling_closed_form(self):
        """V(sigma, A, pi0) = -eps [q^2 - m_a^2 - m_0^2 + m_a^2 m_0^2/(B0 m_u eps^2)]/sqrt(N_a N_0) at q^2 = s_{+-}."""
        e = self.estimate
        sigma = explicit.Leg(ZERO, 0j, 1. + 0j, 'sigma')
        value = e.vertex([sigma, e.legs['A'], e.legs['0']], e.dot_rule(0.))
        alp, pion = e.alp_mass_squared, e.neutral_pion_mass_squared
        closed = (-e.epsilon / math.sqrt(e.normalization(alp) * e.normalization(pion))
                  * (e.pair_invariant - alp - pion + alp * pion / (e.b_zero * e.up_mass * e.epsilon ** 2)))
        self.assertLess(abs(value - closed) / abs(closed), 1e-8)
        charged = e.vertex([sigma, e.legs['+'], e.legs['-']], e.dot_rule(0.))
        self.assertAlmostEqual(charged.real, -e.b_zero * e.up_mass * e.epsilon, places=10)

    def test_vanishing_vertices(self):
        """Vertices behind the vanishing cycles: closed S line (rule R2), cubic vertex of one meson
        and two A lines, the charge-violating radial branches."""
        e = self.estimate
        sigma = explicit.Leg(ZERO, 0j, 1. + 0j, explicit.LOOP)
        alp_slot = e.slot_legs()[10]
        legs = e.legs
        dot = e.dot_rule(e.radial_mass_squared)
        checks = {
            'closed S line on four legs (rule R2)': [legs['A'], legs['+'], legs['-'], legs['0'], sigma, sigma],
            'closed S line with beta (rule R2)': [legs['beta'], legs['A'], legs['0'], sigma, sigma],
            'closed S line with gamma (rule R2)': [legs['gamma'], legs['+'], legs['-'], sigma, sigma],
            'cubic vertex, one meson and two A lines': [explicit.Leg(np.diag([1., 0, 0]).astype(complex), 0j, 0j, 'x'), alp_slot, alp_slot],
            'sigma A pi+': [explicit.Leg(ZERO, 0j, 1. + 0j, 'sigma'), legs['A'], legs['+']],
            'sigma pi0 pi-': [explicit.Leg(ZERO, 0j, 1. + 0j, 'sigma'), legs['0'], legs['-']],
            'four radial': [legs['beta'], legs['gamma'], sigma, sigma],
        }
        for name, vertex_legs in checks.items():
            with self.subTest(name):
                self.assertLess(abs(e.vertex(vertex_legs, dot)), 1e-10)


class CompletenessAndPoleTest(SetupMixin, unittest.TestCase):
    def test_u3_completeness(self):
        """sum_i F(lambda_i) G(lambda_i) = 2 sum_bc F(E_bc) G(E_cb) for vertex functions F, G."""
        e = self.estimate
        dot = e.dot_rule(e.alp_mass_squared)
        alp_slot = e.slot_legs()[10]

        def first(matrix):
            return e.vertex([e.legs['A'], e.legs['+'], alp_slot, explicit.Leg(matrix, 0j, 0j, explicit.LOOP)], dot)

        def second(matrix):
            return e.vertex([e.legs['-'], e.legs['0'], alp_slot, explicit.Leg(matrix, 0j, 0j, explicit.LOOP)], dot)

        generator_sum = sum(first(g) * second(g) for g in gell_mann_generators())
        elementary_sum = 2 * sum(first(explicit.unit_matrix(b, c)) * second(explicit.unit_matrix(c, b))
                                 for b in range(3) for c in range(3))
        self.assertLess(abs(generator_sum - elementary_sum), 1e-9 * abs(generator_sum))

    def test_normalization_is_derivative_of_pole_function(self):
        e = self.estimate
        for mass_squared in (e.alp_mass_squared, e.neutral_pion_mass_squared):
            step = 1e-6 * mass_squared
            derivative = (e.pole_function(mass_squared + step) - e.pole_function(mass_squared - step)) / (2 * step)
            self.assertLess(abs(derivative / e.normalization(mass_squared) - 1), 1e-6)
            self.assertLess(abs(e.pole_function(mass_squared)), 1e-9 * mass_squared)


class DiagramListTest(SetupMixin, unittest.TestCase):
    def test_list_equals_labelled_cycle_sum(self):
        """The explicit list with its symmetry factors equals sum over ordered partitions and line
        types with weight 1/(2n), by class (all line types, including the vanishing ones); the sigma
        self-energy 6.7 is compared separately from the rest of class 6."""
        e = self.estimate
        sigma_self_energy = 'sigma self-energy 6.7'

        def class_name(lines, vertices, branches):
            alp, radial = lines.count('A'), lines.count('S')
            if alp + radial == 0 or alp + radial >= 3:
                return None
            if branches == ('beta', 'gamma'):
                return sigma_self_energy if lines == ('S', 'S') else None
            if alp == 1 and radial == 0:
                return '1' if vertices == 1 else '2'
            if alp == 2:
                return '3'
            if radial == 1 and alp == 0:
                return 'closed S line (rule R2)' if vertices == 1 else '4'
            return '6' if radial == 2 else '5'

        labelled = {}
        for branches in (('A', '+', '-', '0'), ('beta', 'A', '0'), ('gamma', '+', '-'), ('beta', 'gamma')):
            for vertices in range(1, len(branches) + 1):
                for assignment in itertools.product(range(vertices), repeat=len(branches)):
                    blocks = [tuple(b for b, a in zip(branches, assignment) if a == k) for k in range(vertices)]
                    if not all(blocks):
                        continue
                    for lines in itertools.product('AMS', repeat=vertices):
                        name = class_name(lines, vertices, branches)
                        if name is None:
                            continue
                        value = e.cycle(tuple(blocks), lines, e.loop_scale_squared(lines)) / (2 * vertices)
                        labelled[name] = labelled.get(name, 0j) + value
        listed = {}
        for label, weight, blocks, lines in explicit.DIAGRAMS:
            name = sigma_self_energy if label == '6.7' else explicit.diagram_class(label)
            listed[name] = listed.get(name, 0j) + weight * e.cycle(blocks, lines, e.loop_scale_squared(lines))
        for name in set(labelled) | set(listed):
            with self.subTest(name):
                reference = max(abs(labelled.get(name, 0j)), 1e-30)
                self.assertLess(abs(labelled.get(name, 0j) - listed.get(name, 0j)) / reference, 1e-10)


if __name__ == '__main__':
    unittest.main()
