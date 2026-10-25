"""Match the strict ALP--vector--three-pion pole phase.

Run .venv/bin/python scripts/derive_alp_three_pion_pole_phase.py
or append --test. Natural units, GeV, metric (+---), epsilon^0123=+1,
S=1+iT. Only signed on-shell residues and ordinary neutral-boson crossing
are used; an experimental decay-rate normalization has no absolute phase.
"""

from functools import lru_cache
import sys
import unittest

import sympy as sp

import derive_alp_doubly_virtual_tff as tree
import derive_physical_alp_vgamma as radiative


def show(label, expression):
    print(f"{label}:\n{sp.latex(expression)}\n")


@lru_cache(maxsize=1)
def derive():
    """Derive numerator matching, boson crossing, and propagator phases."""
    stored = radiative.derive()
    direct_trace, mass_trace = stored['derivative_trace'], stored['mass_trace']
    direct, pair = tree.ordered_numerators(direct_trace, mass_trace)
    direct = direct.subs({tree.tensor.m: stored['alp_mass'],
                         tree.tensor.MV: stored['operator_mass']})
    pair = pair.subs(tree.tensor.m, stored['alp_mass'])
    electric_charge, mass = stored['electric_charge'], stored['external_mass']
    matching = {}
    numerators = {}
    for name, generator in stored['generators'].items():
        direct_source = direct.subs({
            direct_trace: sp.trace(generator * stored['charge'] * stored['derivative_source']),
            mass_trace: sp.trace(generator * stored['charge'] * stored['mass_source']),
            tree.tensor.s1: mass**2, tree.tensor.s2: 0})
        pair_source = pair.subs({
            direct_trace: sp.trace(generator * stored['response_diagonal'] * stored['derivative_source']),
            mass_trace: sp.trace(generator * stored['response_diagonal'] * stored['mass_source']),
            tree.tensor.s1: mass**2, tree.tensor.s2: 0})
        numerator = sp.expand(-2 * electric_charge * direct_source
                              + sp.sqrt(2) * electric_charge * pair_source)
        numerators[name] = numerator
        matching[name] = sp.simplify(numerator / mass - stored['amplitudes'][name])

    metric = sp.diag(1, -1, -1, -1)
    momentum = sp.Matrix(sp.symbols('k_0:4', real=True))
    incoming_photon = sp.Matrix(sp.symbols('q_0:4', real=True))
    vector_polarization = sp.Matrix(sp.symbols('n_0:4'))
    photon_polarization = sp.Matrix(sp.symbols('epsilon_0:4'))
    current = sp.Matrix(sp.symbols('J_0:4'))

    def lower_epsilon_contraction(first, second, third, fourth):
        # Lower epsilon_0123=-1 in the declared metric. SymPy's built-in
        # LeviCivita is the canonical all-upper symbol in this expression.
        return sp.expand(sum(-sp.LeviCivita(a, b, c, d)
                             * first[a] * second[b] * third[c] * fourth[d]
                             for a in range(4) for b in range(4)
                             for c in range(4) for d in range(4)))

    scattering_tensor = lower_epsilon_contraction(
        vector_polarization, photon_polarization, momentum, incoming_photon)
    decay_tensor = lower_epsilon_contraction(
        vector_polarization, photon_polarization, momentum, -incoming_photon)
    crossing_ratio = sp.cancel(decay_tensor / scattering_tensor)
    zeta = sp.Symbol('zeta_i')
    conventional_phase = sp.solve(sp.Eq(zeta, crossing_ratio), zeta)[0]

    lowered_momentum, lowered_photon = metric * momentum, metric * incoming_photon
    production_tensor = sp.Matrix(4, 4, lambda a, mu: sum(
        sp.LeviCivita(a, mu, r, s) * lowered_momentum[r] * lowered_photon[s]
        for r in range(4) for s in range(4)))
    spin_sum = -metric + momentum * momentum.T / mass**2
    contraction = production_tensor.T * metric * spin_sum * metric * current
    transverse_contraction = -production_tensor.T * metric * current
    projector_residual = (contraction - transverse_contraction).applyfunc(sp.simplify)

    virtuality = sp.Symbol('s', real=True)
    production_numerator = sp.Function('N_i')(mass**2, sp.Symbol('y', real=True))
    decay_function = sp.Symbol('F_i_on')
    scattering_coefficient = conventional_phase * production_numerator / mass
    # M_prod=C_i epsilon(n*,epsilon,k,q), M_dec=i F_i_on n.J.
    # iM_prod * iPi/(s-M_i^2+i0) * iM_dec; divide by i for M.
    amplitude_coefficient = sp.simplify(
        (sp.I * scattering_coefficient) * (sp.I / (virtuality - mass**2))
        * (sp.I * sp.I * decay_function) / sp.I)
    recorded_coefficient = (sp.I * zeta * production_numerator * decay_function
                            / (mass * (mass**2 - virtuality)))
    composition_difference = sp.simplify(
        amplitude_coefficient - recorded_coefficient.subs(zeta, conventional_phase))

    phase = sp.Symbol('varphi_i', real=True)
    normalization_phase = sp.exp(sp.I * phase)
    rate_phase_residual = sp.simplify(normalization_phase * sp.conjugate(normalization_phase) - 1)
    return locals()


class PolePhaseTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = derive()

    def test_signed_real_photon_production(self):
        self.assertEqual(self.result['matching'], {'omega': 0, 'phi': 0})

    def test_crossing_and_ward(self):
        result = self.result
        self.assertEqual(sp.expand(result['decay_tensor'] + result['scattering_tensor']), 0)
        self.assertEqual(result['conventional_phase'], -1)
        tensor, metric = result['production_tensor'], result['metric']
        self.assertEqual((tensor * metric * result['incoming_photon']).applyfunc(sp.expand), sp.zeros(4, 1))
        self.assertEqual((tensor.T * metric * result['momentum']).applyfunc(sp.expand), sp.zeros(4, 1))

    def test_spin_sum_and_propagator_composition(self):
        self.assertEqual(self.result['projector_residual'], sp.zeros(4, 1))
        self.assertEqual(self.result['composition_difference'], 0)

    def test_decay_rate_does_not_set_global_phase(self):
        self.assertEqual(self.result['rate_phase_residual'], 0)


def main():
    result = derive()
    print('GeV, metric (+---), epsilon^0123=+1; S=1+iT; strict pole only.')
    print('Photon crossing uses the same neutral-photon LSZ phase on both amplitudes.')
    for name in ('omega', 'phi'):
        show(f'N_{name}(M_ext^2,0)/M_ext minus stored signed G', result['matching'][name])
    show('Outgoing-photon decay epsilon tensor', result['decay_tensor'])
    show('Incoming-photon scattering epsilon tensor', result['scattering_tensor'])
    show('Crossing ratio derived from the two tensors', result['crossing_ratio'])
    show('Conventional zeta_i for q incoming and M_dec=i F_i_on n.J', result['conventional_phase'])
    show('Spin-one projector contraction residual', result['projector_residual'])
    show('Composed coefficient multiplying E^T g Pi g J', result['amplitude_coefficient'])
    show('Difference from recorded strict pole with matched zeta', result['composition_difference'])
    show('Decay-rate phase invariance residual', result['rate_phase_residual'])
    print('A physical overall decay phase relative to the signed source remains unmeasured by decay rates.')


if __name__ == '__main__':
    if '--test' in sys.argv:
        unittest.main(argv=[sys.argv[0]])
    else:
        main()
