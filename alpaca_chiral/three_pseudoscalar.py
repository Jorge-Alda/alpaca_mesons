"""LO chiral plus leading G_V exchange ALP three-pseudoscalar widths.

Masses, fa and widths are in GeV. Neutral vector poles use the fitted
non-ideal omega-phi mixing and fitted masses. Charged vectors use the
flavour basis and PDG pole masses. The model uses G_V=f_pi/sqrt(2),
constant PDG widths and tensor-to-Proca matching.
No rescattering, loops or additional odd-parity contribution is included.
"""
from collections import Counter
from itertools import permutations
from math import factorial

import numpy as np

from .pseudoscalar_diag import ALPcouplings, _no_mixing_couplings, fit_results
from .vector_inputs import _read_vector_inputs

F_PI = 0.093


def _commutator(first, second):
    """Return the flavour commutator."""
    return first @ second - second @ first


def _pseudoscalars():
    """Return (mass, outgoing-state generator) for each pseudoscalar."""
    charged = {}
    for name, mass, row, column in (
        ('pi+', .13957039, 0, 1), ('pi-', .13957039, 1, 0),
        ('K+', .493677, 0, 2), ('K-', .493677, 2, 0),
        ('K0', .497611, 1, 2), ('K0bar', .497611, 2, 1)):
        generator = np.zeros((3, 3), dtype=complex)
        generator[column, row] = np.sqrt(2)
        charged[name] = mass, generator
    for name, mass, key in (('pi0', .1349768, 'lam_pi0'),
                            ('eta', .547862, 'lam_eta'),
                            ('eta_prime', .95778, 'lam_etap')):
        charged[name] = mass, np.asarray(fit_results[key], dtype=complex)
    epsilon_k = 2.228e-3 * np.exp(1j*np.deg2rad(43.52))
    normalization = np.sqrt(2*(1+abs(epsilon_k)**2))
    kaon_p, kaon_q = (1+epsilon_k)/normalization, (1-epsilon_k)/normalization
    charged['KS'] = .497611, (np.conjugate(kaon_p)*charged['K0'][1]+np.conjugate(kaon_q)*charged['K0bar'][1])
    charged['KL'] = .497611, (np.conjugate(kaon_p)*charged['K0'][1]-np.conjugate(kaon_q)*charged['K0bar'][1])
    return charged


def _vector_poles():
    """Return (mass, width, generator), with Tr(T T dagger)=1.

    Neutral masses and non-ideal generators come from the provisional fit;
    charged states retain the flavour basis and PDG masses. Widths are PDG.
    """
    poles = {}
    for name, mass, width, row, column in (
        ('rho+', .77511, .1491, 0, 1), ('rho-', .77511, .1491, 1, 0),
        ('K*+', .89167, .0514, 0, 2), ('K*-', .89167, .0514, 2, 0),
        ('K*0', .89555, .0473, 1, 2), ('K*0bar', .89555, .0473, 2, 1)):
        generator = np.zeros((3, 3), dtype=complex)
        generator[row, column] = 1
        poles[name] = mass, width, generator
    fitted = _read_vector_inputs()
    for name, width in (('rho0', .1474), ('omega', .00868), ('phi', .004249)):
        poles[name] = (fitted['vector_masses'][name], width,
                       fitted['vector_generators'][name])
    return poles


def three_pseudoscalar_amplitude(ma: float, fa: float, couplings: ALPcouplings,
                                  meson_1: str, meson_2: str, meson_3: str):
    """Return callable M(s12,s23), including G_V exchange and matching.

    The ordered invariants are s12=(p1+p2)^2 and s23=(p2+p3)^2.
    Supported particles: pi+, pi-, pi0, K+, K-, K0, K0bar, KS, KL,
    eta, eta_prime. Outgoing generators are Hermitian conjugated.
    All six labelled permutations retain identical-particle multiplicities.
    Masses, decay constants and widths use GeV. Scalar factors are prepared
    once, with finite-width poles and matching retained exactly.
    """
    particles = _pseudoscalars()
    masses, generators = zip(*(particles[name] for name in (meson_1, meson_2, meson_3)))
    axial, pseudoscalar, _ = _no_mixing_couplings(couplings, ma)
    quark_masses = np.asarray(fit_results['B0mq'], dtype=complex)
    mass_source = pseudoscalar @ quark_masses + quark_masses @ pseudoscalar
    derivative_traces = [0j, 0j, 0j]
    mass_trace = 0j
    for first, second, differentiated in permutations(range(3)):
        derivative_traces[differentiated] += np.trace(
            axial @ _commutator(generators[first],
                                _commutator(generators[second], generators[differentiated])))
        mass_trace += np.trace(generators[first] @ generators[second]
                              @ generators[differentiated] @ mass_source)
    exchange_coefficients = {}
    for spectator in range(3):
        first, second = (index for index in range(3) if index != spectator)
        for mass, width, generator in _vector_poles().values():
            flavour = np.trace(generator.conj().T @ _commutator(generators[first], generators[second]))
            flavour *= np.trace(generator @ _commutator(generators[spectator], axial))
            if flavour != 0:
                key = spectator, first, second, mass, width
                exchange_coefficients.setdefault(key, []).append(flavour)
    # Cancel flavour coefficients sharing the same pole and kinematics first.
    exchanges = []
    for key, coefficients in exchange_coefficients.items():
        flavour = sum(coefficients)
        coefficient_scale = sum(abs(coefficient) for coefficient in coefficients)
        if abs(flavour) > 64*np.finfo(float).eps*coefficient_scale:
            exchanges.append((*key, flavour))

    squared_masses = tuple(float(mass*mass) for mass in masses)
    parent_squared = ma*ma
    total_squared = parent_squared + sum(squared_masses)
    normalization = 12*fa*F_PI
    # Substitute s13=total_squared-s12-s23 in the LO chiral polynomial.
    chiral_constant = complex(-mass_trace/(6*fa*F_PI) + sum(
        derivative_traces[index]*(parent_squared+squared_masses[index])
        for index in range(3))/normalization
        - derivative_traces[1]*total_squared/normalization)
    chiral_first = complex((derivative_traces[1]-derivative_traces[2])/normalization)
    chiral_second = complex((derivative_traces[1]-derivative_traces[0])/normalization)
    invariant_indices = {(0, 1): 0, (1, 2): 1, (0, 2): 2}
    grouped_exchanges = {}
    for spectator, first, second, mass, width, flavour in exchanges:
        pole = complex(mass*mass, -mass*width)
        grouped_exchanges.setdefault((spectator, first, second), []).append(
            (pole, complex(flavour*pole/(8*F_PI*fa*mass*mass))))
    prepared_exchanges = tuple((
        invariant_indices[(first, second)],
        invariant_indices[tuple(sorted((first, spectator)))],
        invariant_indices[tuple(sorted((second, spectator)))],
        (squared_masses[second]-squared_masses[first])
        *(parent_squared-squared_masses[spectator]), tuple(poles))
        for (spectator, first, second), poles in grouped_exchanges.items())

    def evaluate(s12, s23):
        """Evaluate the same coherent amplitude with prepared scalar factors.

        On shell, the pole numerator is s_ij times the matching polynomial.
        Their exact sum is evaluated with the same complex pole and width.
        No zero-width or low-energy expansion is made.
        """
        invariants = s12, s23, total_squared-s12-s23
        result = chiral_constant + chiral_first*s12 + chiral_second*s23
        for pair, first, second, constant, poles in prepared_exchanges:
            pair_invariant = invariants[pair]
            numerator = pair_invariant*(invariants[first]-invariants[second])+constant
            residue = 0j
            for pole, coefficient in poles:
                residue += coefficient/(pole-pair_invariant)
            result += numerator*residue
        return result

    def amplitude(s12, s23):
        """Return the complex scalar amplitude at one Dalitz point."""
        return complex(evaluate(s12, s23))

    amplitude.vectorized = evaluate
    amplitude.invariant_poles = tuple(
        tuple(pole for pair, _, _, _, poles in prepared_exchanges
              if pair == invariant for pole, _ in poles)
        for invariant in range(3))
    return amplitude


def three_pseudoscalar_width(
        ma: float, fa: float, couplings: ALPcouplings,
        meson_1: str, meson_2: str, meson_3: str, *,
        order: int = 16, epsrel: float = 1e-5) -> float:
    """Return Gamma(a -> 3P) in GeV using converged Gauss-Legendre rules.

    Retains the full coherent amplitude, finite widths, physical phase
    space and identical-particle factorials. Inner intervals are split at poles
    and ±2, ±8 mass times width; orders double four times with adaptive
    fallback. The outer grid uses a sine-squared endpoint map and order
    doubling, with adaptive outer fallback. ``epsrel`` sets outer
    relative tolerance and four times the inner convergence tolerance.
    """
    from ._three_body_quadrature import integrate_three_body

    names = meson_1, meson_2, meson_3
    masses = tuple(_pseudoscalars()[name][0] for name in names)
    if ma <= sum(masses):
        return 0.
    amplitude = three_pseudoscalar_amplitude(ma, fa, couplings, *names)
    first, second, third = masses

    def limits(s):
        """Physical s23 boundaries at fixed s12, in GeV squared."""
        center = second**2+third**2+(s+second**2-first**2)*(ma**2-s-third**2)/(2*s)
        product = (s-(first+second)**2)*(s-(first-second)**2)
        product *= ((ma-third)**2-s)*((ma+third)**2-s)
        half_range = np.sqrt(np.maximum(0., product))/(2*s)
        return center-half_range, center+half_range

    integral = integrate_three_body(
        lambda s: lambda t: np.abs(amplitude.vectorized(s, t))**2,
        limits, (first+second)**2, (ma-third)**2,
        ma**2+sum(mass**2 for mass in masses), amplitude.invariant_poles,
        order=order, epsrel=epsrel)
    symmetry = np.prod([factorial(number) for number in Counter(names).values()])
    return float(integral/(256*np.pi**3*ma**3*symmetry))
