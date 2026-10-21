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
from scipy.integrate import quad

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

    def amplitude(s12, s23):
        """Evaluate the invariant amplitude for the requested ordered channel."""
        invariants = {(0, 1): s12, (1, 2): s23,
                      (0, 2): ma**2+sum(mass**2 for mass in masses)-s12-s23}
        result = -mass_trace/(6*fa*F_PI)
        for differentiated in range(3):
            pair = tuple(index for index in range(3) if index != differentiated)
            result += derivative_traces[differentiated]*(
                ma**2+masses[differentiated]**2-invariants[pair])/(12*fa*F_PI)
        return complex(result)
        for spectator, first, second, mass, width, flavour in exchanges:
            pair_invariant = invariants[(first, second)]
            first_invariant = invariants[tuple(sorted((first, spectator)))]
            second_invariant = invariants[tuple(sorted((second, spectator)))]
            first_squared, second_squared, spectator_squared = (
                masses[index]**2 for index in (first, second, spectator))
            pole_numerator = pair_invariant/4*(
                (pair_invariant+second_squared-first_squared)*(first_invariant-first_squared-spectator_squared)
                -(pair_invariant+first_squared-second_squared)*(second_invariant-second_squared-spectator_squared))
            matching = ((first_invariant-first_squared-spectator_squared)*(ma**2+second_squared-first_invariant)
                        -(second_invariant-second_squared-spectator_squared)*(ma**2+first_squared-second_invariant))/4
            result += flavour/(2*F_PI*fa*mass**2)*(
                pole_numerator/(mass**2-pair_invariant-1j*mass*width)+matching)
        return complex(result)
    return amplitude


def three_pseudoscalar_width(ma: float, fa: float, couplings: ALPcouplings,
                              meson_1: str, meson_2: str, meson_3: str) -> float:
    """Return Gamma(a -> P1 P2 P3) in GeV, including vector interference.

    Identical particles contribute 1/prod(n!). Numerical integration uses
    relative tolerance 1e-5. Thresholds use observed pseudoscalar masses.
    """
    names = meson_1, meson_2, meson_3
    particles = _pseudoscalars()
    first_mass, second_mass, third_mass = (particles[name][0] for name in names)
    if ma <= first_mass+second_mass+third_mass:
        return 0.0
    amplitude = three_pseudoscalar_amplitude(ma, fa, couplings, *names)

    def integrate_second_invariant(s12):
        """Integrate over the second Dalitz invariant at fixed s12."""
        center = second_mass**2+third_mass**2+(
            (s12+second_mass**2-first_mass**2)*(ma**2-s12-third_mass**2))/(2*s12)
        product = ((s12-(first_mass+second_mass)**2)*(s12-(first_mass-second_mass)**2)
                   *((ma-third_mass)**2-s12)*((ma+third_mass)**2-s12))
        half_width = np.sqrt(max(0., product))/(2*s12)
        return quad(lambda s23: abs(amplitude(s12, s23))**2,
                    center-half_width, center+half_width, epsabs=0., epsrel=1e-5, limit=200)[0]

    integral = quad(integrate_second_invariant, (first_mass+second_mass)**2,
                    (ma-third_mass)**2, epsabs=0., epsrel=1e-5, limit=200)[0]
    symmetry = np.prod([factorial(number) for number in Counter(names).values()])
    return float(integral/(256*np.pi**3*ma**3*symmetry))
