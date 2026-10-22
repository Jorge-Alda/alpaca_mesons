"""Odd resonance contribution to ALP decays to a photon and two mesons.

Masses, fa, G_V and widths use GeV. The tree amplitude retains the original
seven VJP and four VVP operators, including direct photon insertion, vector
conversion and the charged-source c4 term. Local anomalous contacts, even
radiative terms, loops and rescattering are omitted. Internal vectors use
fitted masses and constant PDG widths; external mesons have physical masses
and vertices have the library's LO Alda:2026nrd masses and generators.
"""
from dataclasses import dataclass
import math

import numpy as np
from scipy.integrate import quad

from .pseudoscalar_diag import ALPcouplings, _no_mixing_couplings
from .three_pseudoscalar import F_PI, _pseudoscalars
from .vector_inputs import _read_vector_inputs
from .vector_two_pseudoscalar import (
    _anticommutator, _commutator, _meson_inputs, _real_components,
    _sum_flavour_coefficients,
)


def _canonicalize_terms(terms):
    """Set a channel-wide relative roundoff cancellation to exact zero."""
    total = sum((constant+linear for _, _, _, constant, linear in terms), 0j)
    scale = sum(abs(constant)+abs(linear)
                for _, _, _, constant, linear in terms)
    if abs(total) <= 64*np.finfo(float).eps*scale:
        return ()
    return tuple(terms)


def _photon_two_pseudoscalar_t_limits(ma, first_mass, second_mass, s):
    """Return the physical t boundaries in GeV^2 at fixed s (GeV^2)."""
    discriminant = (s-(first_mass+second_mass)**2)*(s-(first_mass-second_mass)**2)
    center = second_mass**2 + (ma**2-s)*(s+second_mass**2-first_mass**2)/(2*s)
    half_range = (ma**2-s)*math.sqrt(max(0., discriminant))/(2*s)
    return center-half_range, center+half_range


def _canonical_meson_name(name):
    """Return a supported meson name in the library's canonical spelling."""
    aliases = {'k+': 'K+', 'k-': 'K-', 'k0': 'K0', 'k0bar': 'K0bar',
               'kbar0': 'K0bar', 'anti-k0': 'K0bar', 'ks': 'KS', 'k_s': 'KS',
               'kl': 'KL', 'k_l': 'KL',
               "eta'": 'eta_prime', 'etaprime': 'eta_prime'}
    if isinstance(name, str) and name in _pseudoscalars():
        return name
    if isinstance(name, str):
        canonical = aliases.get(name.lower())
        if canonical is not None:
            return canonical
    raise ValueError(f'unsupported pseudoscalar: {name!r}')


def _prepare_photon_two_pseudoscalar_terms(ma: float, fa: float, couplings: ALPcouplings,
                                           meson_1: str, meson_2: str):
    """Prepare flavour coefficients for the odd scalar amplitude.

    Returns the external-state metadata and ungrouped pole terms used by the
    optimized scalar evaluator. The physics inputs and channel selection are
    unchanged from the former public amplitude constructor.
    """
    first_name, second_name = map(_canonical_meson_name, (meson_1, meson_2))
    state = (first_name, second_name)
    particles = _pseudoscalars()
    masses = tuple(particles[name][0] for name in state)
    charges = {'pi+': 1, 'K+': 1, 'pi-': -1, 'K-': -1}
    if ma <= sum(masses) or sum(charges.get(name, 0) for name in state):
        return ma, masses, state, ()
    axial, pseudoscalar, _ = _no_mixing_couplings(couplings, ma)
    fitted = _read_vector_inputs()
    g_v = fitted['even_parameters']['G_V_GeV']
    generators, tree_masses = _meson_inputs()
    first, second = generators[first_name], generators[second_name]
    chi, odd = fitted['quark_mass_matrix'], fitted['odd_coefficients']
    charge = np.diag([2/3, -1/3, -1/3])
    electric_charge = math.sqrt(4*math.pi*fitted['alpha_em'])
    photon_response = sum(
        fitted['photon_vector_residues'][name]/fitted['vector_masses'][name]**2
        * fitted['vector_generators'][name] for name in ('rho0', 'omega', 'phi'))
    vector_widths = {'rho0': .1474, 'omega': .00868, 'phi': .004249,
                     'rhoplus': .1491, 'Kstarplus': .0514, 'Kstar0': .0473}
    vector_groups = [(fitted['vector_masses'][name], vector_widths[name],
                      (fitted['vector_generators'][name],))
                     for name in ('rho0', 'omega', 'phi')]
    for row, column, name in ((0, 1, 'rhoplus'), (0, 2, 'Kstarplus'), (1, 2, 'Kstar0')):
        mass = fitted['even_parameters']['m_rho_GeV' if (row, column) == (0, 1) else 'm_Kstar_GeV']
        vector_groups.append((mass, vector_widths[name], _real_components(row, column)))
    contributions = {}
    for mass, width, components in vector_groups:
        for internal in components:
            for invariant, source, mass_source, source_squared, emitted, trace, factor in (
                (0, 2*axial/fa, -2*_anticommutator(pseudoscalar, chi)/fa,
                 ma**2, None, np.trace(internal@_commutator(first, second)),
                 -g_v/(math.sqrt(2)*F_PI**2)),
                (1, -second/F_PI, _anticommutator(second, chi)/F_PI,
                 tree_masses[second_name]**2, second,
                 np.trace(internal@_commutator(first, axial)), math.sqrt(2)*g_v/(F_PI*fa)),
                (2, -first/F_PI, _anticommutator(first, chi)/F_PI,
                 tree_masses[first_name]**2, first,
                 np.trace(internal@_commutator(second, axial)), -math.sqrt(2)*g_v/(F_PI*fa))):
                derivative_response = np.trace(_anticommutator(internal, photon_response)@source)/2
                mass_response = np.trace(_anticommutator(internal, photon_response)@mass_source)/2
                derivative_charge = np.trace(_anticommutator(internal, charge)@source)/2
                mass_charge = np.trace(_anticommutator(internal, charge)@mass_source)/2
                constant = math.sqrt(2)*electric_charge*4*(
                    (odd['d1']-odd['d3'])*source_squared*derivative_response-4*odd['d2']*mass_response)
                linear = math.sqrt(2)*electric_charge*4*odd['d3']*derivative_response
                constant -= 8*electric_charge/fitted['operator_mass']*(
                    derivative_charge*source_squared*(odd['c1']+odd['c2']-odd['c5'])/2
                    -2*odd['c3']*mass_charge)
                linear -= 8*electric_charge/fitted['operator_mass']*derivative_charge*(
                    (-odd['c1']+odd['c2']+odd['c5'])/2-odd['c6'])
                if emitted is not None:
                    constant -= 8*electric_charge*odd['c4']/(F_PI*fitted['operator_mass'])*np.trace(
                        internal@_commutator(_commutator(emitted, charge), chi))
                constants, linears = contributions.setdefault((invariant, mass, width), ([], []))
                constants.append(factor*trace*constant)
                linears.append(factor*trace*linear)
    terms = []
    for key, (constants, linears) in contributions.items():
        constant, linear = map(_sum_flavour_coefficients, (constants, linears))
        if constant != 0j or linear != 0j:
            terms.append((*key, constant, linear))
    terms = _canonicalize_terms(terms)
    return ma, masses, state, terms


@dataclass(frozen=True)
class _PhotonTwoPseudoscalarAmplitudeOptimized:
    """Prepared scalar arithmetic for the same odd amplitude, in GeV^-3.

    Each invariant group contains (complex pole, constant, linear), with
    pole = m_V^2 - i m_V Gamma_V. No propagator approximation is made.
    """
    ma: float
    masses: tuple
    final_state: tuple
    terms: tuple
    invariant_terms: tuple
    total_squared: float
    threshold_squared: float
    mass_difference_squared: float
    second_mass_squared: float
    daughter_squared_difference: float

    def t_limits(self, s):
        """Physical t boundaries in GeV^2 using prepared mass combinations."""
        photon_energy_factor = (self.ma*self.ma-s)/(2*s)
        center = self.second_mass_squared + photon_energy_factor*(
            s+self.daughter_squared_difference)
        half_range = photon_energy_factor*math.sqrt(max(
            0., (s-self.threshold_squared)*(s-self.mass_difference_squared)))
        return center-half_range, center+half_range

    def at_fixed_s(self, s):
        """Return the scalar evaluator with all s-only work already done."""
        s_terms, t_terms, u_terms = self.invariant_terms
        scalar_s = 0j
        for pole, constant, linear in s_terms:
            scalar_s += (constant+linear*s)/(s-pole)
        u_offset = self.total_squared-s

        def scalar(t):
            """Coherent scalar coefficient at the prepared s, in GeV^-3."""
            value = scalar_s
            for pole, constant, linear in t_terms:
                value += (constant+linear*t)/(t-pole)
            u = u_offset-t
            for pole, constant, linear in u_terms:
                value += (constant+linear*u)/(u-pole)
            return value

        return scalar

    def __call__(self, s, t):
        """Evaluate the coherent scalar coefficient in GeV^-3."""
        if not self.terms:
            return 0j
        return self.at_fixed_s(s)(t)

    def squared(self, s, t):
        """Return the photon polarization sum at physical Dalitz coordinates."""
        if not self.terms:
            return 0.
        lower, upper = self.t_limits(s)
        return max(0., s*(t-lower)*(upper-t)/4)*abs(self(s, t))**2


def photon_two_pseudoscalar_amplitude(
        ma: float, fa: float, couplings: ALPcouplings,
        meson_1: str, meson_2: str):
    """Return the optimized callable odd scalar coefficient in GeV^-3.

    Inputs, units, retained contributions and conventions are those of
    the retained odd $a\to\gamma PP$ model. The amplitude is
    epsilon(eps_gamma,p1,p2,k) times this coefficient, with metric (+---) and
    epsilon^0123=+1. The callable takes s=(p1+p2)^2 and t=(k+p2)^2 in GeV^2.
    `at_fixed_s(s)` reuses the s-channel contribution in inner integrations.
    Threshold and charge-forbidden pairs return an empty amplitude.
    """
    ma, masses, state, terms = _prepare_photon_two_pseudoscalar_terms(
        ma, fa, couplings, meson_1, meson_2)
    first_mass, second_mass = masses
    groups = tuple(tuple((complex(mass*mass, -mass*width),
                          complex(constant), complex(linear))
                         for index, mass, width, constant, linear in terms
                         if index == invariant) for invariant in range(3))
    return _PhotonTwoPseudoscalarAmplitudeOptimized(
        ma, masses, state, terms,
        groups, ma*ma+first_mass*first_mass+second_mass*second_mass,
        (first_mass+second_mass)**2, (first_mass-second_mass)**2,
        second_mass*second_mass, second_mass*second_mass-first_mass*first_mass)


def photon_two_pseudoscalar_width(
        ma: float, fa: float, couplings: ALPcouplings,
        meson_1: str, meson_2: str, *, order: int = 16,
        epsrel: float = 1e-5) -> float:
    """Return the odd photon PP width in GeV with hybrid quadrature.

    Equal-mass pairs use vectorized Gauss-Legendre quadrature in both Dalitz
    variables, with a sine map in s to smooth the threshold. Unequal-mass
    pairs retain adaptive outer integration because crossed-pole entry can
    converge slowly with a global rule. Both variables are split at poles
    and at two and eight m_V Gamma_V on either side. Orders n and 2n are
    compared, doubling up to 8n before adaptive fallback. Inner convergence
    uses epsrel/4. The amplitude, finite widths, physical Dalitz domain and
    symmetry factor are unchanged. Only numerical quadrature is approximated.

    ``order`` is the initial positive integer quadrature order; ``epsrel``
    controls the outer relative tolerance and inner convergence requirement.
    """
    amplitude = photon_two_pseudoscalar_amplitude(
        ma, fa, couplings, meson_1, meson_2)
    if not amplitude.terms:
        return 0.
    from scipy.special import roots_legendre

    rules = [roots_legendre(order*factor) for factor in (1, 2, 4, 8)]
    polynomial_nodes, polynomial_weights = roots_legendre(2)
    polynomial_integral = float(polynomial_weights@(1-polynomial_nodes**2))
    s_terms, t_terms, u_terms = amplitude.invariant_terms

    def integrate_t(s):
        """Integrate t with explicit finite-width resonance subdivisions."""
        lower, upper = amplitude.t_limits(s)
        if upper <= lower:
            return 0.
        scalar = amplitude.at_fixed_s(s)
        if not t_terms and not u_terms:
            # The polarization factor is quadratic in t; two nodes are exact.
            half_range = (upper-lower)/2
            return s/4*half_range**3*polynomial_integral*abs(scalar((lower+upper)/2))**2
        boundaries = {lower, upper}
        for reflected, terms in ((False, t_terms), (True, u_terms)):
            for pole, _, _ in terms:
                center = amplitude.total_squared-s-pole.real if reflected else pole.real
                for offset in (-8., -2., 0., 2., 8.):
                    point = center+offset*abs(pole.imag)
                    if lower < point < upper:
                        boundaries.add(point)
        boundaries = np.array(sorted(boundaries))
        centers = (boundaries[1:]+boundaries[:-1])/2
        half_ranges = (boundaries[1:]-boundaries[:-1])/2
        previous = None
        for nodes, weights in rules:
            coordinates = centers[:, None]+half_ranges[:, None]*nodes
            integrand = (s/4*(coordinates-lower)*(upper-coordinates)
                         * np.abs(scalar(coordinates))**2)
            integral = float(np.sum(half_ranges*(integrand@weights)))
            if previous is not None and abs(integral-previous) <= epsrel/4*abs(integral):
                return integral
            previous = integral

        def squared(t):
            """Scalar evaluation for the unresolved inner integral."""
            return max(0., s/4*(t-lower)*(upper-t))*abs(scalar(t))**2

        return quad(squared, lower, upper, points=boundaries[1:-1],
                    epsabs=0., epsrel=epsrel/4, limit=150)[0]

    lower, upper = amplitude.threshold_squared, ma*ma
    points = sorted({pole.real for pole, _, _ in s_terms if lower < pole.real < upper})
    s_boundaries = {lower, upper, *points}
    for pole, _, _ in s_terms:
        for offset in (-8., -2., 2., 8.):
            point = pole.real+offset*abs(pole.imag)
            if lower < point < upper:
                s_boundaries.add(point)
    s_boundaries = np.array(sorted(s_boundaries))
    previous = None
    outer_rules = rules if amplitude.masses[0] == amplitude.masses[1] else ()
    for nodes, weights in outer_rules:
        # Evaluate all outer nodes together, including the crossed-pole splits.
        angles = (nodes+1)*math.pi/4
        fractions = np.sin(angles)**2
        s_intervals = np.diff(s_boundaries)
        s_values = (s_boundaries[:-1, None]+s_intervals[:, None]*fractions).ravel()
        s_weights = (s_intervals[:, None]*weights*math.pi/4*np.sin(2*angles)).ravel()
        photon_factor = (ma*ma-s_values)/(2*s_values)
        t_centers = amplitude.second_mass_squared+photon_factor*(
            s_values+amplitude.daughter_squared_difference)
        t_half_ranges = photon_factor*np.sqrt(np.maximum(
            0., (s_values-amplitude.threshold_squared)*
            (s_values-amplitude.mass_difference_squared)))
        t_lower, t_upper = t_centers-t_half_ranges, t_centers+t_half_ranges
        t_boundaries = [t_lower, t_upper]
        for reflected, terms in ((False, t_terms), (True, u_terms)):
            for pole, _, _ in terms:
                center = amplitude.total_squared-s_values-pole.real if reflected else pole.real
                for offset in (-8., -2., 0., 2., 8.):
                    t_boundaries.append(np.clip(
                        np.broadcast_to(center+offset*abs(pole.imag), s_values.shape),
                        t_lower, t_upper))
        t_boundaries = np.sort(np.array(t_boundaries).T, axis=1)
        t_ranges = np.diff(t_boundaries, axis=1)/2
        t_values = ((t_boundaries[:, 1:]+t_boundaries[:, :-1])[:, :, None]/2
                    +t_ranges[:, :, None]*nodes)
        s_grid = s_values[:, None, None]
        squared = (s_grid/4*(t_values-t_lower[:, None, None])*
                   (t_upper[:, None, None]-t_values)*np.abs(amplitude(np.broadcast_to(s_grid, t_values.shape), t_values))**2)
        inner_integrals = np.sum(t_ranges*(squared@weights), axis=1)
        integral = float(s_weights@inner_integrals)
        if previous is not None and abs(integral-previous) <= epsrel/4*abs(integral):
            break
        previous = integral
    else:
        integral = quad(integrate_t, lower, upper, points=points,
                        epsabs=0., epsrel=epsrel, limit=150)[0]
    symmetry = 2 if amplitude.final_state[0] == amplitude.final_state[1] else 1
    return integral/(256*math.pi**3*ma**3*symmetry)
