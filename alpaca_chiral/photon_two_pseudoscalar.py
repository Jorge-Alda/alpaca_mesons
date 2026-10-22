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


@dataclass(frozen=True)
class _PhotonTwoPseudoscalarAmplitude:
    """One prepared ordered meson pair and its scalar amplitude in GeV^-3.

    Terms contain (invariant index, pole mass, pole width, constant, linear).
    Indices 0,1,2 denote s=(p1+p2)^2,t=(k+p2)^2,u=(k+p1)^2.
    Coefficients of identical propagators are summed before evaluation.
    """
    ma: float
    masses: tuple
    final_state: tuple
    terms: tuple

    def __call__(self, s, t):
        """Evaluate the coherent scalar; an empty model returns exact zero."""
        if not self.terms:
            return 0j
        invariants = (s, t, self.ma**2 + sum(m*m for m in self.masses)-s-t)
        return sum(((constant+linear*invariants[index]) /
                    (invariants[index]-mass**2+1j*mass*width)
                    for index, mass, width, constant, linear in self.terms), 0j)

    def squared(self, s, t):
        """Photon-polarization sum, using the stable physical Dalitz form."""
        if not self.terms:
            return 0.
        lower, upper = _photon_two_pseudoscalar_t_limits(self.ma, *self.masses, s)
        # det Gram = s (t-t_minus)(t_plus-t)/4 in the physical region.
        determinant = max(0., s*(t-lower)*(upper-t)/4)
        return determinant*abs(self(s, t))**2


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


def photon_two_pseudoscalar_amplitude(ma: float, fa: float, couplings: ALPcouplings,
                                       meson_1: str, meson_2: str):
    """Return the callable odd scalar coefficient in GeV^-3 for a -> gamma PP.

    The amplitude is epsilon(eps_gamma,p1,p2,k) times this coefficient,
    with metric (+---) and epsilon^0123=+1. The callable takes
    s=(p1+p2)^2 and t=(k+p2)^2 in GeV^2 and supplies the photon polarization
    sum through squared(s,t). Meson names follow three_pseudoscalar_amplitude.
    KS and KL use coherent outgoing states with the library's complex
    epsilon_K. ALPcouplings is rotated to the no-mixing basis with
    pseudoscalar_diag.py. Sources are assumed to commute with electric charge.
    Threshold and charge-forbidden pairs return an empty amplitude.
    """
    first_name, second_name = map(_canonical_meson_name, (meson_1, meson_2))
    state = (first_name, second_name)
    particles = _pseudoscalars()
    masses = tuple(particles[name][0] for name in state)
    charges = {'pi+': 1, 'K+': 1, 'pi-': -1, 'K-': -1}
    if ma <= sum(masses) or sum(charges.get(name, 0) for name in state):
        return _PhotonTwoPseudoscalarAmplitude(ma, masses, state, ())
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
    return _PhotonTwoPseudoscalarAmplitude(ma, masses, state, terms)


def photon_two_pseudoscalar_width(ma: float, fa: float, couplings: ALPcouplings,
                                   meson_1: str, meson_2: str) -> float:
    """Return the retained odd Gamma(a -> gamma P1 P2) in GeV.

    Nested Dalitz quadrature uses relative tolerance 1e-5 and subdivisions
    at internal vector poles. The identical-meson factor is included.
    Closed, charge-forbidden and exactly cancelled channels return zero;
    SciPy quadrature convergence warnings propagate to the caller.
    """
    amplitude = photon_two_pseudoscalar_amplitude(ma, fa, couplings, meson_1, meson_2)
    epsabs, epsrel, limit = 0., 1e-5, 150
    if not amplitude.terms:
        return 0.
    ma, (first_mass, second_mass) = amplitude.ma, amplitude.masses
    total_squared = ma**2+first_mass**2+second_mass**2

    def integrate_t(s):
        """Integrate the photon polarization sum over t at fixed s."""
        lower, upper = _photon_two_pseudoscalar_t_limits(ma, first_mass, second_mass, s)
        points = set()
        for index, mass, _, _, _ in amplitude.terms:
            pole = mass**2 if index == 1 else total_squared-s-mass**2
            if index in (1, 2) and lower < pole < upper:
                points.add(pole)
        return quad(lambda t: amplitude.squared(s, t), lower, upper,
                    points=sorted(points), epsabs=epsabs, epsrel=epsrel, limit=limit)[0]

    lower, upper = (first_mass+second_mass)**2, ma**2
    points = sorted({mass**2 for index, mass, _, _, _ in amplitude.terms
                     if index == 0 and lower < mass**2 < upper})
    integral = quad(integrate_t, lower, upper, points=points,
                    epsabs=epsabs, epsrel=epsrel, limit=limit)[0]
    symmetry = 2 if amplitude.final_state[0] == amplitude.final_state[1] else 1
    return integral/(256*math.pi**3*ma**3*symmetry)
