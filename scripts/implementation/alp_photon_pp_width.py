"""Selected a -> gamma P P odd resonance widths, in GeV.

Implements memory/alp/three-meson/alp-photon-pp-amplitude.md and its width
record, with no-mixing c_q^A,c_q^P supplied directly. The default inputs are
the preliminary + G_V fit. This is only the retained odd resonance piece;
local anomalous contacts, even radiative terms, loops and rescattering are
outside the calculation. Internal vectors have constant-width denominators
r^2-M_R^2+i M_R Gamma_R. Vertex masses are LO; external masses are physical.

Example::

    from scripts.implementation.alp_photon_pp_width import alp_photon_pp_width
    width = alp_photon_pp_width(1., [1., 0., 0.], [0., 0., 0.],
                                ('pi+', 'pi-'), 1000.)
"""
from dataclasses import dataclass
import math

import numpy as np
from scipy.integrate import quad

from .alp_vpp_width import (
    _anticommutator, _commutator, _couplings, _meson_inputs,
    _positive_number, _real_components,
)
from .alp_vv_width import read_preliminary_fit_parameters
from .three_meson_width import (
    EPSILON_K, _MESONS, _canonical_name, _kaon_mixing_coefficients,
)
from .vector_meson_exchange import _aggregate_flavor_coefficients


@dataclass(frozen=True)
class PhotonPPAmplitude:
    """One prepared ordered meson pair and its scalar amplitude in GeV^-3.

    Terms contain (invariant index, pole mass, pole width, constant, linear).
    Indices 0,1,2 denote s=(p1+p2)^2,t=(k+p2)^2,u=(k+p1)^2.
    Coefficients of identical propagators are summed before evaluation.
    """
    alp_mass: float
    masses: tuple
    final_state: tuple
    terms: tuple

    def __call__(self, s, t):
        """Evaluate the coherent scalar; an empty model returns exact zero."""
        if not self.terms:
            return 0j
        invariants = (s, t, self.alp_mass**2 + sum(m*m for m in self.masses)-s-t)
        return sum(((constant+linear*invariants[index]) /
                    (invariants[index]-mass**2+1j*mass*width)
                    for index, mass, width, constant, linear in self.terms), 0j)

    def squared(self, s, t):
        """Photon-polarization sum, using the stable physical Dalitz form."""
        if not self.terms:
            return 0.
        lower, upper = photon_pp_t_limits(self.alp_mass, *self.masses, s)
        # det Gram = s (t-t_minus)(t_plus-t)/4 in the physical region.
        determinant = max(0., s*(t-lower)*(upper-t)/4)
        return determinant*abs(self(s, t))**2


def photon_pp_t_limits(alp_mass, first_mass, second_mass, s):
    """Return the physical t boundaries in GeV^2 at fixed s (GeV^2)."""
    discriminant = (s-(first_mass+second_mass)**2)*(s-(first_mass-second_mass)**2)
    center = second_mass**2 + (alp_mass**2-s)*(s+second_mass**2-first_mass**2)/(2*s)
    half_range = (alp_mass**2-s)*math.sqrt(max(0., discriminant))/(2*s)
    return center-half_range, center+half_range


def _flavor_projection(name, epsilon_k):
    """Outgoing bra coefficients, retaining the complete complex epsilon_K."""
    if name not in ('KS', 'KL'):
        return ((name, 1.+0j),)
    p_mixing, q_mixing = _kaon_mixing_coefficients(epsilon_k)
    return (('K0', p_mixing.conjugate()),
            ('K0bar', (1 if name == 'KS' else -1)*q_mixing.conjugate()))


def build_photon_pp_amplitude(alp_mass, c_q_A, c_q_P, final_state, f_a, *,
                               fit_parameters=None, g_v=None, epsilon_k=EPSILON_K):
    """Prepare one pair (P1,P2), with dimensionless no-mixing source inputs.

    Source matrices may be diagonal three-vectors or Hermitian 3x3 matrices
    commuting with electric charge. f_a, optional G_V and masses use GeV.
    Default G_V is the conditional value in the preliminary fit. Physical
    K_S/K_L states are projected coherently before pole coefficients are
    aggregated. Only the requested pair is constructed.
    """
    alp_mass = _positive_number(alp_mass, 'alp_mass')
    f_a = _positive_number(f_a, 'f_a')
    axial, pseudoscalar = _couplings(c_q_A), _couplings(c_q_P)
    first_name, second_name = map(_canonical_name, final_state)
    state = (first_name, second_name)
    masses = tuple(_MESONS[name].mass for name in state)
    if alp_mass <= sum(masses) or sum(_MESONS[name].charge for name in state):
        return PhotonPPAmplitude(alp_mass, masses, state, ())
    fitted = fit_parameters if fit_parameters is not None else read_preliminary_fit_parameters()
    g_v = _positive_number(fitted.even_parameters['G_V_GeV'] if g_v is None else g_v, 'g_v')
    generators, tree_masses = _meson_inputs(epsilon_k)
    pion_constant = 0.093
    chi, odd = fitted.quark_mass_matrix, fitted.odd_coefficients
    charge = np.diag([2/3, -1/3, -1/3])
    electric_charge = math.sqrt(4*math.pi*fitted.alpha_em)
    photon_response = sum((fitted.photon_vector_residues[name]/fitted.vector_masses[name]**2
                           *fitted.vector_generators[name] for name in ('rho0', 'omega', 'phi')))
    vector_groups = [(fitted.vector_masses[name], fitted.shared_total_widths[name],
                      (fitted.vector_generators[name],)) for name in ('rho0', 'omega', 'phi')]
    for row, column, width_name in ((0, 1, 'rhoplus'), (0, 2, 'Kstarplus'), (1, 2, 'Kstar0')):
        mass = fitted.even_parameters['m_rho_GeV' if (row, column) == (0, 1) else 'm_Kstar_GeV']
        vector_groups.append((mass, fitted.shared_total_widths[width_name], _real_components(row, column)))
    contributions = {}
    for first_flavor, first_weight in _flavor_projection(first_name, epsilon_k):
        for second_flavor, second_weight in _flavor_projection(second_name, epsilon_k):
            first, second = generators[first_flavor].conj().T, generators[second_flavor].conj().T
            projection = first_weight*second_weight
            for mass, width, components in vector_groups:
                for internal in components:
                    for index, source, mass_source, source_squared, emitted, trace, factor in (
                        (0, 2*axial/f_a, -2*_anticommutator(pseudoscalar, chi)/f_a,
                         alp_mass**2, None, np.trace(internal@_commutator(first, second)),
                         -g_v/(math.sqrt(2)*pion_constant**2)),
                        (1, -second/pion_constant, _anticommutator(second, chi)/pion_constant,
                         tree_masses[second_flavor]**2, second,
                         np.trace(internal@_commutator(first, axial)), math.sqrt(2)*g_v/(pion_constant*f_a)),
                        (2, -first/pion_constant, _anticommutator(first, chi)/pion_constant,
                         tree_masses[first_flavor]**2, first,
                         np.trace(internal@_commutator(second, axial)), -math.sqrt(2)*g_v/(pion_constant*f_a))):
                        derivative_response = np.trace(_anticommutator(internal, photon_response)@source)/2
                        mass_response = np.trace(_anticommutator(internal, photon_response)@mass_source)/2
                        derivative_charge = np.trace(_anticommutator(internal, charge)@source)/2
                        mass_charge = np.trace(_anticommutator(internal, charge)@mass_source)/2
                        constant = math.sqrt(2)*electric_charge*4*(
                            (odd['d1']-odd['d3'])*source_squared*derivative_response-4*odd['d2']*mass_response)
                        linear = math.sqrt(2)*electric_charge*4*odd['d3']*derivative_response
                        constant -= 8*electric_charge/fitted.operator_mass*(
                            derivative_charge*source_squared*(odd['c1']+odd['c2']-odd['c5'])/2
                            -2*odd['c3']*mass_charge)
                        linear -= 8*electric_charge/fitted.operator_mass*derivative_charge*(
                            (-odd['c1']+odd['c2']+odd['c5'])/2-odd['c6'])
                        if emitted is not None:
                            constant -= 8*electric_charge*odd['c4']/(pion_constant*fitted.operator_mass)*np.trace(
                                internal@_commutator(_commutator(emitted, charge), chi))
                        constants, linears = contributions.setdefault((index, mass, width), ([], []))
                        constants.append(projection*factor*trace*constant)
                        linears.append(projection*factor*trace*linear)
    terms = []
    for key, (constants, linears) in contributions.items():
        constant, linear = map(_aggregate_flavor_coefficients, (constants, linears))
        if constant != 0j or linear != 0j:
            terms.append((*key, constant, linear))
    return PhotonPPAmplitude(alp_mass, masses, state, tuple(terms))


def alp_photon_pp_width(alp_mass, c_q_A, c_q_P, final_state, f_a, *,
                        fit_parameters=None, g_v=None, epsilon_k=EPSILON_K,
                        epsabs=0., epsrel=1e-5, limit=150):
    """Return the selected odd resonance partial width in GeV.

    Forbidden or exactly cancelled channels return 0 before quadrature.
    Nested scipy quad integrates s,t, subdividing at internal vector poles.
    epsabs applies to each unnormalized integral; epsrel is dimensionless.
    SciPy convergence warnings propagate to the caller.
    """
    amplitude = build_photon_pp_amplitude(alp_mass, c_q_A, c_q_P, final_state, f_a,
        fit_parameters=fit_parameters, g_v=g_v, epsilon_k=epsilon_k)
    if not amplitude.terms:
        return 0.
    alp_mass, (first_mass, second_mass) = amplitude.alp_mass, amplitude.masses
    total_squared = alp_mass**2+first_mass**2+second_mass**2

    def integrate_t(s):
        lower, upper = photon_pp_t_limits(alp_mass, first_mass, second_mass, s)
        points = set()
        for index, mass, _, _, _ in amplitude.terms:
            pole = mass**2 if index == 1 else total_squared-s-mass**2
            if index in (1, 2) and lower < pole < upper:
                points.add(pole)
        return quad(lambda t: amplitude.squared(s, t), lower, upper,
                    points=sorted(points), epsabs=epsabs, epsrel=epsrel, limit=limit)[0]

    lower, upper = (first_mass+second_mass)**2, alp_mass**2
    points = sorted({mass**2 for index, mass, _, _, _ in amplitude.terms
                     if index == 0 and lower < mass**2 < upper})
    integral = quad(integrate_t, lower, upper, points=points,
                    epsabs=epsabs, epsrel=epsrel, limit=limit)[0]
    symmetry = 2 if amplitude.final_state[0] == amplitude.final_state[1] else 1
    return integral/(256*math.pi**3*alp_mass**3*symmetry)
