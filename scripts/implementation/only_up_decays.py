#!/usr/bin/env python3
"""Only-up phase model: LO finite-mixing decay amplitudes.

The electromagnetic coupling is fixed by alpha_em = 1/137. Outputs are
LO tree amplitudes in the conventions of the imported calculators. The
radial scalar, resonance exchange, loops and NLO terms are omitted.
See only_up_inputs.py for the quark masses and source derivatives.
Exception: three_pseudoscalar and three_pseudoscalar_width include the
radial-mode (sigma) exchange when radial_mass_squared is given; by default
they omit it, as before.

"""
import math
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import quad

_IMPLEMENTATION_DIRECTORY = str(Path(__file__).resolve().parent)
if _IMPLEMENTATION_DIRECTORY not in sys.path:
    sys.path.insert(0, _IMPLEMENTATION_DIRECTORY)
from only_up_inputs import only_up_inputs
from calculate_alp_meson_masses_finite import calculate_masses
from alp_photon_pp_wzw_finite_mixing import external_legs, wzw_coefficient


def spectrum(inputs):
    return calculate_masses(inputs['f_a'], inputs['m_a0'], inputs['hat_S'],
                            inputs['hat_P'], inputs['hat_Theta'], B_0=inputs['B_0'],
                            f_pi=inputs['f_pi'], m_0_squared=inputs['m_0_squared'])


def two_photon_coefficient(inputs):
    """Return C with L = alpha C a_phys F Ftilde/(4 pi f_a)."""
    legs, _ = external_legs('pi+', 'pi-', inputs['m_a0'], inputs['f_a'],
                            inputs['hat_S'], inputs['hat_P'], inputs['hat_Theta'],
                            f_pi=inputs['f_pi'], B_0=inputs['B_0'],
                            m_0_squared=inputs['m_0_squared'])
    charges_squared = np.diag([4/9, 1/9, 1/9])
    return float(np.real(3*inputs['f_a']/inputs['f_pi']*np.trace(charges_squared @ legs[0])))


def photon_two_pseudoscalar(inputs, mesons):
    """Return the WZW amplitude coefficient in GeV^-3."""
    return wzw_coefficient(*mesons, inputs['m_a0'], inputs['f_a'], inputs['hat_S'],
                           inputs['hat_P'], inputs['hat_Theta'], f_pi=inputs['f_pi'],
                           B_0=inputs['B_0'], m_0_squared=inputs['m_0_squared'], e=math.sqrt(4*math.pi/137))


def two_photon_width(inputs):
    """Return the LO WZW two-photon width in GeV.

    The coefficient convention is
    L = alpha C a_phys F Ftilde/(4 pi f_a), with alpha = e^2/(4 pi).
    Vector-meson contributions interfere with the WZW amplitude.
    """
    mass = spectrum(inputs)['a']
    coefficient = two_photon_coefficient(inputs)
    alpha = 1/137
    return alpha**2*coefficient**2*mass**3/(64*math.pi**3*inputs['f_a']**2)


def photon_two_pseudoscalar_width(inputs, mesons, *, epsrel=1e-8):
    """Return the LO WZW gamma PP width in GeV using mixed pole masses.

    Integrate the photon polarization sum analytically over the second
    Dalitz invariant and use adaptive one-dimensional quadrature for s.
    The constant WZW amplitude interferes with vector-meson amplitudes;
    the sum of their separate widths is not the complete width.
    Neutral kaons are degenerate in this diagonal only-up model.
    """
    if len(mesons) != 2:
        raise ValueError('provide two meson labels')
    masses = spectrum(inputs)
    kaon_names = ('K0', 'K0bar', 'K_S', 'K_L')
    def meson_mass(name):
        return masses['K0_1'] if name in kaon_names else masses[name]
    first_mass, second_mass = (meson_mass(name) for name in mesons)
    parent_mass = masses['a']
    if parent_mass <= first_mass + second_mass:
        return 0.
    coefficient = photon_two_pseudoscalar(inputs, mesons)
    if coefficient == 0:
        return 0.
    lower = (first_mass+second_mass)**2
    upper = parent_mass**2
    def integrand(s):
        kallen = (s-(first_mass+second_mass)**2)*(s-(first_mass-second_mass)**2)
        return (upper-s)**3*max(0., kallen)**1.5/s**2
    integral = quad(integrand, lower, upper, epsabs=0., epsrel=epsrel)[0]
    symmetry = 2 if mesons[0] == mesons[1] else 1
    return abs(coefficient)**2*integral/(6144*math.pi**3*parent_mass**3*symmetry)


class RadialExchange:
    """Radial-mode exchange in a_phys -> P_1 P_2 P_3 of the only-up model, polar coordinates.

    Value at (s12, s23) in GeV^2: sum over the channels (r|uv) of
    g_ar(s_uv) g_uv(s_uv)/(m_sigma^2 - s_uv), with the off-shell vertex of
    memory/alp/three-pseudoscalar/alp-three-pseudoscalar-radial-exchange.md,
        f_a g_xy(q^2) = (B0/2)<S_10 {L_x, L_y}> - (B0 f_pi/f_a)<P_11 (xi_x L_y + xi_y L_x)>
                        - (2 B0 f_pi^2/f_a^2) xi_x xi_y <S_12> + xi_x xi_y (q^2 - m_x^2 - m_y^2),
    on the original polar columns (L_x, xi_x) of the legs.  The columns are those of the
    amplitude without sigma (amplitude.external), so the relative sign is fixed.  In the
    only-up model the radial field multiplies the up-quark entry of the source,
    M(a, sigma) = M(a, 0) + (sigma/f_a) e_1 M(a, 0), e_1 = diag(1, 0, 0), hence
    S_10 = e_1 S_00, P_11 = e_1 P_01 and S_12 = e_1 S_02.  The meson and ALP exchanges of this
    model vanish (memory/alp/three-pseudoscalar/murayama-cubic-vertex-parity-rule.md).
    """

    def __init__(self, inputs, amplitude, radial_mass_squared):
        projector = np.diag([1., 0., 0.])
        hat_S, hat_P = inputs['hat_S'], inputs['hat_P']
        self.radial_mass_squared = radial_mass_squared
        self.scalar_radial = projector @ hat_S[0]
        self.pseudoscalar_radial_angular = projector @ hat_P[1]
        self.scalar_radial_angular_squared = projector @ hat_S[2]/2
        self.inputs, self.legs = inputs, amplitude.external
        self.masses_squared = [amplitude.masses['a']**2] + [amplitude.masses[name]**2 for name in amplitude.mesons]
        self.channels = []
        for spectator, first, second in ((1, 2, 3), (2, 1, 3), (3, 1, 2)):
            if self.coupling(0, spectator, 1.) != 0 and self.coupling(first, second, 1.) != 0:
                self.channels.append((spectator, first, second))

    def coupling(self, x, y, momentum_squared):
        """g_xy(q^2) in GeV for legs x, y (0 = parent, 1..3 = daughters)."""
        B_0, f_pi, f_a = self.inputs['B_0'], self.inputs['f_pi'], self.inputs['f_a']
        (first, first_alp), (second, second_alp) = self.legs[x], self.legs[y]
        value = (B_0/2*np.trace(self.scalar_radial @ (first @ second + second @ first))
                 - B_0*f_pi/f_a*np.trace(self.pseudoscalar_radial_angular @ (first_alp*second + second_alp*first))
                 - 2*B_0*f_pi**2/f_a**2*first_alp*second_alp*np.trace(self.scalar_radial_angular_squared)
                 + first_alp*second_alp*(momentum_squared - self.masses_squared[x] - self.masses_squared[y]))
        return value/f_a

    def __call__(self, s12, s23):
        s13 = sum(self.masses_squared) - s12 - s23
        invariants = {(1, 2): s12, (2, 3): s23, (1, 3): s13}
        value = 0.
        for spectator, first, second in self.channels:
            momentum_squared = invariants[(first, second)]
            value = value + (self.coupling(0, spectator, momentum_squared)*self.coupling(first, second, momentum_squared)
                             / (self.radial_mass_squared - momentum_squared))
        return value


def three_pseudoscalar(inputs, mesons, radial_mass_squared=None):
    """Prepare the amplitude as a function of s12 and s23 in GeV^2.

    K_S and K_L are outgoing states projected with the complex PDG 2024
    epsilon_K approximation in memory/kaon/neutral-kaon-amplitude-projection.md.
    Both neutral kaons use the common leading-order flavour-state mass.

    Without radial_mass_squared (default) the radial mode is omitted. With
    radial_mass_squared = m_sigma^2 in GeV^2 the sigma exchange (RadialExchange) is
    added; in this model m_sigma^2 = B0 m_u f_pi^2/f_a^2 + 2 lambda f_a^2 with the
    quartic coupling lambda > 0 of DiLuzio:2026poh.
    """
    from calculate_alp_three_pseudoscalar_finite import prepare_amplitude

    # Charge conjugation equates the two flavor orderings for a neutral,
    # C-even spectator in the only-up model, including radial exchange.
    if len(mesons) == 3 and 'K_S' in mesons and 'K_L' in mesons and any(
            name in ('pi0', 'eta', 'eta_prime') for name in mesons):
        def vanishing_amplitude(s12, s23):
            return 0j
        vanishing_amplitude.components = []
        return vanishing_amplitude

    if any(name in ('K_S', 'K_L') for name in mesons):
        from itertools import product

        # PDG 2024 magnitude and phenomenological phase; see the kaon mixing note.
        epsilon_k = 2.228e-3 * np.exp(1j * np.deg2rad(43.52))
        normalization = np.sqrt(2 * (1 + abs(epsilon_k)**2))
        p_conjugate = np.conjugate((1 + epsilon_k) / normalization)
        q_conjugate = np.conjugate((1 - epsilon_k) / normalization)
        choices = [((('K0_1', (p_conjugate + (q_conjugate if name == 'K_S' else -q_conjugate))/np.sqrt(2)),
                     ('K0_2', 1j*(-p_conjugate + (q_conjugate if name == 'K_S' else -q_conjugate))/np.sqrt(2)))
                    if name in ('K_S', 'K_L') else ((name, 1),)) for name in mesons]
        terms = []
        for components in product(*choices):
            flavour_labels, coefficients = zip(*components)
            terms.append((np.prod(coefficients), three_pseudoscalar(inputs, flavour_labels, radial_mass_squared)))
        def projected(s12, s23):
            return sum(coefficient * amplitude(s12, s23) for coefficient, amplitude in terms)

        projected.components = terms
        return projected

    amplitude = prepare_amplitude(inputs['f_a'], inputs['m_a0'], inputs['hat_S'], inputs['hat_P'],
                                  inputs['hat_Theta'], mesons, f_pi=inputs['f_pi'],
                                  B_0=inputs['B_0'], m_0_squared=inputs['m_0_squared'])
    if radial_mass_squared is None:
        return amplitude
    radial = RadialExchange(inputs, amplitude, radial_mass_squared)
    total = lambda s12, s23: amplitude(s12, s23) + radial(s12, s23)
    total.without_radial, total.radial = amplitude, radial
    return total


def three_pseudoscalar_width(inputs, mesons, radial_mass_squared=None, *, epsrel=1e-5):
    """Return the LO three-body width in GeV.

    Without radial_mass_squared (default) the radial mode is omitted, as in
    calculate_alp_three_pseudoscalar_finite.decay_width. With radial_mass_squared
    = m_sigma^2 in GeV^2 the sigma exchange is included coherently (see
    three_pseudoscalar). A sigma or meson pole inside the Dalitz region of a
    channel with nonzero couplings needs a width and raises ValueError.
    """
    from calculate_alp_three_pseudoscalar_finite import decay_width

    physical_kaons = any(name in ('K_S', 'K_L') for name in mesons)
    if radial_mass_squared is None and not physical_kaons:
        return decay_width(inputs['f_a'], inputs['m_a0'], inputs['hat_S'], inputs['hat_P'],
                           inputs['hat_Theta'], mesons, f_pi=inputs['f_pi'],
                           B_0=inputs['B_0'], m_0_squared=inputs['m_0_squared'])
    from collections import Counter
    from alpaca_chiral._three_body_quadrature import integrate_three_body

    charge = {'pi+': 1, 'pi-': -1, 'K+': 1, 'K-': -1}
    if sum(charge.get(name, 0) for name in mesons) != 0:
        return 0.
    masses = spectrum(inputs)
    parent = masses['a']
    daughters = [masses['K0_1' if name in ('K_S', 'K_L') else name] for name in mesons]
    if parent <= sum(daughters):
        return 0.
    total = three_pseudoscalar(inputs, mesons, radial_mass_squared)
    components = total.components if physical_kaons else [(1, total)]
    if not components:
        return 0.
    pole_groups = {1: [], 2: [], 3: []}
    for coefficient, component in components:
        amplitude = component.without_radial if radial_mass_squared is not None else component
        radial = component.radial if radial_mass_squared is not None else None
        for spectator, residues in amplitude.exchanges:
          others = [daughters[i - 1] for i in (1, 2, 3) if i != spectator]
          low, high = sum(others)**2, (parent - daughters[spectator - 1])**2
          poles = []
          for pole in sorted(set(amplitude.poles)):
              if sum(residues[amplitude.poles == pole]) != 0:
                  poles.append(complex(pole))
          if radial is not None and any(spectator == channel[0] for channel in radial.channels):
              poles.append(complex(radial_mass_squared))
          if any(low <= pole.real <= high for pole in poles):
              raise ValueError('accessible internal pole: stable tree-level width diverges')
          pole_groups[spectator].extend(poles)
    first, second, third = daughters

    def kallen(x, y, z):
        return (x - y - z)**2 - 4*y*z

    def limits(s):
        root = np.sqrt(np.maximum(0., kallen(s, first**2, second**2)*kallen(parent**2, s, third**2)))
        midpoint = second**2 + third**2 + (s + second**2 - first**2)*(parent**2 - s - third**2)/(2*s)
        return midpoint - root/(2*s), midpoint + root/(2*s)

    integral = integrate_three_body(lambda s: (lambda t: np.abs(total(s, t))**2), limits,
                                    (first + second)**2, (parent - third)**2,
                                    parent**2 + first**2 + second**2 + third**2,
                                    (pole_groups[3], pole_groups[1], pole_groups[2]), epsrel=epsrel)
    symmetry = np.prod([math.factorial(count) for count in Counter(mesons).values()])
    return integral/(256*math.pi**3*parent**3*symmetry)



def sigma_two_pseudoscalar_width(kappa_u, quartic_coupling, mesons):
    """Return the tree-level sigma -> P_i P_j width in GeV.

    The inputs are kappa_u, the dimensionless scalar coupling lambda, and
    two physical field labels: 'a', 'pi0', 'eta', 'eta_prime', 'pi+',
    'pi-', 'K+', 'K-', 'K0', 'K0bar', 'K_S', or 'K_L'. The calculation
    uses the Cartesian vertex and the repository's LO finite-mixing mass
    calculator and neutral pole vectors. DiLuzio:2026poh supplies the model,
    not a generic formula for the diagonalized masses.
    """
    from calculate_alp_three_pseudoscalar_finite import prepare_amplitude

    if len(mesons) != 2:
        raise ValueError('provide two pseudoscalar labels')
    inputs = only_up_inputs(kappa_u)
    allowed_neutral = {'a', 'pi0', 'eta', 'eta_prime'}
    charged = {'pi+', 'pi-', 'K+', 'K-'}
    neutral_kaons = {'K0', 'K0bar', 'K_S', 'K_L'}
    if not set(mesons) <= allowed_neutral | charged | neutral_kaons:
        raise ValueError('unknown pseudoscalar label')
    # The up-quark projector has no neutral-kaon vertex. All remaining
    # pairs outside the neutral sector and charged conjugate pairs violate
    # electric charge or strangeness conservation.
    if not (set(mesons) <= allowed_neutral or set(mesons) in ({'pi+', 'pi-'}, {'K+', 'K-'})):
        return 0.
    labels = [name if name != 'a' else 'pi0' for name in mesons]
    amplitude = prepare_amplitude(inputs['f_a'], inputs['m_a0'], inputs['hat_S'],
                                  inputs['hat_P'], inputs['hat_Theta'], labels + ['pi0'],
                                  f_pi=inputs['f_pi'], B_0=inputs['B_0'],
                                  m_0_squared=inputs['m_0_squared'])
    legs = [amplitude.external[0] if name == 'a' else amplitude.external[index+1]
            for index, name in enumerate(mesons)]
    masses = amplitude.masses
    first_mass, second_mass = (masses[name] for name in mesons)
    f_a, B_0, f_pi = inputs['f_a'], inputs['B_0'], inputs['f_pi']
    up_mass = inputs['hat_S'][0][0, 0]
    sigma_squared = B_0*up_mass*f_pi**2/f_a**2 + 2*quartic_coupling*f_a**2
    sigma_mass = math.sqrt(sigma_squared)
    if sigma_mass <= first_mass + second_mass:
        return 0.
    (first, first_angular), (second, second_angular) = legs
    scalar_source = np.diag([up_mass, 0., 0.])
    coupling = (B_0/(2*f_a)*np.trace(scalar_source @ (first @ second + second @ first))
                + 2*quartic_coupling*f_a*first_angular*second_angular)
    phase_space = math.sqrt((sigma_squared-(first_mass+second_mass)**2)
                            *(sigma_squared-(first_mass-second_mass)**2))
    return float(abs(coupling)**2*phase_space/(16*math.pi*sigma_mass**3*(1+(mesons[0] == mesons[1]))))
