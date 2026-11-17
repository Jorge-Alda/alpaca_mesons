#!/usr/bin/env python3
"""Only-up phase model: LO finite-mixing decay amplitudes.

The electromagnetic coupling is fixed by alpha_em = 1/137. Outputs are
LO tree amplitudes in the conventions of the imported calculators. The
radial scalar, resonance exchange, loops and NLO terms are omitted.
See only_up_inputs.py for the quark masses and source derivatives.

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
from calculate_alp_three_pseudoscalar_finite import prepare_amplitude, decay_width
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


def three_pseudoscalar(inputs, mesons):
    """Prepare the amplitude as a function of s12 and s23 in GeV^2."""
    return prepare_amplitude(inputs['f_a'], inputs['m_a0'], inputs['hat_S'], inputs['hat_P'],
                             inputs['hat_Theta'], mesons, f_pi=inputs['f_pi'],
                             B_0=inputs['B_0'], m_0_squared=inputs['m_0_squared'])


def three_pseudoscalar_width(inputs, mesons):
    """Return the LO three-body width in GeV."""
    return decay_width(inputs['f_a'], inputs['m_a0'], inputs['hat_S'], inputs['hat_P'],
                       inputs['hat_Theta'], mesons, f_pi=inputs['f_pi'],
                       B_0=inputs['B_0'], m_0_squared=inputs['m_0_squared'])

