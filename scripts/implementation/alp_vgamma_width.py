"""No-mixing ALP decays to omega gamma and phi gamma.

Masses, f_a and widths are in GeV; amplitudes are in GeV^-1. The real
diagonal c_q^A and c_q^P inputs are in u,d,s order in Alda:2026nrd's
no-mixing basis. The calculation uses the original eleven odd operators,
epsilon^0123=+1, tree level, leading electromagnetic order and narrow
external vectors. It includes direct VJP and photon-converted VVP terms.
It adds no ALP--meson pole or independent local anomalous a V gamma term.

Run from the repository root, for example::

    .venv/bin/python -m scripts.implementation.alp_vgamma_width \
        --axial 1 0 0 --pseudoscalar 0 0 0 --ma 1.5 --fa 1000

Only fit loading and the photon--vector response use other project modules.
The default fit is the conditional negative-strange single-trace solution.
"""

import argparse
import json
import math

import numpy as np

from .alp_vv_width import DEFAULT_FIT_LABEL, read_preliminary_fit_parameters
from .vector_sector import VectorSector


def _positive_number(value, name):
    """Read a finite positive mass or decay constant."""
    if isinstance(value, (bool, complex)) or not np.isscalar(value):
        raise ValueError(name+' must be a positive finite real number')
    value = float(value)
    if not math.isfinite(value) or value <= 0:
        raise ValueError(name+' must be a positive finite real number')
    return value


def _quark_couplings(values, name):
    """Read three finite real diagonal quark couplings in u,d,s order."""
    couplings = np.asarray(values)
    if (couplings.shape != (3,) or np.iscomplexobj(couplings)
            or not np.all(np.isfinite(couplings))):
        raise ValueError(name+' must contain three finite real entries (u,d,s)')
    return couplings.astype(float)


def alp_vector_photon_amplitudes(axial_quark_couplings,
                                 pseudoscalar_quark_couplings,
                                 alp_mass, alp_decay_constant, *,
                                 fit_parameters=None):
    """Return omega gamma and phi gamma scalar amplitudes in GeV^-1.

    The tensor amplitude is G epsilon_{mu nu rho sigma} eps_V^mu
    eps_gamma^nu q_V^rho q_gamma^sigma. Couplings are dimensionless arrays
    of three real diagonal entries. All other inputs come from
    ``read_preliminary_fit_parameters`` unless its result is supplied.
    The fit reader retains fixed central inputs and its selected mass and
    source hypothesis; it currently accepts the two single-trace choices.
    Amplitudes are also returned below threshold as analytic expressions.
    """
    alp_mass = _positive_number(alp_mass, 'alp_mass')
    alp_decay_constant = _positive_number(alp_decay_constant, 'alp_decay_constant')
    fitted = (fit_parameters if fit_parameters is not None
              else read_preliminary_fit_parameters())
    axial = np.diag(_quark_couplings(axial_quark_couplings, 'axial_quark_couplings'))
    pseudoscalar = np.diag(_quark_couplings(
        pseudoscalar_quark_couplings, 'pseudoscalar_quark_couplings'))
    mass_anticommutator = (pseudoscalar @ fitted.quark_mass_matrix
                          + fitted.quark_mass_matrix @ pseudoscalar)
    charges = np.diag([2/3, -1/3, -1/3])
    vector_names = ('rho0', 'omega', 'phi')
    vectors = VectorSector(
        tuple(fitted.vector_masses[name] for name in vector_names),
        tuple(fitted.vector_generators[name] for name in vector_names),
        tuple(fitted.photon_vector_residues[name] for name in vector_names))
    # W_gamma(0) = sum_i R_i T_i / M_i^2, with distinct physical poles.
    photon_response = vectors.response(0.).real
    electric_charge = math.sqrt(4*math.pi*fitted.alpha_em)
    odd = fitted.odd_coefficients
    amplitudes = {}
    for name in ('omega', 'phi'):
        vector_mass = fitted.vector_masses[name]
        generator = fitted.vector_generators[name]
        direct_axial_trace = np.trace(generator @ charges @ axial)
        direct_mass_trace = np.trace(generator @ charges @ mass_anticommutator)
        converted_axial_trace = np.trace(generator @ photon_response @ axial)
        converted_mass_trace = np.trace(
            generator @ photon_response @ mass_anticommutator)
        # c4,c7,d4 have zero neutral real-photon contractions.
        direct = -8/(fitted.operator_mass*alp_decay_constant)*(
            direct_axial_trace*(
                odd['c1']*(alp_mass**2-vector_mass**2)
                + odd['c2']*(alp_mass**2+vector_mass**2)
                - odd['c5']*(alp_mass**2-vector_mass**2)
                - 2*odd['c6']*vector_mass**2)
            + 4*odd['c3']*direct_mass_trace)
        # The VVP expression already includes both vector assignments.
        converted = 8*math.sqrt(2)/alp_decay_constant*(
            converted_axial_trace*((odd['d1']-odd['d3'])*alp_mass**2
                                   + odd['d3']*vector_mass**2)
            + 4*odd['d2']*converted_mass_trace)
        amplitudes[name+'gamma'] = float(
            electric_charge/vector_mass*(direct+converted))
    return amplitudes


def alp_vector_photon_widths(axial_quark_couplings,
                             pseudoscalar_quark_couplings,
                             alp_mass, alp_decay_constant, *,
                             fit_parameters=None):
    """Return scalar amplitudes in GeV^-1 and partial widths in GeV.

    Gamma = |G|^2 (m_a^2-M_ext^2)^3 / (32 pi m_a^3), above threshold;
    the width is zero at and below threshold. There is no identical-particle
    factor or initial spin average. Daughter widths and fit uncertainties
    are not integrated or propagated.
    """
    alp_mass = _positive_number(alp_mass, 'alp_mass')
    fitted = (fit_parameters if fit_parameters is not None
              else read_preliminary_fit_parameters())
    amplitudes = alp_vector_photon_amplitudes(
        axial_quark_couplings, pseudoscalar_quark_couplings,
        alp_mass, alp_decay_constant, fit_parameters=fitted)
    widths = {}
    for name in ('omega', 'phi'):
        channel = name+'gamma'
        vector_mass = fitted.vector_masses[name]
        widths[channel] = (abs(amplitudes[channel])**2
                           * (alp_mass**2-vector_mass**2)**3
                           / (32*math.pi*alp_mass**3)
                           if alp_mass > vector_mass else 0.)
    return {'amplitudes_GeV_inv': amplitudes, 'widths_GeV': widths}


def main():
    """Read no-mixing inputs and print both amplitudes and widths as JSON."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--axial', nargs=3, type=float, required=True,
                        metavar=('CU', 'CD', 'CS'))
    parser.add_argument('--pseudoscalar', nargs=3, type=float, required=True,
                        metavar=('CU', 'CD', 'CS'))
    parser.add_argument('--ma', type=float, required=True, help='ALP mass in GeV')
    parser.add_argument('--fa', type=float, required=True, help='f_a in GeV')
    parser.add_argument('--fit-path', help='preliminary-fit JSON path')
    parser.add_argument('--fit-label', default=DEFAULT_FIT_LABEL)
    arguments = parser.parse_args()
    fitted = read_preliminary_fit_parameters(
        arguments.fit_path, fit_label=arguments.fit_label)
    result = alp_vector_photon_widths(
        arguments.axial, arguments.pseudoscalar, arguments.ma, arguments.fa,
        fit_parameters=fitted)
    print(json.dumps({'fit_label': fitted.label, **result}, indent=2))


if __name__ == '__main__':
    main()
