"""ALP decays to omega omega, omega phi and phi phi, with masses in GeV.

Only NumPy, the standard library and the preliminary fit JSON are needed.
The quark couplings are supplied directly in the no-mixing basis, in u,d,s
order. No transformation of these couplings is performed.

Example (from the repository root)::

    from scripts.implementation.alp_vv_width import (
        read_preliminary_fit_parameters, alp_vector_pair_widths)
    fit = read_preliminary_fit_parameters()
    result = alp_vector_pair_widths(
        axial_quark_couplings=[1., 0., 0.],
        pseudoscalar_quark_couplings=[0., 0., 0.],
        alp_mass=2.5, alp_decay_constant=1000., fit_parameters=fit)

The result contains scalar amplitudes in GeV^-1 and widths in GeV. The
calculation retains the original eleven odd operators, the fitted one-trace
even electromagnetic source and both photon-to-vector conversion diagrams.
Daughter vectors are stable on-shell particles in the phase space.
"""
from dataclasses import dataclass
import json
import math
from pathlib import Path

import numpy as np


DEFAULT_FIT_PATH = (Path(__file__).resolve().parents[2]
                    / 'data/vector_sector/odd_preliminary_fit.json')
DEFAULT_FIT_LABEL = 'central_single_trace_negative'
VECTOR_PAIRS = (('omega', 'omega'), ('omega', 'phi'), ('phi', 'phi'))

# Fixed central inputs of the preliminary fit, not additional fit parameters.
# B0 m_q: Alda:2026nrd, data/vector_sector/alda_arxiv_v1.json, arXiv v1.
# chi_plus is the exact-isospin even-sector value 4 B0 m_q; the odd
# insertion retains the separate u and d quark masses, as in the fit.
# alpha_em: data/vector_sector/stage3_results.json, assumptions.alpha_em.
# M_V: fixed odd-operator normalization, odd-preliminary-fit.md (GeV).
FIXED_CALCULATION_INPUTS = {
    'B0mq_GeV2': (0.006568, 0.011740, 0.235921),
    'chi_plus_diagonal_GeV2': (0.036615999999999996,
                              0.036615999999999996, 0.943684),
    'alpha_em': 0.0072973525692838015,
    'operator_mass_GeV': 0.770,
}


@dataclass(frozen=True)
class PreliminaryFitParameters:
    """Inputs for the amplitudes, with physical quantities named explicitly.

    ``odd_coefficients`` contains dimensionless c1,...,c7,d1,...,d4 in the
    canonical epsilon^0123=+1 convention. ``even_parameters`` contains the
    fitted mass and electromagnetic-source coefficients with the units in
    their names. Vector flavour generators obey Tr(T_i T_j)=delta_ij;
    they are derived matrices, rather than fitted nuisance parameters.
    ``photon_vector_residues`` are signed electromagnetic-current matrix
    elements in GeV. The quark mass matrix is chi=2 B0 m_q, in GeV^2.
    """
    label: str
    odd_coefficients: dict
    even_parameters: dict
    vector_masses: dict
    vector_generators: dict
    photon_vector_residues: dict
    quark_mass_matrix: np.ndarray
    alpha_em: float
    operator_mass: float
    shared_total_widths: dict


def _positive_number(value, name):
    if isinstance(value, (bool, complex)) or not np.isscalar(value):
        raise ValueError(name + ' must be a positive finite real number')
    value = float(value)
    if not math.isfinite(value) or value <= 0:
        raise ValueError(name + ' must be a positive finite real number')
    return value


def read_preliminary_fit_parameters(fit_path=None, *, fit_label=DEFAULT_FIT_LABEL):
    """Read one one-trace preliminary-fit solution and its fixed inputs.

    Fixed central Alda quark masses, the exact-isospin even mass source,
    alpha_em and M_V are pinned in ``FIXED_CALCULATION_INPUTS`` with their
    provenance. An optional JSON ``calculation_inputs`` object can override
    these explicit inputs. The result gives original odd coefficients, named even fit
    parameters, the derived vector masses/generators and electromagnetic
    residues. Shared measured total widths are retained separately; these
    widths do not enter the stable-daughter two-body formula.
    """
    document = json.loads(Path(fit_path or DEFAULT_FIT_PATH).read_text())
    selected = [entry for entry in document['results'] if entry['label'] == fit_label]
    if len(selected) != 1:
        raise ValueError('fit_label must select exactly one preliminary fit')
    selected = selected[0]
    if selected.get('model') != 'original_11':
        raise ValueError('the fit must use the original eleven odd operators')
    if selected['scenario']['branch'] != 'near_ideal':
        raise ValueError('this reader requires the near-ideal neutral vector mass branch')
    hypothesis = selected['scenario']['source_hypothesis']
    if hypothesis not in ('single_trace_positive', 'single_trace_negative'):
        raise ValueError('this calculation requires a one-trace electromagnetic source fit')
    fixed = dict(FIXED_CALCULATION_INPUTS)
    fixed.update(document.get('calculation_inputs', {}))
    even_parameters = {name: float(value) for name, value
                       in selected['even_parameters_by_name'].items()}
    rho_mass = even_parameters['m_rho_GeV']
    kstar_mass = even_parameters['m_Kstar_GeV']
    singlet_mass_coefficient = even_parameters['Delta_V_GeV2']
    singlet_mass_breaking_coefficient = even_parameters['kappa_V']
    chi_plus_diagonal = np.asarray(fixed['chi_plus_diagonal_GeV2'], dtype=float)
    if chi_plus_diagonal.shape != (3,) or not np.all(np.isfinite(chi_plus_diagonal)):
        raise ValueError('chi_plus_diagonal_GeV2 must have three finite entries')
    if chi_plus_diagonal[0] != chi_plus_diagonal[1]:
        raise ValueError('the fitted even vector sector uses exact isospin')
    light_mass_source, strange_mass_source = chi_plus_diagonal[[0, 2]]
    # Eliminate the common mass and octet splitting with the rho/Kstar masses.
    neutral_mass_squared_matrix = np.array([
        [rho_mass**2 + 2*singlet_mass_coefficient/3
         + 2*singlet_mass_breaking_coefficient*light_mass_source,
         math.sqrt(2)*singlet_mass_coefficient/3
         + singlet_mass_breaking_coefficient*(light_mass_source+strange_mass_source)/math.sqrt(2)],
        [math.sqrt(2)*singlet_mass_coefficient/3
         + singlet_mass_breaking_coefficient*(light_mass_source+strange_mass_source)/math.sqrt(2),
         2*kstar_mass**2-rho_mass**2 + singlet_mass_coefficient/3
         + 2*singlet_mass_breaking_coefficient*strange_mass_source]])
    vector_mass_squared = np.linalg.eigvalsh(neutral_mass_squared_matrix)
    if np.any(vector_mass_squared <= 0):
        raise ValueError('fitted neutral vector mass matrix must be positive')
    vector_mixing_angle = 0.5*math.atan2(
        2*neutral_mass_squared_matrix[0, 1],
        neutral_mass_squared_matrix[1, 1]-neutral_mass_squared_matrix[0, 0])
    ideal_omega_generator = np.diag([1., 1., 0.])/math.sqrt(2)
    ideal_phi_generator = np.diag([0., 0., 1.])
    vector_generators = {
        'rho0': np.diag([1., -1., 0.])/math.sqrt(2),
        'omega': math.cos(vector_mixing_angle)*ideal_omega_generator
                 - math.sin(vector_mixing_angle)*ideal_phi_generator,
        'phi': math.sin(vector_mixing_angle)*ideal_omega_generator
               + math.cos(vector_mixing_angle)*ideal_phi_generator}
    vector_masses = dict(zip(('rho0', 'omega', 'phi'),
                             (rho_mass, *np.sqrt(vector_mass_squared))))
    electric_charges = np.diag([2/3, -1/3, -1/3])
    even_electromagnetic_source = (
        even_parameters['F_V_GeV']*electric_charges
        + 2*even_parameters['z_V_GeV_inverse']*np.diag(chi_plus_diagonal)@electric_charges)
    photon_vector_residues = {
        name: float(np.trace(generator@even_electromagnetic_source))
        for name, generator in vector_generators.items()}
    quark_masses = np.asarray(fixed['B0mq_GeV2'], dtype=float)
    if quark_masses.shape != (3,) or not np.all(np.isfinite(quark_masses)):
        raise ValueError('B0mq_GeV2 must have three finite entries')
    odd_coefficients = {name: float(value) for name, value
                        in selected['canonical_representative_coefficients'].items()}
    return PreliminaryFitParameters(
        fit_label, odd_coefficients, even_parameters, vector_masses,
        vector_generators, photon_vector_residues, np.diag(2*quark_masses),
        _positive_number(fixed['alpha_em'], 'alpha_em'),
        _positive_number(fixed['operator_mass_GeV'], 'operator_mass_GeV'),
        dict(selected['shared_total_widths_gev']))


def _quark_couplings(values, name):
    couplings = np.asarray(values)
    if (couplings.shape != (3,) or np.iscomplexobj(couplings)
            or not np.all(np.isfinite(couplings))):
        raise ValueError(name + ' must contain three finite real entries (u,d,s)')
    return couplings.astype(float)


def alp_vector_pair_amplitudes(axial_quark_couplings, pseudoscalar_quark_couplings,
                               alp_mass, alp_decay_constant, *, fit_parameters=None):
    """Return G_ij for M_ij=G_ij epsilon eps_i eps_j q_i q_j in GeV^-1.

    Couplings are the dimensionless c_q^A,c_q^P in the no-mixing basis,
    each in u,d,s order. Mass and decay constant are in GeV. ``fit_parameters``
    is returned by ``read_preliminary_fit_parameters``; omitting it reads the default
    negative-strange one-trace fit. The formula exposes VVP and VJP operator
    contractions directly and retains both photon-conversion assignments.

    With s=M_i^2, t=M_j^2 and q_i.q_j=(m_a^2-s-t)/2, the VVP numerator is
    4[d1*m_a^2-d3*(m_a^2-s-t)] Tr(T_i T_j derivative_source)
    -16*d2*Tr(T_i T_j mass_source). The ordered VJP numerator is
    (4/M_V){[c1*(t+q_i.q_j)+c2*(s+q_i.q_j)-c5*q_i.q_j-c6*s]
    Tr(T_i T_j derivative_source)-2*c3*Tr(T_i T_j mass_source)}.
    These contractions assume diagonal neutral matrices. c4,c7,d4 vanish.
    """
    alp_mass = _positive_number(alp_mass, 'alp_mass')
    alp_decay_constant = _positive_number(alp_decay_constant, 'alp_decay_constant')
    fitted = fit_parameters if fit_parameters is not None else read_preliminary_fit_parameters()
    axial_quark_couplings = _quark_couplings(axial_quark_couplings, 'axial_quark_couplings')
    pseudoscalar_quark_couplings = _quark_couplings(
        pseudoscalar_quark_couplings, 'pseudoscalar_quark_couplings')
    alp_derivative_source_matrix = np.diag(2*axial_quark_couplings/alp_decay_constant)
    pseudoscalar_coupling_matrix = np.diag(pseudoscalar_quark_couplings)
    alp_mass_source_matrix = -2/alp_decay_constant*(
        pseudoscalar_coupling_matrix@fitted.quark_mass_matrix
        + fitted.quark_mass_matrix@pseudoscalar_coupling_matrix)
    odd = fitted.odd_coefficients
    electric_charge_squared = 4*math.pi*fitted.alpha_em
    amplitudes = {}
    for first_vector, second_vector in VECTOR_PAIRS:
        first_mass = fitted.vector_masses[first_vector]
        second_mass = fitted.vector_masses[second_vector]
        vector_generator_product = (fitted.vector_generators[first_vector]
                                    @fitted.vector_generators[second_vector])
        derivative_source_trace = float(np.trace(
            vector_generator_product@alp_derivative_source_matrix))
        mass_source_trace = float(np.trace(vector_generator_product@alp_mass_source_matrix))
        momentum_scalar_product = (alp_mass**2-first_mass**2-second_mass**2)/2
        # Original d1,d2,d3 VVP terms; d4 vanishes for neutral on-shell vectors.
        two_vector_contribution = (
            4*(odd['d1']*alp_mass**2
               - 2*odd['d3']*momentum_scalar_product)*derivative_source_trace
            - 16*odd['d2']*mass_source_trace)

        def vector_photon_contribution(vector_mass_squared, photon_mass_squared):
            # c4 is a neutral commutator; c7 has zero on-shell contraction.
            return 4/fitted.operator_mass*(derivative_source_trace*(
                odd['c1']*(photon_mass_squared+momentum_scalar_product)
                + odd['c2']*(vector_mass_squared+momentum_scalar_product)
                - odd['c5']*momentum_scalar_product
                - odd['c6']*vector_mass_squared) - 2*odd['c3']*mass_source_trace)

        photon_conversion_contribution = -2*math.sqrt(2)*electric_charge_squared*(
            fitted.photon_vector_residues[second_vector]
            * vector_photon_contribution(first_mass**2, second_mass**2)
            + fitted.photon_vector_residues[first_vector]
            * vector_photon_contribution(second_mass**2, first_mass**2))
        amplitudes[first_vector+second_vector] = (
            two_vector_contribution+photon_conversion_contribution)/(first_mass*second_mass)
    return amplitudes


def alp_vector_pair_widths(axial_quark_couplings, pseudoscalar_quark_couplings,
                           alp_mass, alp_decay_constant, *, fit_parameters=None):
    """Return amplitudes in GeV^-1 and narrow-vector partial widths in GeV.

    Below and at threshold the width is zero. The omega omega and phi phi
    channels include 1/2! for identical particles. No daughter spectral
    integration or uncertainty propagation is included.
    """
    alp_mass = _positive_number(alp_mass, 'alp_mass')
    fitted = fit_parameters if fit_parameters is not None else read_preliminary_fit_parameters()
    amplitudes = alp_vector_pair_amplitudes(
        axial_quark_couplings, pseudoscalar_quark_couplings,
        alp_mass, alp_decay_constant, fit_parameters=fitted)
    widths = {}
    for first_vector, second_vector in VECTOR_PAIRS:
        channel = first_vector+second_vector
        first_mass = fitted.vector_masses[first_vector]
        second_mass = fitted.vector_masses[second_vector]
        if alp_mass <= first_mass+second_mass:
            widths[channel] = 0.
            continue
        phase_space_polynomial = ((alp_mass**2-(first_mass+second_mass)**2)
                                  *(alp_mass**2-(first_mass-second_mass)**2))
        identical_particle_factor = 2 if first_vector == second_vector else 1
        widths[channel] = abs(amplitudes[channel])**2*phase_space_polynomial**1.5/(
            32*math.pi*alp_mass**3*identical_particle_factor)
    return {'amplitudes_GeV_inv': amplitudes, 'widths_GeV': widths}
