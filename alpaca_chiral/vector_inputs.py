"""Central RChT inputs from the original-eleven-operator preliminary fit.

Masses are in GeV. The single-trace negative-strange solution, its fixed
inputs, and its limited calibration scope follow odd-preliminary-fit.md.
"""
import json
import math
from functools import lru_cache
from pathlib import Path
import numpy as np

DEFAULT_FIT_PATH = Path(__file__).with_name('odd_preliminary_fit.json')
DEFAULT_FIT_LABEL = 'central_single_trace_negative'

FIXED_CALCULATION_INPUTS = {
    'B0mq_GeV2': (0.006568, 0.011740, 0.235921),
    'chi_plus_diagonal_GeV2': (0.036615999999999996,
                              0.036615999999999996, 0.943684),
    'alpha_em': 0.0072973525692838015,
    'operator_mass_GeV': 0.770,
}


@lru_cache(maxsize=1)
def _read_vector_inputs():
    """Reconstruct masses, normalized generators, and signed source residues."""
    fit_path = DEFAULT_FIT_PATH
    fit_label = DEFAULT_FIT_LABEL
    document = json.loads(Path(fit_path or DEFAULT_FIT_PATH).read_text())
    selected = [entry for entry in document['results'] if entry['label'] == fit_label]
    selected = selected[0]
    hypothesis = selected['scenario']['source_hypothesis']
    fixed = dict(FIXED_CALCULATION_INPUTS)
    fixed.update(document.get('calculation_inputs', {}))
    even_parameters = {name: float(value) for name, value
                       in selected['even_parameters_by_name'].items()}
    rho_mass = even_parameters['m_rho_GeV']
    kstar_mass = even_parameters['m_Kstar_GeV']
    singlet_mass_coefficient = even_parameters['Delta_V_GeV2']
    singlet_mass_breaking_coefficient = even_parameters['kappa_V']
    chi_plus_diagonal = np.asarray(fixed['chi_plus_diagonal_GeV2'], dtype=float)
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
    odd_coefficients = {name: float(value) for name, value
                        in selected['canonical_representative_coefficients'].items()}
    return dict(zip(
        ('label', 'odd_coefficients', 'even_parameters', 'vector_masses',
         'vector_generators', 'photon_vector_residues', 'quark_mass_matrix',
         'alpha_em', 'operator_mass', 'shared_total_widths'),
        (fit_label, odd_coefficients, even_parameters, vector_masses,
         vector_generators, photon_vector_residues, np.diag(2*quark_masses),
         fixed['alpha_em'], fixed['operator_mass_GeV'],
         dict(selected['shared_total_widths_gev']))))
