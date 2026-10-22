"""Conditional LO RChT G_V fit with the provisional parameters held fixed.

Run: .venv/bin/python scripts/fit_vector_gv.py measurements.json output.json
The output contains only metadata.bibtex_references and parameters.G_V_GeV.
Regenerate this file before running scripts/summarize_odd_preliminary_fit.py.
Input rows specify vector, channel, observable and value; covariance uses
GeV for partial widths and dimensionless branching fractions. K*pi sums modes.
Physical external masses enter phase space; fitted tensor masses enter the
on-shell tensor-to-Proca conversion. Natural units, F=0.093 GeV, no loops,
no electromagnetic corrections and no short-distance constraint are assumed.
Widths determine |G_V| only. The reported positive branch is a convention.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
FIT_PATH = ROOT / 'data/vector_sector/odd_preliminary_fit.json'
ALDA_PATH = ROOT / 'data/vector_sector/alda_2026nrd_arxiv_v1.json'
DECAY_CONSTANT_GEV = 0.093


def matrix_unit(row, column):
    """Return a flavor matrix with one unit entry."""
    result = sp.zeros(3)
    result[row, column] = 1
    return result


def flavor_factors(mixing_angle):
    """Return squared LO VPP flavor traces, summed over stated final states."""
    pion_plus = matrix_unit(0, 1)
    pion_minus = matrix_unit(1, 0)
    kaon_plus = matrix_unit(0, 2)
    kaon_minus = matrix_unit(2, 0)
    kaon_zero = matrix_unit(1, 2)
    antikaon_zero = matrix_unit(2, 1)
    neutral_inputs = json.loads(ALDA_PATH.read_text())
    pion_zero = sp.diag(*neutral_inputs['published_lambda_diagonals']['pi0'])/sp.sqrt(2)
    kstar_minus = matrix_unit(2, 0)
    antikstar_zero = matrix_unit(2, 1)
    rho = sp.diag(1, -1, 0) / sp.sqrt(2)
    omega_ideal = sp.diag(1, 1, 0) / sp.sqrt(2)
    phi_ideal = sp.diag(0, 0, 1)
    phi = sp.sin(mixing_angle)*omega_ideal + sp.cos(mixing_angle)*phi_ideal
    omega = sp.cos(mixing_angle)*omega_ideal - sp.sin(mixing_angle)*phi_ideal
    def squared_trace(vector, first, second):
        return float(sp.N(sp.trace(vector*(first*second-second*first))**2))
    return {
        'kstar_charged_K+_pi0': squared_trace(kstar_minus, kaon_plus, pion_zero),
        'kstar_charged_K0_pi+': squared_trace(kstar_minus, kaon_zero, pion_plus),
        'kstar_neutral_K0_pi0': squared_trace(antikstar_zero, kaon_zero, pion_zero),
        'kstar_neutral_K+_pi-': squared_trace(antikstar_zero, kaon_plus, pion_minus),
        'phi_pi+_pi-': squared_trace(phi, pion_plus, pion_minus),
        'omega_pi+_pi-': squared_trace(omega, pion_plus, pion_minus),
        'rho0_pi+_pi-': squared_trace(rho, pion_plus, pion_minus),
        'phi_K+_K-': squared_trace(phi, kaon_plus, kaon_minus),
        'phi_K0_K0bar': squared_trace(phi, kaon_zero, antikaon_zero),
    }


def width_coefficient(physical_mass, tensor_mass, first_mass, second_mass,
                      squared_flavor_trace, decay_constant=DECAY_CONSTANT_GEV):
    """Return Gamma/G_V² in GeV^-1 for an allowed two-body channel.

    The derivative vertex gives g_VPP=G_V M_physical² trace/(sqrt(2) F²
    M_tensor). The spin average gives Gamma=g_VPP² p³/(6 pi M_physical²).
    """
    if physical_mass <= first_mass + second_mass:
        return 0.0
    momentum = math.sqrt((physical_mass**2-(first_mass+second_mass)**2)*
                         (physical_mass**2-(first_mass-second_mass)**2))/(2*physical_mass)
    return squared_flavor_trace*physical_mass**2*momentum**3/(12*math.pi*
            decay_constant**4*tensor_mass**2)


def fit_squared_coupling(coefficients, widths, covariance):
    """Fit nonnegative G_V² analytically using a fixed Gaussian covariance.

    Return unscaled experimental curvature and exact Delta chi²=1 limits
    on the nonnegative G_V² branch. No fit-quality rescaling is performed.
    """
    coefficients, widths = np.asarray(coefficients), np.asarray(widths)
    inverse_times_coefficient = np.linalg.solve(covariance, coefficients)
    curvature = float(coefficients @ inverse_times_coefficient)
    unconstrained = float(widths @ inverse_times_coefficient / curvature)
    fitted = max(0.0, unconstrained)
    residual = coefficients*fitted-widths
    chi_squared = float(residual @ np.linalg.solve(covariance, residual))
    radius = math.sqrt(1/curvature + (fitted-unconstrained)**2)
    lower, upper = max(0., unconstrained-radius), unconstrained+radius
    return dict(G_V_GeV=math.sqrt(fitted), G_V_delta_chi2_one_GeV=[math.sqrt(lower),math.sqrt(upper)],
                G_V_squared_GeV2=fitted, chi2=chi_squared,
                nominal_dof=len(widths)-1, predictions=(coefficients*fitted).tolist(),
                residuals=residual.tolist())


def fixed_parameters():
    """Read the selected provisional fit and reconstruct its neutral masses."""
    document = json.loads(FIT_PATH.read_text())
    selected = document
    parameters = selected['parameters']
    rho_mass, kaon_mass = parameters['m_rho_GeV'], parameters['m_Kstar_GeV']
    singlet_shift, breaking = parameters['Delta_V_GeV2'], parameters['kappa_V']
    light, strange = 0.036616, 0.943684
    mass_matrix = np.array([[rho_mass**2+2*singlet_shift/3+2*breaking*light,
        math.sqrt(2)*singlet_shift/3+breaking*(light+strange)/math.sqrt(2)],
        [math.sqrt(2)*singlet_shift/3+breaking*(light+strange)/math.sqrt(2),
         2*kaon_mass**2-rho_mass**2+singlet_shift/3+2*breaking*strange]])
    mixing_angle = parameters['omega_phi_mixing_angle_rad']
    masses = dict(zip(('rho0','omega','phi'),[rho_mass,*np.sqrt(np.linalg.eigvalsh(mass_matrix))]))
    masses['kstar_charged'] = kaon_mass
    masses['kstar_neutral'] = kaon_mass
    return selected, masses, mixing_angle


def main():
    """Fit explicitly sourced measurements without modifying provisional inputs."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('measurements',type=Path)
    parser.add_argument('output',type=Path)
    arguments = parser.parse_args()
    measurements = json.loads(arguments.measurements.read_text())
    selected, tensor_masses, mixing_angle = fixed_parameters()
    factors = flavor_factors(mixing_angle)
    for source in measurements['sources'].values():
        assert hashlib.sha256((ROOT/source['local_path']).read_bytes()).hexdigest()==source['sha256']
    neutral_inputs = json.loads(ALDA_PATH.read_text())
    physical_masses = neutral_inputs['physical_masses_gev']
    def measurement_coefficient(row):
        # At negligible CP violation, K0 K0bar projects into distinct KL KS
        # with the same width; no identical-particle factor is inserted.
        daughter_modes = {
            'kstar_charged_Kpi': [('K_charged','pi0','kstar_charged_K+_pi0'),
                                 ('K_neutral','pi_charged','kstar_charged_K0_pi+')],
            'kstar_neutral_Kpi': [('K_neutral','pi0','kstar_neutral_K0_pi0'),
                                 ('K_charged','pi_charged','kstar_neutral_K+_pi-')],
            'phi_K+_K-': [('K_charged','K_charged','phi_K+_K-')],
            'phi_KL_KS': [('K_neutral','K_neutral','phi_K0_K0bar')],
        }[row['channel']]
        coefficient = sum(width_coefficient(row['physical_vector_mass_GeV'],
            tensor_masses[row['vector']], physical_masses[first],physical_masses[second],
            factors[channel]) for first,second,channel in daughter_modes)
        if row['observable']=='branching_fraction':
            coefficient /= selected['parameters']['Gamma_total_phi_GeV']
        return coefficient
    coefficients = np.array([measurement_coefficient(row) for row in measurements['rows']])
    values = np.array([row['value'] for row in measurements['rows']])
    covariance = np.array(measurements['covariance'])
    result = fit_squared_coupling(coefficients,values,covariance)
    output = dict(
        metadata=dict(bibtex_references=[
            'Ecker:1988te', 'Ecker:1989yg', 'ParticleDataGroup:2024cfk', 'Alda:2026nrd']),
        parameters=dict(G_V_GeV=result['G_V_GeV']))
    arguments.output.write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
    print(arguments.output)


if __name__=='__main__':
    main()
