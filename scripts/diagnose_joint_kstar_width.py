"""Compute the K* width diagnostic at the accepted joint-fit endpoint.

Run: .venv/bin/python scripts/diagnose_joint_kstar_width.py
Masses and widths are in GeV; chi-square is dimensionless.
"""
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import assess_joint_vector_fit_adequacy as adequacy
import audit_shared_even_gv as audit
import fit_babar_kkpi_global_vector_sector as joint_fit

ROOT = Path(__file__).resolve().parents[1]
COMPARISON = ROOT / adequacy.COMPARISON
ENDPOINT = 'basin_probes/even_shared_widths'


def common_width_fit(observed_widths, covariance):
    """Fit a single unconstrained width to two Gaussian observations."""
    observed_widths = np.asarray(observed_widths, dtype=float)
    covariance = np.asarray(covariance, dtype=float)
    unit_vector = np.ones(2)
    inverse_unit = np.linalg.solve(covariance, unit_vector)
    fitted_width = float(observed_widths @ inverse_unit / (unit_vector @ inverse_unit))
    difference_variance = covariance[0, 0] + covariance[1, 1] - 2 * covariance[0, 1]
    minimum_chi2 = float((observed_widths[0] - observed_widths[1])**2 / difference_variance)
    return fitted_width, minimum_chi2


def diagnostic():
    """Reconstruct the endpoint and compare its width prescriptions."""
    comparison = json.loads(COMPARISON.read_text())
    candidate = adequacy.load_accepted_endpoints(comparison, [ENDPOINT])[0]
    parameters = np.asarray(candidate['endpoint'], dtype=float)
    model = joint_fit.GlobalLikelihood(energy_order=2, pair_order=16, angle_order=16)
    point, _, tensor_mass, coupling, widths = model.unpack(parameters)
    rows = model.gv_data['rows'][:2]
    physical_masses = [row['physical_vector_mass_GeV'] for row in rows]
    physical_coefficients = audit.kstar_coefficients(rows, model.alda['physical_masses_gev'],
                                                      model.kaon_flavors, physical_masses)
    physical_widths = coupling**2 * physical_coefficients
    common_input = (widths['Kstarplus'] + widths['Kstar0']) / 2
    reference = dict(model.preliminary, m_rho_GeV=point[0], m_Kstar_GeV=tensor_mass,
                     Gamma_total_omega_GeV=widths['omega'], Gamma_total_phi_GeV=widths['phi'])
    propagator = joint_fit.running_propagator('kstar', tensor_mass**2, 0., tensor_mass,
                                               common_input, reference)
    on_shell_width = float(-(1 / propagator).imag / tensor_mass)
    observed = np.array([model.odd_data['total_widths'][name]['value_gev']
                         for name in ('Kstarplus', 'Kstar0')])
    errors = np.array([model.odd_data['total_widths'][name]['sigma_gev']
                       for name in ('Kstarplus', 'Kstar0')])
    covariance = np.diag(errors**2)
    fitted_width, minimum_chi2 = common_width_fit(observed, covariance)
    return {'endpoint': ENDPOINT, 'units': 'GeV', 'G_V_GeV': float(coupling),
            'fitted_tensor_mass_GeV': float(tensor_mass),
            'baseline_widths_GeV': {name: float(widths[name]) for name in ('Kstarplus', 'Kstar0')},
            'propagator_on_shell_width_GeV': on_shell_width,
            'physical_mass_tree_widths_GeV': dict(zip(('Kstarplus', 'Kstar0'), map(float, physical_widths))),
            'physical_mass_tree_mean_GeV': float(np.mean(physical_widths)),
            'observed_widths_GeV': observed.tolist(), 'covariance_GeV2': covariance.tolist(),
            'best_common_width_GeV': fitted_width, 'common_width_minimum_chi2': minimum_chi2,
            'common_width_nominal_dof': 1}


if __name__ == '__main__':
    print(json.dumps(diagnostic(), indent=2, allow_nan=False))
