"""Inspect the correlated phi-to-kaon residuals at the accepted joint-fit point.

Run: .venv/bin/python scripts/diagnose_joint_phi_kaon_residuals.py
Branching fractions and residuals are dimensionless.
"""
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
ENDPOINT = 'basin_probes/even_shared_widths'


def channel_residuals(cholesky_residuals, covariance):
    """Return physical differences and symmetric-whitened channel components.

    The positive symmetric inverse square root fixes the otherwise arbitrary
    orthogonal rotation between whitened vectors.
    """
    cholesky_residuals = np.asarray(cholesky_residuals, dtype=float)
    covariance = np.asarray(covariance, dtype=float)
    cholesky = np.linalg.cholesky(covariance)
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    inverse_square_root = (eigenvectors / np.sqrt(eigenvalues)) @ eigenvectors.T
    rotation = inverse_square_root @ cholesky
    differences = cholesky @ cholesky_residuals
    return differences, rotation @ cholesky_residuals, rotation


def diagnostic():
    """Read the high-order residual ledger and the likelihood covariance."""
    ledger = json.loads((ROOT / 'data/vector_sector/joint_vector_fit_adequacy.json').read_text())
    endpoint = next(entry for entry in ledger['endpoints'] if entry['id'] == ENDPOINT)
    measurements = json.loads((ROOT / 'data/vector_sector/gv_measurements.json').read_text())
    covariance = np.asarray(measurements['covariance'], dtype=float)[2:, 2:]
    observed = np.array([row['value'] for row in measurements['rows'][2:]])
    selected = [row for row in endpoint['standardized_residuals'] if row['block'] == 'phi_KK']
    if [row['index'] for row in selected] != [26, 27]:
        raise ValueError('Unexpected phi branching-fraction residual ordering')
    cholesky_residuals = np.array([row['pull'] for row in selected])
    differences, rotated, rotation = channel_residuals(cholesky_residuals, covariance)
    chi2 = float(rotated @ rotated)
    stored_chi2 = endpoint['blocks']['phi_KK']['chi2']
    if not np.isclose(chi2, stored_chi2, rtol=0., atol=1e-10):
        raise ValueError('Rotated residuals do not reproduce the stored block chi-square')
    predictions = observed + differences
    # Sum and ratio pulls use the covariance of measured fractions. The ratio
    # error is propagated to first order at the observed central values.
    sum_direction = np.array([1., 1.])
    ratio_gradient = np.array([1 / observed[1], -observed[0] / observed[1]**2])
    sum_sigma = float(np.sqrt(sum_direction @ covariance @ sum_direction))
    ratio_sigma = float(np.sqrt(ratio_gradient @ covariance @ ratio_gradient))
    return {
        'endpoint': ENDPOINT,
        'source_quadrature_orders': [4, 128, 128],
        'channel_order': ['K+K-', 'K0K0bar'],
        'covariance': covariance.tolist(),
        'observed_branching_fractions': observed.tolist(),
        'predicted_branching_fractions': predictions.tolist(),
        'prediction_minus_observation': differences.tolist(),
        'cholesky_whitened_residuals': cholesky_residuals.tolist(),
        'rotation_to_symmetric_whitening': rotation.tolist(),
        'symmetric_whitened_channel_components': rotated.tolist(),
        'rotated_chi2': chi2,
        'stored_block_chi2': stored_chi2,
        'sum_residual_in_standard_deviations': float(sum(differences) / sum_sigma),
        'predicted_to_observed_ratio_difference_in_standard_deviations_linearized':
            float((predictions[0] / predictions[1] - observed[0] / observed[1]) / ratio_sigma),
    }


if __name__ == '__main__':
    print(json.dumps(diagnostic(), indent=2, allow_nan=False))
