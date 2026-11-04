"""Reproduce the shared even-sector and G_V likelihood diagnosis.

Run from the repository root:
  .venv/bin/python scripts/audit_shared_even_gv.py --output data/vector_sector/shared_even_gv_audit.json
  .venv/bin/python scripts/audit_shared_even_gv.py --self-test

Natural units are used: masses and widths are in GeV, branching fractions are
dimensionless, and covariance entries use the corresponding squared units.
The stored joint fit is a nonconverged snapshot. Diagnostic subset fits and
the conditional G_V fit reuse its measurements and are not extra likelihoods.
"""

import argparse
import hashlib
import json
from pathlib import Path
import sys
import unittest

import numpy as np

import fit_babar_kkpi_global_vector_sector as joint_fit
import fit_vector_gv


ROOT = Path(__file__).resolve().parents[1]
BASELINE_PATH = ROOT / 'data/vector_sector/babar_kkpi_global_vector_sector_fit.json'
GV_DATA_PATH = ROOT / 'data/vector_sector/gv_measurements.json'
PDG_DATA_PATH = ROOT / 'data/vector_sector/pdg2024.json'
ODD_DATA_PATH = ROOT / 'data/odd_parity/odd_low_energy_measurements.json'
ALDA_DATA_PATH = ROOT / 'data/vector_sector/alda_2026nrd_arxiv_v1.json'
CONDITIONAL_FIT_PATH = ROOT / 'data/vector_sector/gv_conditional_fit.json'
PDG_EXTRACT_PATH = ROOT / 'memory/extracts/ParticleDataGroup:2024cfk/vector-to-two-pseudoscalar-inputs.md'
BLOCK_BOUNDARIES = (
    ('odd', 0, 12), ('masses_and_leptonic', 12, 19),
    ('shared_widths', 19, 26), ('phi_KK', 26, 28),
    ('BaBar', 28, 73), ('normalization_constraints', 73, 75),
)
KSTAR_MODES = (
    (('K_charged', 'pi0', 'kstar_charged_K+_pi0'),
     ('K_neutral', 'pi_charged', 'kstar_charged_K0_pi+')),
    (('K_neutral', 'pi0', 'kstar_neutral_K0_pi0'),
     ('K_charged', 'pi_charged', 'kstar_neutral_K+_pi-')),
)


def sha256(path):
    """Return the SHA-256 of an input file."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path):
    return json.loads(path.read_text())


def pair_fit(coefficients, observations, covariance):
    """Minimize a Gaussian pair likelihood over one nonnegative scale."""
    coefficients = np.asarray(coefficients, dtype=float)
    observations = np.asarray(observations, dtype=float)
    covariance = np.asarray(covariance, dtype=float)
    inverse_coefficients = np.linalg.solve(covariance, coefficients)
    scale = max(0.0, float(observations @ inverse_coefficients /
                           (coefficients @ inverse_coefficients)))
    differences = scale * coefficients - observations
    whitened = np.linalg.solve(np.linalg.cholesky(covariance), differences)
    return {
        'coefficients': coefficients.tolist(),
        'observations': observations.tolist(),
        'covariance': covariance.tolist(),
        'best_scale': scale,
        'predictions': (scale * coefficients).tolist(),
        'differences': differences.tolist(),
        'whitened_residuals': whitened.tolist(),
        'chi2': float(whitened @ whitened),
        'nominal_degrees_of_freedom': len(observations) - 1,
    }


def block_contributions(residuals):
    """Sum the squared residuals over the stored likelihood boundaries."""
    residuals = np.asarray(residuals, dtype=float)
    if len(residuals) != 75:
        raise ValueError('The saved joint likelihood does not have 75 residuals')
    return {name: {'count': end - start,
                   'chi2': float(residuals[start:end] @ residuals[start:end])}
            for name, start, end in BLOCK_BOUNDARIES}


def kstar_coefficients(rows, physical_masses, flavor_factors, tensor_masses):
    """Calculate charged and neutral summed K-pion widths per G_V squared."""
    coefficients = []
    for row, daughter_modes, tensor_mass in zip(rows, KSTAR_MODES, tensor_masses):
        coefficients.append(sum(
            fit_vector_gv.width_coefficient(
                row['physical_vector_mass_GeV'], tensor_mass,
                physical_masses[first], physical_masses[second],
                flavor_factors[flavor])
            for first, second, flavor in daughter_modes))
    return np.asarray(coefficients)


def assemble_alternative_residuals(model, parameters, state, charged_width, neutral_width):
    """Propagate changed shared K* widths through the fixed joint candidate."""
    point, sector, tensor_mass, coupling, total_widths = state
    alternate_widths = dict(total_widths, Kstarplus=charged_width, Kstar0=neutral_width)
    alternate_state = (point, sector, tensor_mass, coupling, alternate_widths)
    old_residuals, _, _ = model.old_blocks(parameters, alternate_state)
    # The unblocked evaluator caches propagators using parameters alone. Clear
    # that cache because this sensitivity changes widths at fixed parameters.
    model.babar_cache_key = None
    model.babar_cached = {}
    babar_residuals, _ = model.babar_blocks(parameters, alternate_state)
    return np.r_[old_residuals, babar_residuals, parameters[30:32]]


def audit():
    baseline = read_json(BASELINE_PATH)
    measurements = read_json(GV_DATA_PATH)
    pdg = read_json(PDG_DATA_PATH)
    odd = read_json(ODD_DATA_PATH)
    alda = read_json(ALDA_DATA_PATH)
    recorded_hashes = baseline['input_sha256']
    input_hashes = {
        name: {'recorded': recorded, 'current': sha256(ROOT / name),
               'matches': sha256(ROOT / name) == recorded}
        for name, recorded in recorded_hashes.items()
    }
    mismatches = [name for name, result in input_hashes.items() if not result['matches']]
    if mismatches:
        raise ValueError('Baseline input SHA-256 mismatch: ' + ', '.join(mismatches))

    quadrature = baseline['quadrature']
    model = joint_fit.GlobalLikelihood(
        energy_order=quadrature['energy'],
        pair_order=quadrature['pair_per_segment'],
        angle_order=quadrature['angle_per_segment'], bin_block_size=4)
    parameters = np.asarray(baseline['parameter_values'], dtype=float)
    stored_residuals = np.asarray(baseline['residuals'], dtype=float)
    recomputed_residuals = model.residual(parameters)
    if not np.allclose(recomputed_residuals, stored_residuals, rtol=1e-8, atol=1e-8):
        raise ValueError('The saved joint residuals are not reproduced by the current inputs and code')
    recomputed_blocks = block_contributions(recomputed_residuals)
    for name, block in recomputed_blocks.items():
        if not np.isclose(block['chi2'], baseline['blocks'][name]['chi2'], rtol=1e-10, atol=1e-10):
            raise ValueError('Saved block differs from recomputation: ' + name)

    stored_diagnostic = baseline['shared_even_gv_diagnostic']
    recomputed_diagnostic = joint_fit.even_gv_residual_diagnostic(model, parameters)
    stored_diagnostic_residuals = np.asarray(stored_diagnostic['residuals'])
    if not np.isclose(stored_diagnostic_residuals @ stored_diagnostic_residuals,
                      stored_diagnostic['chi2'], rtol=1e-12, atol=1e-12):
        raise ValueError('Stored diagnostic residuals do not reproduce its chi-squared')
    if not np.isclose(recomputed_diagnostic['chi2'], stored_diagnostic['chi2'],
                      rtol=1e-8, atol=1e-8):
        raise ValueError('The recomputed shared even-sector diagnostic differs from the stored value')

    rows = measurements['rows']
    covariance = np.asarray(measurements['covariance'], dtype=float)
    physical_masses = alda['physical_masses_gev']
    state = model.unpack(parameters)
    fitted_tensor_mass = float(state[2])
    coupling = float(state[3])
    flavor_factors = model.kaon_flavors
    baseline_coefficients = kstar_coefficients(rows[:2], physical_masses, flavor_factors,
                                              [fitted_tensor_mass] * 2)
    observations = np.array([row['value'] for row in rows])
    kstar_fit = pair_fit(baseline_coefficients, observations[:2], covariance[:2, :2])
    phi_coefficients = np.array([
        fit_vector_gv.width_coefficient(row['physical_vector_mass_GeV'], 1.0,
                                        physical_masses[daughter], physical_masses[daughter], 1.0)
        for row, daughter in zip(rows[2:], ('K_charged', 'K_neutral'))])
    phi_fit = pair_fit(phi_coefficients, observations[2:], covariance[2:, 2:])
    for result, key in ((kstar_fit, 'exact_Kstar_pair_minimum_chi2'),
                        (phi_fit, 'exact_phi_pair_minimum_chi2')):
        if not np.isclose(result['chi2'], stored_diagnostic[key], rtol=1e-11, atol=1e-11):
            raise ValueError('Pair minimum differs from the stored diagnostic: ' + key)

    provisional, provisional_masses, provisional_angle = fit_vector_gv.fixed_parameters()
    provisional_factors = fit_vector_gv.flavor_factors(provisional_angle)
    conditional_coefficients = list(kstar_coefficients(
        rows[:2], physical_masses, provisional_factors,
        [provisional_masses['kstar_charged'], provisional_masses['kstar_neutral']]))
    for row, daughter, flavor in zip(rows[2:], ('K_charged', 'K_neutral'),
                                     ('phi_K+_K-', 'phi_K0_K0bar')):
        conditional_coefficients.append(fit_vector_gv.width_coefficient(
            row['physical_vector_mass_GeV'], provisional_masses['phi'],
            physical_masses[daughter], physical_masses[daughter],
            provisional_factors[flavor]) / provisional['parameters']['Gamma_total_phi_GeV'])
    conditional_fit = fit_vector_gv.fit_squared_coupling(
        conditional_coefficients, observations, covariance)
    stored_conditional_coupling = read_json(CONDITIONAL_FIT_PATH)['parameters']['G_V_GeV']
    if not np.isclose(conditional_fit['G_V_GeV'], stored_conditional_coupling,
                      rtol=1e-12, atol=1e-12):
        raise ValueError('Recomputed conditional coupling differs from the stored fit')
    conditional_fit.update({
        'scope': 'Four-row fit at fixed provisional masses, mixing angle and phi total width; its data overlap the joint likelihood.',
        'overlaps_joint_likelihood': True,
        'added_to_joint_chi2': False,
        'coefficients': conditional_coefficients,
        'observations': observations.tolist(),
        'covariance': covariance.tolist(),
        'stored_G_V_GeV': stored_conditional_coupling,
    })

    # Sensitivity to using each measured on-shell K* mass in the external
    # tensor normalization, while preserving physical daughter masses and
    # the baseline charge-mean internal propagator prescription.
    physical_vector_masses = [row['physical_vector_mass_GeV'] for row in rows[:2]]
    alternate_coefficients = kstar_coefficients(rows[:2], physical_masses, flavor_factors,
                                               physical_vector_masses)
    alternate_pair = pair_fit(alternate_coefficients, observations[:2], covariance[:2, :2])
    alternate_widths = coupling**2 * alternate_coefficients
    alternate_residuals = assemble_alternative_residuals(model, parameters, state,
                                                        *alternate_widths)
    alternate_blocks = block_contributions(alternate_residuals)
    baseline_chi2 = float(recomputed_residuals @ recomputed_residuals)
    alternate_chi2 = float(alternate_residuals @ alternate_residuals)
    alternate_result = {
        'scope': 'Fixed joint candidate sensitivity, not a refit or a complete isospin-breaking model.',
        'prescription': 'Each physical K* mass replaces the common fitted tensor mass only in its external K-pion width normalization; the charge-mean internal propagator prescription is retained and its shared width is recomputed.',
        'fixed': 'All stored joint parameter coordinates, physical phase-space masses, data and covariance.',
        'free_in_pair_minimum': 'One nonnegative common G_V squared; no joint parameters are refitted.',
        'tensor_masses_GeV': physical_vector_masses,
        'shared_widths_GeV': alternate_widths.tolist(),
        'joint_chi2_at_fixed_candidate': alternate_chi2,
        'delta_joint_chi2_at_fixed_candidate': alternate_chi2 - baseline_chi2,
        'blocks_at_fixed_candidate': alternate_blocks,
        'delta_blocks_at_fixed_candidate': {name: alternate_blocks[name]['chi2'] - block['chi2']
                                            for name, block in recomputed_blocks.items()},
        'kstar_pair_minimum': alternate_pair,
        'delta_kstar_pair_minimum_chi2': alternate_pair['chi2'] - kstar_fit['chi2'],
    }

    even_rows = []
    for category, names in (('masses', ('rho_charged', 'kstar_charged', 'omega', 'phi')),
                            ('leptonic_widths', ('rho', 'omega', 'phi'))):
        for name in names:
            even_rows.append({'category': category, 'name': name, **pdg[category][name]})
    source_hashes = {}
    for source in measurements['sources'].values():
        path = ROOT / source['local_path']
        source_hashes[source['local_path']] = {
            'recorded': source['sha256'], 'current': sha256(path) if path.exists() else None,
            'matches': sha256(path) == source['sha256'] if path.exists() else None,
        }
    provenance_paths = [BASELINE_PATH, GV_DATA_PATH, CONDITIONAL_FIT_PATH, PDG_DATA_PATH, ODD_DATA_PATH,
                        ALDA_DATA_PATH, PDG_EXTRACT_PATH,
                        ROOT / 'data/vector_sector/odd_preliminary_fit.json',
                        ROOT / 'data/vector_sector/stage3_results.json',
                        ROOT / 'scripts/fit_vector_gv.py',
                        ROOT / 'scripts/fit_babar_kkpi_global_vector_sector.py', Path(__file__).resolve()]
    return {
        'scope': 'Reproducible diagnostic of one stored nonconverged joint-fit snapshot; subset fits are not additional likelihood data.',
        'units': 'Masses and widths GeV; branching fractions dimensionless; covariance uses corresponding squared units.',
        'bibtex_keys': baseline['bibtex_keys'] + ['Ecker:1988te', 'Ecker:1989yg'],
        'input_sha256': {str(path.relative_to(ROOT)): sha256(path) for path in provenance_paths},
        'baseline_recorded_input_sha256': input_hashes,
        'source_archive_sha256': source_hashes,
        'input_inventory': {
            'even_rows': even_rows,
            'even_sources': pdg['sources'],
            'shared_total_width_rows': odd['total_widths'],
            'shared_width_sources': odd['sources'],
            'conditional_G_V_rows': rows,
            'conditional_G_V_sources': measurements['sources'],
            'conditional_G_V_covariance': measurements['covariance'],
            'covariance_assumption': measurements['covariance_assumption'],
            'physical_daughter_masses_GeV': physical_masses,
            'Kstar_daughter_modes': KSTAR_MODES,
            'fitted_common_Kstar_tensor_mass_GeV': fitted_tensor_mass,
            'physical_vector_masses_GeV': [row['physical_vector_mass_GeV'] for row in rows],
            'joint_candidate_G_V_GeV': coupling,
            'joint_shared_widths_GeV': state[4],
            'joint_Kstar_internal_width_GeV': (state[4]['Kstarplus'] + state[4]['Kstar0']) / 2,
            'duplicate_policy': 'Kstar Kpi row values and errors match shared total-width rows; each appears once in shared_widths. Phi branching rows occur once in phi_KK. The four-row conditional diagnostic overlaps these and is not added.',
            'saturation': 'The joint model sets each Kstar total width equal to its predicted K-pion partial width; the data summaries are not counted twice.',
        },
        'baseline': {
            'chi2': baseline['chi2'], 'dof': baseline['dof'],
            'success': baseline['success'], 'message': baseline['message'], 'nfev': baseline['nfev'],
            'quadrature': quadrature, 'assumptions': baseline['assumptions'],
            'stored_blocks': baseline['blocks'], 'recomputed_blocks': recomputed_blocks,
            'maximum_residual_difference': float(np.max(np.abs(recomputed_residuals - stored_residuals))),
            'shared_even_diagnostic_stored': stored_diagnostic,
            'shared_even_diagnostic_recomputed': recomputed_diagnostic,
            'diagnostic_maximum_residual_difference': float(np.max(np.abs(
                np.asarray(recomputed_diagnostic['residuals']) - stored_diagnostic_residuals))),
        },
        'pair_minima': {
            'scope': 'Separate one-scale diagnostic minima; they overlap the joint likelihood and are not additional chi-squared contributions.',
            'Kstar': {'free_scale': 'nonnegative G_V squared, with the common tensor mass fixed', **kstar_fit},
            'phi': {'free_scale': 'nonnegative common normalization containing G_V squared, phi flavor factor, tensor mass and total width', **phi_fit},
        },
        'conditional_four_row_fit': conditional_fit,
        'comparisons': {
            'physical_tensor_mass_sensitivity': alternate_result,
            'non_saturated_Kstar_likelihood': {
                'status': 'unavailable',
                'source': str(PDG_EXTRACT_PATH.relative_to(ROOT)),
                'branching_fractions': {
                    'Kstarplus': {'value': 0.99902, 'sigma': 0.00009, 'central_deficit_from_one': 0.00098},
                    'Kstar0': {'value': 0.99754, 'sigma': 0.00021, 'central_deficit_from_one': 0.00246},
                },
                'reason': 'These source fractions and the partial/total width summaries arise from overlapping constrained fits; their cross-covariance is unavailable. They are not independent added measurements.',
            },
        },
    }


class PairFitTests(unittest.TestCase):
    """Check the two-observation Gaussian minimum against its closed form."""

    def test_correlated_two_observation_minimum(self):
        coefficients = np.array([2.0, 3.0])
        observations = np.array([1.0, 2.0])
        covariance = np.array([[4.0, 1.0], [1.0, 9.0]])
        result = pair_fit(coefficients, observations, covariance)
        # For a two-row Gaussian, profiling a common scale gives this
        # contrast divided by its propagated variance.
        contrast = coefficients[1] * observations[0] - coefficients[0] * observations[1]
        variance = (coefficients[1]**2 * covariance[0, 0]
                    + coefficients[0]**2 * covariance[1, 1]
                    - 2 * coefficients[0] * coefficients[1] * covariance[0, 1])
        self.assertAlmostEqual(result['chi2'], contrast**2 / variance, places=14)

    def test_nonnegative_scale_boundary(self):
        result = pair_fit([1.0, 2.0], [-1.0, -1.0], np.eye(2))
        self.assertEqual(result['best_scale'], 0.0)
        self.assertAlmostEqual(result['chi2'], 2.0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, help='Path for the new audit JSON')
    parser.add_argument('--self-test', action='store_true', help='Run the analytic pair-fit unit checks')
    arguments = parser.parse_args()
    if arguments.self_test:
        result = unittest.main(argv=[__file__], exit=False)
        if not result.result.wasSuccessful():
            sys.exit(1)
        return
    if arguments.output is None:
        parser.error('--output is required unless --self-test is supplied')
    forbidden_paths = {BASELINE_PATH, GV_DATA_PATH, CONDITIONAL_FIT_PATH, PDG_DATA_PATH, ODD_DATA_PATH,
                       ALDA_DATA_PATH, PDG_EXTRACT_PATH, Path(__file__).resolve(),
                       ROOT / 'data/vector_sector/odd_preliminary_fit.json',
                       ROOT / 'data/vector_sector/stage3_results.json',
                       ROOT / 'scripts/fit_vector_gv.py',
                       ROOT / 'scripts/fit_babar_kkpi_global_vector_sector.py'}
    baseline = read_json(BASELINE_PATH)
    forbidden_paths.update(ROOT / name for name in baseline['input_sha256'])
    forbidden_paths.update(ROOT / source['local_path'] for source in read_json(GV_DATA_PATH)['sources'].values())
    if arguments.output.resolve() in {path.resolve() for path in forbidden_paths}:
        parser.error('The output path must not overwrite an input or the stored fit')
    result = audit()
    arguments.output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(arguments.output)


if __name__ == '__main__':
    main()
