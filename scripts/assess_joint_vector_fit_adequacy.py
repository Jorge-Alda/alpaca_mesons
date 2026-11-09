"""Goodness of fit at the accepted joint vector-sector minima.

Run: .venv/bin/python scripts/assess_joint_vector_fit_adequacy.py --output data/vector_sector/joint_vector_fit_adequacy.json

Every accepted endpoint of data/vector_sector/joint_vector_endpoint_comparison.json
is evaluated at quadrature orders (energy, pair, angle) = (4,128,128) with
four-bin blocks. Units: GeV, nb; residuals are standardized (dimensionless).

Degrees of freedom (stated assumptions):
  * total: number of residuals (73 observed summaries + 2 normalization
    constraints = 75) minus the 32 free optimizer coordinates = 43;
  * per block: the nominal p-value uses the block row count as degrees of
    freedom, ignoring the parameters fitted to that block. These per-block
    p-values are therefore not calibrated; they only rank blocks.
  * no p-value is corrected for bounded or nonlinear parameters.
"""
import argparse
import hashlib
import json
import shlex
import sys
import time
from pathlib import Path

import numpy as np
from scipy.stats import chi2 as chi2_distribution

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fit_babar_kkpi_global_vector_sector as joint_fit
import fit_joint_vector_multistart as multi
from fit_babar_kkpi_g_coefficients import CHANNELS

ROOT = joint_fit.ROOT
COMPARISON = 'data/vector_sector/joint_vector_endpoint_comparison.json'
ORDERS = (4, 128, 128)
BLOCK_SIZE = 4
BLOCKS = multi.BLOCKS
TOTAL_DEGREES_OF_FREEDOM = 75-32
LARGEST_PULL_COUNT = 15
CHI2_MATCH_TOLERANCE = 1e-8
# Energy sub-ranges for the BaBar residual pattern, in GeV (fixed a priori).
ENERGY_RANGES = ((0., 1.50), (1.50, 1.60), (1.60, 1.70), (1.70, 1.78))
EVEN_ROW_LABELS = ('mass rho charged', 'mass Kstar charged', 'mass omega', 'mass phi',
                   'electronic width rho', 'electronic width omega', 'electronic width phi')


def load_accepted_endpoints(comparison, requested=None):
    """Return the accepted candidates; refuse any requested non-accepted one."""
    candidates = {candidate['id']: candidate for candidate in comparison['candidates']}
    accepted = [identifier for identifier in comparison['ranking'] if candidates[identifier]['accepted']]
    if requested is None:
        return [candidates[identifier] for identifier in accepted]
    for identifier in requested:
        if identifier not in candidates:
            raise ValueError(f'{identifier}: not in the endpoint comparison ledger')
        if candidates[identifier]['accepted'] is not True:
            raise ValueError(f'{identifier}: endpoint is not accepted; no goodness-of-fit result is reported')
    return [candidates[identifier] for identifier in requested]


def row_labels(model):
    """Describe each of the 75 residual rows."""
    labels = []
    for row in model.observations:
        labels.append(dict(block='odd', label=f'{row["channel"]} {row["observable"]}'))
    labels += [dict(block='masses_and_leptonic', label=label) for label in EVEN_ROW_LABELS]
    labels += [dict(block='shared_widths', label=f'total width {name}') for name in model.total_order]
    labels.append(dict(block='phi_KK', label='phi K+K- branching fraction (Cholesky-whitened component 1)'))
    labels.append(dict(block='phi_KK', label='phi KLKS branching fraction given K+K- (Cholesky-whitened component 2)'))
    for channel in CHANNELS:
        data = model.babar[channel]
        for (lower, upper), observed, error in zip(data['edges'], data['observed'], data['errors']):
            labels.append(dict(block='BaBar', label=f'BaBar {channel} {lower:.2f}-{upper:.2f} GeV',
                               channel=channel, energy_low_gev=float(lower),
                               energy_high_gev=float(upper),
                               energy_center_gev=float((lower+upper)/2),
                               observed_nb=float(observed), error_nb=float(error)))
    labels += [dict(block='normalization_constraints', label=f'normalization pull {channel}')
               for channel in CHANNELS]
    assert len(labels) == 75
    return labels


def block_summary(residual, blocks=BLOCKS, mask=None):
    """Per-block chi2, row count and nominal p-value (row count as dof)."""
    summary = {}
    for name, start, stop in blocks:
        rows = np.arange(start, stop)
        if mask is not None:
            rows = rows[mask[start:stop]]
        values = residual[rows]
        value = float(values @ values)
        summary[name] = dict(chi2=value, count=int(len(rows)),
                             nominal_p_value_row_count_dof=(float(chi2_distribution.sf(value, len(rows)))
                                                            if len(rows) else None))
    return summary


def runs_statistics(signs):
    """Wald-Wolfowitz runs count of residual signs, ordered in energy."""
    signs = np.asarray(signs)
    positive = int(np.sum(signs > 0))
    negative = int(np.sum(signs < 0))
    runs = int(1+np.sum(signs[1:] != signs[:-1])) if len(signs) else 0
    total = positive+negative
    if positive == 0 or negative == 0:
        return dict(positive=positive, negative=negative, runs=runs, expected_runs=None,
                    standard_deviation=None, standardized_runs=None)
    expected = 1+2*positive*negative/total
    variance = 2*positive*negative*(2*positive*negative-total)/(total**2*(total-1))
    return dict(positive=positive, negative=negative, runs=runs, expected_runs=float(expected),
                standard_deviation=float(np.sqrt(variance)),
                standardized_runs=float((runs-expected)/np.sqrt(variance)) if variance > 0 else None)


def babar_energy_pattern(residual, labels):
    """Per channel: signs ordered in energy, runs, and chi2 in energy sub-ranges."""
    pattern = {}
    for channel in CHANNELS:
        rows = [index for index, label in enumerate(labels) if label.get('channel') == channel]
        values = residual[rows]
        centers = np.array([labels[index]['energy_center_gev'] for index in rows])
        signs = np.sign(values)
        ranges = []
        for lower, upper in ENERGY_RANGES:
            selected = (centers >= lower) & (centers < upper)
            part = values[selected]
            ranges.append(dict(energy_low_gev=lower, energy_high_gev=upper, count=int(selected.sum()),
                               chi2=float(part @ part), mean_pull=float(part.mean()) if len(part) else None,
                               positive=int(np.sum(part > 0)), negative=int(np.sum(part < 0))))
        pattern[channel] = dict(energy_center_gev=centers.tolist(), pulls=values.tolist(),
                                sign_sequence=''.join('+' if sign > 0 else '-' for sign in signs),
                                chi2=float(values @ values), count=len(rows),
                                runs=runs_statistics(signs), energy_ranges=ranges)
    return pattern


def local_errors(jacobian):
    """Linearized errors sqrt(diag((J^T J)^-1)), as in the fit script."""
    column_norms = np.linalg.norm(jacobian, axis=0)
    _, singular, right = np.linalg.svd(jacobian/column_norms, full_matrices=False)
    covariance = ((right.T/np.maximum(singular, 1e-14)**2)@right)/column_norms[:, None]/column_norms[None, :]
    return np.sqrt(np.diag(covariance)), singular


def assess_endpoint(model, candidate, labels, with_local_errors=True):
    """Goodness-of-fit record for one accepted endpoint."""
    started = time.perf_counter()
    values = np.asarray(candidate['endpoint'], dtype=float)
    residual = model.residual(values)
    if residual.shape != (75,) or not np.all(np.isfinite(residual)) or np.all(
            residual == 1e8+np.linalg.norm(values)):
        raise ValueError(f'{candidate["id"]}: invalid residual at the common order')
    total = float(residual @ residual)
    difference = abs(total-candidate['common_chi2'])
    ranked = np.argsort(-np.abs(residual), kind='stable')
    record = dict(id=candidate['id'], endpoint=values.tolist(), chi2=total,
                  comparison_ledger_chi2=candidate['common_chi2'],
                  absolute_chi2_difference_to_comparison=difference,
                  chi2_matches_comparison=bool(difference <= CHI2_MATCH_TOLERANCE),
                  residual_count=75, free_coordinate_count=32,
                  degrees_of_freedom=TOTAL_DEGREES_OF_FREEDOM,
                  nominal_p_value=float(chi2_distribution.sf(total, TOTAL_DEGREES_OF_FREEDOM)),
                  chi2_per_degree_of_freedom=total/TOTAL_DEGREES_OF_FREEDOM,
                  blocks=block_summary(residual),
                  standardized_residuals=[dict(index=index, pull=float(residual[index]), **labels[index])
                                          for index in range(75)],
                  largest_pulls=[dict(rank=rank+1, index=int(index), pull=float(residual[index]),
                                      chi2_contribution=float(residual[index]**2), **labels[index])
                                 for rank, index in enumerate(ranked[:LARGEST_PULL_COUNT])],
                  babar_energy_pattern=babar_energy_pattern(residual, labels))
    if with_local_errors:
        jacobian = model.jacobian(values)
        errors, singular = local_errors(jacobian)
        record['local_errors'] = dict(zip(joint_fit.PARAMETER_NAMES, errors.tolist()))
        record['local_error_vector'] = errors.tolist()
        record['column_normalized_singular_values'] = singular.tolist()
        record['local_error_definition'] = ('sqrt(diag((J^T J)^-1)) from the (4,128,128) Jacobian '
                                            '(analytic odd columns, forward differences otherwise); '
                                            'linearized, not profile intervals')
    record['elapsed_seconds'] = time.perf_counter()-started
    return record


def file_hash(path):
    """SHA-256 of a file given relative to the repository root."""
    return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()


def assess(model, comparison, requested=None, with_local_errors=True):
    """Assess every requested accepted endpoint at the common order."""
    if model.bin_block_size != BLOCK_SIZE or model.phase.orders != ORDERS:
        raise ValueError('assessment requires (4,128,128) with four-bin blocks')
    candidates = load_accepted_endpoints(comparison, requested)
    labels = row_labels(model)
    return dict(endpoints=[assess_endpoint(model, candidate, labels, with_local_errors)
                           for candidate in candidates],
                row_labels=labels)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--output', required=True)
    parser.add_argument('--comparison', default=COMPARISON)
    parser.add_argument('--endpoint', action='append', help='accepted endpoint id; default: all accepted')
    parser.add_argument('--without-local-errors', action='store_true')
    args = parser.parse_args()
    comparison_bytes = (ROOT/args.comparison).read_bytes()
    comparison = json.loads(comparison_bytes)
    model = joint_fit.GlobalLikelihood(*ORDERS, bin_block_size=BLOCK_SIZE)
    result = assess(model, comparison, args.endpoint, not args.without_local_errors)
    best = min(result['endpoints'], key=lambda record: record['chi2'])
    result.update(
        scope='Goodness of fit at accepted, quadrature-qualified minima of the baseline joint likelihood.',
        quadrature=dict(energy=ORDERS[0], pair_per_segment=ORDERS[1], angle_per_segment=ORDERS[2],
                        bin_block_size=BLOCK_SIZE),
        degrees_of_freedom_counting=dict(
            total='75 residuals (73 observed summaries + 2 normalization constraints) minus 32 free coordinates = 43',
            per_block='nominal p-value uses the block row count as dof (assumption; parameters not subtracted)',
            correction='no correction for bounded or nonlinear parameters'),
        energy_ranges_gev=ENERGY_RANGES,
        chi2_match_tolerance=CHI2_MATCH_TOLERANCE,
        all_chi2_match_comparison=all(record['chi2_matches_comparison'] for record in result['endpoints']),
        best_endpoint=best['id'],
        parameter_names=list(joint_fit.PARAMETER_NAMES),
        input_sha256={args.comparison: hashlib.sha256(comparison_bytes).hexdigest(),
                      **{path: file_hash(path) for path in multi.INPUTS}},
        source_sha256=multi.source_hashes(),
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        command=shlex.join([sys.executable, *sys.argv]))
    multi.write_ledger(ROOT/args.output, result)
    print(args.output)


if __name__ == '__main__':
    main()
