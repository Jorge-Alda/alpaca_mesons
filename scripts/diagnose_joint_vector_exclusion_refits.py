"""Refit without eta-prime radiative observations and report mass-block exclusion.

Run: .venv/bin/python scripts/diagnose_joint_vector_exclusion_refits.py \
  --output data/vector_sector/joint_vector_diagnostic_exclusion_refits.json
Units: GeV and nb; residuals and chi-square are dimensionless. Reuses the
previously run mass-and-leptonic exclusion at identical fit/evaluation orders.
"""
import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import assess_joint_vector_fit_adequacy as adequacy
import fit_joint_vector_multistart as multi
import joint_vector_likelihood_options as options
import refit_joint_vector_sensitivity as refits

ROOT = adequacy.ROOT
EXCLUSION_LEDGER = 'data/vector_sector/joint_vector_block_exclusion.json'
ETA_PRIME_CHANNELS = ('etap_omega_gamma', 'etap_gamma_gamma')


def exclude_eta_prime_radiative_rows(model):
    """Remove the two retained eta-prime radiative observations by channel."""
    indices = [index for index, row in enumerate(model.observations)
               if row['channel'] in ETA_PRIME_CHANNELS]
    channels = [model.observations[index]['channel'] for index in indices]
    if sorted(channels) != sorted(ETA_PRIME_CHANNELS):
        raise ValueError(f'Unexpected eta-prime radiative rows: {channels}')
    model.row_mask[indices] = False
    model.all_rows = False
    return indices


def row_changes(evaluation, baseline, labels):
    """Describe each signed pull and its change from the accepted endpoint."""
    current = np.asarray(evaluation['full_residuals'])
    previous = np.asarray([row['pull'] for row in baseline['standardized_residuals']])
    return [dict(index=index, **label, baseline_pull=float(previous[index]),
                 refit_pull=float(current[index]), pull_change=float(current[index]-previous[index]))
            for index, label in enumerate(labels)]


def summarize_refit(record, model, baseline):
    """Add physical eta-prime width and all row changes to an existing refit."""
    high_order = record['high_order']
    if not high_order['valid']:
        raise ValueError('high-order evaluation failed')
    values = np.asarray(record['refit']['endpoint'])
    width = model.unpack(values)[4]['etap']
    return dict(record=record, eta_prime_total_width_gev=float(width),
                row_changes=row_changes(high_order, baseline, adequacy.row_labels(model)))


def mass_exclusion_record():
    """Reuse the previously recorded mass-and-leptonic exclusion fit."""
    ledger = json.loads((ROOT / EXCLUSION_LEDGER).read_text())
    if ledger['best_accepted_endpoint'] != refits.BEST:
        raise ValueError('Mass exclusion used a different starting endpoint')
    expected_fit = {'energy': 2, 'pair_per_segment': 16, 'angle_per_segment': 16}
    expected_high = {'energy': 4, 'pair_per_segment': 128, 'angle_per_segment': 128}
    for key, expected in [('fit_quadrature', expected_fit), ('evaluation_quadrature', expected_high)]:
        if any(ledger[key][name] != value for name, value in expected.items()):
            raise ValueError(f'Unexpected {key} in mass exclusion ledger')
    return next(entry for entry in ledger['entries'] if entry['name'] == 'without_masses_and_leptonic')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    context = refits.accepted_context()
    baseline_ledger = json.loads((ROOT / refits.ADEQUACY).read_text())
    baseline = next(entry for entry in baseline_ledger['endpoints'] if entry['id'] == refits.BEST)
    low, high = refits.make_models({}, (refits.FIT_ORDERS, refits.HIGH_ORDERS))
    indices = exclude_eta_prime_radiative_rows(low)
    if exclude_eta_prime_radiative_rows(high) != indices:
        raise ValueError('Inconsistent removed row ordering across quadratures')
    item = dict(name='without_eta_prime_radiative_rows', excluded_blocks=(),
                definition='Remove the two retained eta-prime radiative observations.')
    prefix = str(ROOT / (args.output + '.eta_prime'))
    eta_record = refits.run_exclusion(item, context, models=(low, high), checkpoint_prefix=prefix)
    if eta_record['refit']['state'] != 'finished' or not eta_record['high_order']['valid']:
        multi.write_ledger(ROOT / args.output, {'eta_prime_exclusion': eta_record})
        raise RuntimeError('Eta-prime exclusion did not yield a valid high-order evaluation')
    mass_record = mass_exclusion_record()
    result = dict(best_accepted_endpoint=refits.BEST, eta_prime_removed_indices=indices,
                  eta_prime_exclusion=summarize_refit(eta_record, high, baseline),
                  mass_and_leptonic_exclusion=summarize_refit(mass_record, high, baseline),
                  mass_exclusion_source=EXCLUSION_LEDGER,
                  fit_orders=list(refits.FIT_ORDERS), evaluation_orders=list(refits.HIGH_ORDERS),
                  optimizer_settings=refits.OPTIMIZER_SETTINGS,
                  input_sha256={path: adequacy.file_hash(path) for path in
                                (adequacy.COMPARISON, refits.ADEQUACY, EXCLUSION_LEDGER, *multi.INPUTS)},
                  source_sha256=multi.source_hashes(),
                  p_value_convention='Nominal block p-values use retained row count as degrees of freedom; '
                                     'parameters are shared and are not subtracted.')
    multi.write_ledger(ROOT / args.output, result)
    print(args.output, flush=True)


if __name__ == '__main__':
    os.environ.setdefault('OMP_NUM_THREADS', '1')
    main()
