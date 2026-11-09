"""Linear chi^2 screen of the 36 fixed-pole shifts at the best accepted endpoint.

Run: .venv/bin/python scripts/screen_joint_vector_pole_shifts.py --output data/vector_sector/joint_vector_pole_shift_screen.json

For each shift the likelihood is evaluated at the fixed best accepted endpoint
(basin_probes/even_shared_widths) at the fit quadrature (2,16,16), four-bin
blocks. With r0 the baseline residual, r the shifted residual and J the
baseline Jacobian (all 32 coordinates, forward differences of
GlobalLikelihood.jacobian), the screen reports
  fixed-parameter Delta chi^2 = |r|^2 - |r0|^2,
  linear refit Delta chi^2   = min_d |r + J d|^2 - min_d |r0 + J d|^2,
and the linear parameter shift d(r) - d(r0) in baseline local errors.
APPROXIMATION: linearization in the parameters and in the pole shift, no
parameter bounds. The screen is used only to order the refits and as a
cross-check; it never removes a variation.

Units: GeV for pole masses and widths; chi^2 dimensionless.
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import assess_joint_vector_fit_adequacy as adequacy
import fit_babar_kkpi_global_vector_sector as joint_fit
import fit_joint_vector_multistart as multi
import joint_vector_likelihood_options as options_module
import refit_joint_vector_sensitivity as refit_module

ROOT = joint_fit.ROOT


def linear_minimum(residual, scaled_jacobian):
    """Least-squares step d (scaled coordinates) and min |r + J d|^2."""
    step = np.linalg.lstsq(scaled_jacobian, -residual, rcond=None)[0]
    remainder = residual+scaled_jacobian@step
    return step, float(remainder @ remainder)


def screen(context, orders=refit_module.FIT_ORDERS):
    """Fixed-parameter and linear-refit Delta chi^2 for every pole shift."""
    endpoint = np.asarray(context['accepted'][0]['endpoint'])
    scales = np.asarray(context['scales'])
    local_errors = np.asarray(context['local_errors'])
    baseline_model = options_module.VariedLikelihood(*orders, bin_block_size=refit_module.BLOCK_SIZE)
    baseline_residual = baseline_model.residual(endpoint)
    scaled_jacobian = baseline_model.jacobian(endpoint)*scales
    baseline_step, baseline_linear = linear_minimum(baseline_residual, scaled_jacobian)
    baseline_chi2 = float(baseline_residual @ baseline_residual)
    entries = []
    for variation in options_module.pole_shift_variations():
        model = options_module.VariedLikelihood(*orders, bin_block_size=refit_module.BLOCK_SIZE,
                                                **variation['options'])
        residual = model.residual(endpoint)
        step, linear = linear_minimum(residual, scaled_jacobian)
        shift = (step-baseline_step)*scales/local_errors
        largest = int(np.argmax(abs(shift)))
        entries.append(dict(name=variation['name'], pole_shift=variation['options']['pole_shift'],
                            fixed_parameter_delta_chi2=float(residual @ residual)-baseline_chi2,
                            linear_refit_delta_chi2=linear-baseline_linear,
                            linear_maximum_shift_in_local_errors=float(abs(shift[largest])),
                            linear_maximum_shift_parameter=joint_fit.PARAMETER_NAMES[largest]))
    order = sorted(entries, key=lambda entry: -abs(entry['linear_refit_delta_chi2']))
    return dict(baseline_chi2_fit_order=baseline_chi2, baseline_linear_minimum=baseline_linear,
                entries=entries, refit_order=[entry['name'] for entry in order])


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    started = time.perf_counter()
    context = refit_module.accepted_context()
    ledger = dict(endpoint=refit_module.BEST,
                  quadrature=dict(zip(('energy', 'pair_per_segment', 'angle_per_segment'), refit_module.FIT_ORDERS),
                                  bin_block_size=refit_module.BLOCK_SIZE),
                  approximation='linear in parameters and pole shift; bounds ignored; ordering only',
                  **screen(context),
                  input_sha256={path: adequacy.file_hash(path) for path in (options_module.PDG_POLES, *multi.INPUTS)},
                  source_sha256=multi.source_hashes(), elapsed_seconds=time.perf_counter()-started)
    multi.write_ledger(ROOT/args.output, ledger)
    print(args.output)


if __name__ == '__main__':
    main()
