"""Per-evaluation cost of the isospin-breaking VPP option against the baseline likelihood.

Timing and call counting only: no fit, no optimizer call, no change of the
likelihood numerics. The likelihood code is instrumented from outside by
temporary counting/timing wrappers that call the original functions.

Run:
  OMP_NUM_THREADS=1 .venv/bin/python scripts/benchmark_isospin_breaking_vpp_cost.py \
      --output data/vector_sector/isospin_breaking_vpp_cost_benchmark.json

Point: the endpoint of run 2 of entry 'scan_point' of
data/vector_sector/joint_vector_isospin_breaking_vpp_fit.json (selected by
fit_joint_vector_isospin_breaking_vpp.restart_run: lowest (4,128,128) chi2).
Quadrature (2,16,16) with four-bin blocks (the fit order of the driver).

Configurations:
  baseline: joint_fit.GlobalLikelihood with the first 32 endpoint coordinates;
  option_r_zero: IsospinBreakingVppLikelihood (option on, default switches,
      the class used for the original scan_point fit) with the 32 coordinates
      and r = (0, 0, 0);
  option_restart: the same class at the full 35-coordinate endpoint;
  reuse_r_zero, reuse_restart: the driver's FitLikelihood (subclass of
      IsospinBreakingVppLikelihoodWithReuse, bitwise-equal outputs) at the
      same two points.

Timing method (memory/rcht/odd-parity/kkpi-global-fit-channel-jacobian.md):
one warm-up residual and Jacobian per configuration, then REPETITIONS rounds in
which the three configurations alternate; medians of perf_counter wall times.
Residual rounds and Jacobian rounds are separate loops.

Call counts (one residual, one Jacobian per configuration; separate pass with
wrappers so the timings above carry no wrapper overhead):
  full residuals (GlobalLikelihood.residual), split into those made inside the
  inherited GlobalLikelihood.jacobian and those made by the option's jacobian
  override; channel-restricted BaBar evaluations (babar_blocks with one
  channel); full-channel and derivative babar_blocks calls; per-block
  integrations (_babar_blocks_full); old_blocks; running-propagator calls
  (split into the zero-excitation calls added by the option's mode_scaled and
  the others); reference_propagator_pair; reference propagator (pole
  function) calls; s_channel_light_correction; vertex-factor evaluations.
Accumulated wall time inside each wrapper is reported as a share estimate
(it includes the wrapper overhead and nested calls; labelled as approximate).
"""
import os

os.environ.setdefault('OMP_NUM_THREADS', '1')

import argparse  # noqa: E402
import contextlib  # noqa: E402
import json  # noqa: E402
import platform  # noqa: E402
import statistics  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from collections import Counter, defaultdict  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import scipy  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fit_babar_kkpi_global_vector_sector as joint_fit  # noqa: E402
import fit_joint_vector_isospin_breaking_vpp as driver  # noqa: E402
import joint_vector_isospin_breaking_vpp_option as vpp_option  # noqa: E402
import joint_vector_isospin_breaking_vpp_reuse as reuse_module  # noqa: E402
import joint_vector_likelihood_options as options_module  # noqa: E402

LEDGER = 'data/vector_sector/joint_vector_isospin_breaking_vpp_fit.json'
ORDERS = driver.FIT_ORDERS
REPETITIONS = 5
CONFIGURATIONS = ('baseline', 'option_r_zero', 'option_restart', 'reuse_r_zero', 'reuse_restart')


def restart_point(ledger_path=LEDGER):
    """(run index, 35-coordinate endpoint) of the restart run of the ledger."""
    ledger = json.loads((driver.ROOT/ledger_path).read_bytes())
    index, run = driver.restart_run(ledger)
    return index, np.asarray(run['endpoint'], dtype=float)


def configurations(endpoint, orders=ORDERS):
    """{name: (model, point)} of the three compared evaluations."""
    baseline = joint_fit.GlobalLikelihood(*orders, bin_block_size=driver.BLOCK_SIZE)
    option = vpp_option.IsospinBreakingVppLikelihood(*orders, bin_block_size=driver.BLOCK_SIZE,
                                                     **driver.MODEL_OPTIONS)
    reuse = driver.FitLikelihood(*orders, bin_block_size=driver.BLOCK_SIZE, **driver.MODEL_OPTIONS)
    base_point = endpoint[:driver.BASELINE_COUNT].copy()
    return dict(baseline=(baseline, base_point),
                option_r_zero=(option, option.extend(base_point)),
                option_restart=(option, endpoint.copy()),
                reuse_r_zero=(reuse, reuse.extend(base_point)),
                reuse_restart=(reuse, endpoint.copy()))


class Instrument:
    """Counting and timing wrappers around the likelihood functions (patched on enter)."""

    def __init__(self):
        self.counts = Counter()
        self.seconds = defaultdict(float)
        self.jacobian_depth = 0
        self.patches = []

    def timed(self, label, function, classify=None):
        """Wrap function: count calls (label or classify(*args)) and accumulate wall time."""
        def wrapper(*args, **kwargs):
            key = classify(*args, **kwargs) if classify else label
            self.counts[key] += 1
            started = time.perf_counter()
            try:
                return function(*args, **kwargs)
            finally:
                self.seconds[key] += time.perf_counter()-started
        return wrapper

    def patch(self, owner, name, replacement):
        self.patches.append((owner, name, getattr(owner, name)))
        setattr(owner, name, replacement)

    def __enter__(self):
        likelihood = joint_fit.GlobalLikelihood
        original_residual = likelihood.residual
        original_jacobian = likelihood.jacobian

        def residual_label(*args, **kwargs):
            return 'full_residual_inside_baseline_jacobian' if self.jacobian_depth else 'full_residual_outside'

        def jacobian(model, *args, **kwargs):
            self.jacobian_depth += 1
            try:
                return original_jacobian(model, *args, **kwargs)
            finally:
                self.jacobian_depth -= 1

        def babar_label(model, values, state, derivatives=False, channels=joint_fit.CHANNELS):
            if derivatives:
                return 'babar_blocks_derivative_pass'
            return 'babar_blocks_all_channels' if tuple(channels) == tuple(joint_fit.CHANNELS) \
                else 'babar_blocks_channel_restricted'

        def running_label(name, invariant, beta, *args):
            ground = isinstance(beta, tuple) and beta == (0., 0.)
            return f'running_propagator_{name}' + ('_zero_excitation' if ground else '')

        self.patch(likelihood, 'residual', self.timed(None, original_residual, residual_label))
        self.patch(likelihood, 'jacobian', jacobian)
        self.patch(likelihood, 'babar_blocks', self.timed(None, likelihood.babar_blocks, babar_label))
        self.patch(likelihood, '_babar_blocks_full', self.timed('block_integration', likelihood._babar_blocks_full))
        self.patch(likelihood, 'old_blocks', self.timed('old_blocks', likelihood.old_blocks))
        running = self.timed(None, vpp_option.RUNNING_PROPAGATOR, running_label)
        pair = self.timed('reference_propagator_pair', vpp_option.REFERENCE_PROPAGATOR_PAIR)
        pole = self.timed('reference_propagator', joint_fit.reference_propagator)
        for owner, name, replacement in ((joint_fit, 'running_propagator', running),
                                         (vpp_option, 'RUNNING_PROPAGATOR', running),
                                         (joint_fit, 'reference_propagator_pair', pair),
                                         (vpp_option, 'REFERENCE_PROPAGATOR_PAIR', pair),
                                         (joint_fit, 'reference_propagator', pole),
                                         (options_module, 'REFERENCE_PROPAGATOR', pole)):
            self.patch(owner, name, replacement)
        def memory_label(memory, name, invariant, beta, *args):
            return 'reuse_' + running_label(name, invariant, beta)

        memory = reuse_module.GroundPropagatorMemory
        self.patch(memory, '__call__', self.timed(None, memory.__call__, memory_label))
        self.patch(reuse_module, 'S_CHANNEL_FACTORS', self.timed('s_channel_factors', reuse_module.S_CHANNEL_FACTORS))
        option = vpp_option.IsospinBreakingVppLikelihood
        self.patch(option, 's_channel_light_correction',
                   self.timed('s_channel_light_correction', option.s_channel_light_correction))
        self.patch(vpp_option, 'vertex_factors', self.timed('vertex_factors', vpp_option.vertex_factors))
        self.patch(vpp_option, 's_channel_factors', self.timed('s_channel_factors', vpp_option.s_channel_factors))
        return self

    def __exit__(self, *exception):
        for owner, name, original in reversed(self.patches):
            setattr(owner, name, original)
        self.patches = []


def count_calls(model, point, method):
    """Counts and wrapper seconds of one residual or Jacobian evaluation."""
    with Instrument() as instrument:
        getattr(model, method)(point)
    return dict(sorted(instrument.counts.items())), {key: round(value, 6)
                                                     for key, value in sorted(instrument.seconds.items())}


def time_evaluations(cases, repetitions=REPETITIONS):
    """Warm-up, then alternating median wall times of residual and Jacobian."""
    for model, point in cases.values():
        model.residual(point)
        model.jacobian(point)
    samples = {method: {name: [] for name in cases} for method in ('residual', 'jacobian')}
    for method in ('residual', 'jacobian'):
        for _ in range(repetitions):
            for name, (model, point) in cases.items():
                started = time.perf_counter()
                getattr(model, method)(point)
                samples[method][name].append(time.perf_counter()-started)
    return {method: {name: dict(median_seconds=statistics.median(values), samples_seconds=values)
                     for name, values in by_name.items()} for method, by_name in samples.items()}


def outputs_unchanged_by_instrument(cases):
    """Residual and Jacobian bitwise equal with and without the wrappers (wrapper correctness)."""
    result = {}
    for name, (model, point) in cases.items():
        plain = (model.residual(point), model.jacobian(point))
        with Instrument():
            wrapped = (model.residual(point), model.jacobian(point))
        result[name] = bool(all(np.array_equal(a, b) for a, b in zip(plain, wrapped)))
    return result


def reuse_bitwise_equal(cases):
    """Bit patterns of residual and Jacobian: reuse class against the option class at both points."""
    result = {}
    for point_name in ('r_zero', 'restart'):
        option, point = cases[f'option_{point_name}']
        reuse, reuse_point = cases[f'reuse_{point_name}']
        result[point_name] = {method: bool(np.array_equal(getattr(option, method)(point).view(np.uint64),
                                                          getattr(reuse, method)(reuse_point).view(np.uint64)))
                              for method in ('residual', 'jacobian')}
    return result


def wall_time_estimate(timings, nfev, njev, configuration='option_restart'):
    """nfev residuals + njev Jacobians at the measured median costs (optimizer linear algebra neglected)."""
    residual = timings['residual'][configuration]['median_seconds']
    jacobian = timings['jacobian'][configuration]['median_seconds']
    return dict(nfev=nfev, njev=njev, seconds=nfev*residual+njev*jacobian,
                assumption='each optimizer iteration costs one Jacobian plus the residual evaluations it '
                           'counts in nfev; trf linear algebra, callbacks and checkpoint writes neglected; '
                           'costs at the restart point taken as representative of the whole run')


def run(output=None, repetitions=REPETITIONS):
    index, endpoint = restart_point()
    cases = configurations(endpoint)
    counts = {name: {method: dict(zip(('calls', 'wrapper_seconds'), count_calls(model, point, method)))
                     for method in ('residual', 'jacobian')} for name, (model, point) in cases.items()}
    timings = time_evaluations(cases, repetitions)
    base = timings['jacobian']['baseline']['median_seconds']
    ledger = json.loads((driver.ROOT/LEDGER).read_bytes())
    scan_runs = ledger['entries'][0]['runs']
    result = dict(
        point=dict(ledger=LEDGER, entry='scan_point', run_index=index, coefficients_gev_minus2=endpoint[32:].tolist()),
        orders=list(ORDERS), bin_block_size=driver.BLOCK_SIZE, repetitions=repetitions,
        omp_num_threads=os.environ.get('OMP_NUM_THREADS'),
        machine=dict(platform=platform.platform(), processor=platform.processor(), python=platform.python_version(),
                     numpy=np.__version__, scipy=scipy.__version__, cpu_count=os.cpu_count()),
        timings=timings,
        jacobian_ratio_to_baseline={name: timings['jacobian'][name]['median_seconds']/base for name in cases},
        residual_ratio_to_baseline={name: timings['residual'][name]['median_seconds']
                                    / timings['residual']['baseline']['median_seconds'] for name in cases},
        counts=counts,
        instrument_preserves_outputs=outputs_unchanged_by_instrument(cases),
        reuse_bitwise_equal_to_option=reuse_bitwise_equal(cases),
        recorded_scan_point_runs=[dict(nfev=run_['nfev'], njev=run_['njev'], refit_seconds=run_['refit_seconds'],
                                       predicted_seconds=wall_time_estimate(timings, run_['nfev'],
                                                                            run_['njev'])['seconds'])
                                  for run_ in scan_runs],
        restart_estimate_at_max_nfev_5000={
            name: wall_time_estimate(timings, 5000, scan_runs[index]['njev'], name)
            for name in ('option_restart', 'reuse_restart')})
    if output:
        (driver.ROOT/output).write_text(json.dumps(result, indent=1)+'\n')
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--output', default=None)
    parser.add_argument('--repetitions', type=int, default=REPETITIONS)
    args = parser.parse_args(argv)
    print(json.dumps(run(args.output, args.repetitions), indent=1))


if __name__ == '__main__':
    main()
