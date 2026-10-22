"""Profile charged photon PP production widths; all masses and widths use GeV.

Run: .venv/bin/python -m scripts.profile_photon_pp_gauss_legendre --repeats 7
The deterministic source is coupling_inputs()[0] from the quadrature benchmark.
Timing includes amplitude preparation; instrumentation runs separately.
"""
import argparse
import cProfile
import importlib
import inspect
import json
import pstats
import statistics
import time
from unittest.mock import patch

import numpy as np

from scripts.benchmark_photon_pp_quadrature import coupling_inputs

photon_module = importlib.import_module('alpaca_chiral.photon_two_pseudoscalar')


def profile_channel(first, second, repeats, order, epsrel):
    """Return warmed timings, vectorized work counts, and Python hotspots."""
    arguments = (1.5, 100., coupling_inputs()[0], first, second)
    controls = dict(order=order, epsrel=epsrel)
    width_function = photon_module.photon_two_pseudoscalar_width
    width_function(*arguments, **controls)
    durations = []
    for _ in range(repeats):
        start = time.perf_counter()
        width = width_function(*arguments, **controls)
        durations.append(time.perf_counter()-start)

    amplitude = photon_module.photon_two_pseudoscalar_amplitude(*arguments)
    amplitude_type = type(amplitude)
    original_call = amplitude_type.__call__
    original_quad = photon_module.quad
    work = dict(vectorized_passes=[], adaptive_quad_calls=0,
                adaptive_integrand_calls=0)

    def counted_amplitude(self, s, t):
        """Count allocated Dalitz grid points, including zero-width panels."""
        if np.ndim(t):
            work['vectorized_passes'].append(dict(
                shape=list(np.shape(t)), grid_points=int(np.size(t))))
        return original_call(self, s, t)

    def counted_quad(function, *args, **kwargs):
        """Count fallback integrals and their integrand calls independently."""
        work['adaptive_quad_calls'] += 1

        def counted_integrand(coordinate):
            work['adaptive_integrand_calls'] += 1
            return function(coordinate)

        return original_quad(counted_integrand, *args, **kwargs)

    with patch.object(amplitude_type, '__call__', counted_amplitude), \
            patch.object(photon_module, 'quad', counted_quad):
        counted_width = width_function(*arguments, **controls)
    work['grid_points_total'] = sum(
        entry['grid_points'] for entry in work['vectorized_passes'])

    profiler = cProfile.Profile()
    profiled_width = profiler.runcall(width_function, *arguments, **controls)
    statistics_profile = pstats.Stats(profiler)
    hotspots = []
    for (filename, line, function), values in sorted(
            statistics_profile.stats.items(), key=lambda item: item[1][3],
            reverse=True)[:15]:
        primitive_calls, calls, own_seconds, cumulative_seconds, _ = values
        hotspots.append(dict(function=function, file=filename, line=line,
                             calls=calls, primitive_calls=primitive_calls,
                             own_seconds=own_seconds,
                             cumulative_seconds=cumulative_seconds))
    poles = [dict(invariant=('s', 't', 'u')[invariant], mass_GeV=mass,
                  width_GeV=pole_width,
                  real_pole_in_s_domain=bool(amplitude.threshold_squared < mass**2 < 1.5**2)
                  if invariant == 0 else None)
             for invariant, mass, pole_width, _, _ in amplitude.terms]
    return dict(pair=[first, second], width_GeV=width,
                counted_width_GeV=counted_width, profiled_width_GeV=profiled_width,
                median_seconds=statistics.median(durations),
                seconds=durations, pole_terms=poles, work=work,
                profiled_total_calls=statistics_profile.total_calls,
                hotspots=hotspots)


def benchmark(repeats=7, order=16, epsrel=1e-5):
    """Compare the two charged channels at identical input masses and sources."""
    records = [profile_channel(first, second, repeats, order, epsrel)
               for first, second in [('pi+', 'pi-'), ('K+', 'K-')]]
    signature = inspect.signature(photon_module.photon_two_pseudoscalar_width)
    return dict(ma_GeV=1.5, fa_GeV=100.,
                source='scripts.benchmark_photon_pp_quadrature.coupling_inputs()[0]',
                repeats=repeats, order=order, epsrel=epsrel,
                production_defaults=dict(order=signature.parameters['order'].default,
                                         epsrel=signature.parameters['epsrel'].default),
                kaon_to_pion_time_ratio=records[1]['median_seconds']/records[0]['median_seconds'],
                results=records)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repeats', type=int, default=7)
    parser.add_argument('--order', type=int, default=16)
    parser.add_argument('--epsrel', type=float, default=1e-5)
    options = parser.parse_args()
    print(json.dumps(benchmark(options.repeats, options.order, options.epsrel), indent=2))
