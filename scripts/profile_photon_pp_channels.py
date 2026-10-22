"""Profile the adaptive reference widths for charged pions and kaons.

Run from the repository root::

    .venv/bin/python -m scripts.profile_photon_pp_channels --repeats 3

Uses m_a=1.5 GeV, f_a=100 GeV, dimensionless no-mixing axial source
diag(0.8, 0, 0), and a zero pseudoscalar source. Fit inputs are loaded once
and shared. Existing resonance and integration prescriptions are retained.
Output consists of one JSON record per channel. Wall times exclude import,
fit loading, warm-up, and the separately instrumented cProfile evaluation.
"""

import argparse
import cProfile
import json
import pstats
from statistics import median
from time import perf_counter

from scripts.implementation.alp_photon_pp_width import (
    alp_photon_pp_width, build_photon_pp_amplitude,
)
from scripts.implementation.alp_vv_width import read_preliminary_fit_parameters


def profile_channel(final_state, fitted_parameters, repeats, hotspot_count=12):
    """Return timing and cumulative profiling statistics for one width."""
    def calculate_width():
        return alp_photon_pp_width(
            1.5, [0.8, 0.0, 0.0], [0.0, 0.0, 0.0], final_state, 100.0,
            fit_parameters=fitted_parameters,
        )

    calculate_width()
    amplitude = build_photon_pp_amplitude(
        1.5, [0.8, 0.0, 0.0], [0.0, 0.0, 0.0], final_state, 100.0,
        fit_parameters=fitted_parameters,
    )
    elapsed_seconds = []
    widths = []
    for _ in range(repeats):
        start = perf_counter()
        widths.append(calculate_width())
        elapsed_seconds.append(perf_counter() - start)
    profiler = cProfile.Profile()
    profiled_width = profiler.runcall(calculate_width)
    statistics = pstats.Stats(profiler)
    hotspots = []
    for (filename, line_number, function_name), values in sorted(
        statistics.stats.items(), key=lambda item: item[1][3], reverse=True
    )[:hotspot_count]:
        primitive_calls, total_calls, own_seconds, cumulative_seconds, _ = values
        hotspots.append({
            "function": f"{filename}:{line_number}({function_name})",
            "primitive_calls": primitive_calls,
            "total_calls": total_calls,
            "own_seconds": own_seconds,
            "cumulative_seconds": cumulative_seconds,
        })
    return {
        "final_state": ["gamma", *final_state],
        "alp_mass_GeV": 1.5,
        "f_a_GeV": 100.0,
        "axial_diagonal": [0.8, 0.0, 0.0],
        "pseudoscalar_diagonal": [0.0, 0.0, 0.0],
        "epsrel": 1e-5,
        "epsabs": 0.0,
        "limit": 150,
        "poles": [
            {"invariant": ("s", "t", "u")[index], "mass_GeV": mass,
             "width_GeV": width,
             "constant": [constant.real, constant.imag],
             "linear": [linear.real, linear.imag]}
            for index, mass, width, constant, linear in amplitude.terms
        ],
        "repeats": repeats,
        "median_wall_seconds": median(elapsed_seconds),
        "wall_seconds": elapsed_seconds,
        "widths_GeV": widths,
        "profiled_width_GeV": profiled_width,
        "profile_total_calls": statistics.total_calls,
        "profile_total_seconds": statistics.total_tt,
        "cumulative_hotspots": hotspots,
    }


def main():
    """Load one fit and print the pion and kaon profiling records."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repeats", type=int, default=3)
    arguments = parser.parse_args()
    if arguments.repeats < 1:
        parser.error("--repeats must be at least 1")
    fitted_parameters = read_preliminary_fit_parameters()
    for final_state in (("pi+", "pi-"), ("K+", "K-")):
        print(json.dumps(profile_channel(final_state, fitted_parameters, arguments.repeats)))


if __name__ == "__main__":
    main()
