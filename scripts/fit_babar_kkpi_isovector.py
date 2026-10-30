"""Conditional BaBar:2007ceh Table 7 spectrum fits in nb.

Run: .venv/bin/python scripts/fit_babar_kkpi_isovector.py --output /tmp/kkpi.json

Equations (11),(12),(15),(18),(19) of arXiv:0710.4451v3 are used directly.
This is a stable-K* quasi-two-body spectrum model, with no RChT coefficient
interpretation. The 19 available bins are uniformly averaged in energy.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import least_squares
from scipy.stats import chi2


ROOT = Path(__file__).resolve().parents[1]
GEV_MINUS_TWO_TO_NB = 389379.3656
KAON_MASS_GEV = 0.493677
PION_MASS_GEV = 0.13957039
KSTAR_MASS_GEV = 0.89166
BACKGROUND_REFERENCE_GEV = 1.709
PARAMETER_NAMES = ("rho_prime_mass_GeV", "rho_prime_width_GeV",
                   "rho_prime_ee_width_times_branching_eV",
                   "signed_sqrt_background_cross_section_nb",
                   "rho_prime_four_pion_branching_fraction")


def vector_pseudoscalar_phase_space(squared_energy):
    """BaBar Eq. (12), in GeV^3, for a stable K*(892) and kaon."""
    squared_energy = np.asarray(squared_energy)
    discriminant = ((squared_energy-KSTAR_MASS_GEV**2-KAON_MASS_GEV**2)**2
                    -4*KSTAR_MASS_GEV**2*KAON_MASS_GEV**2)
    return np.where(squared_energy > (KSTAR_MASS_GEV+KAON_MASS_GEV)**2,
                    np.maximum(discriminant/squared_energy, 0)**1.5, 0)


def four_pion_phase_space(squared_energy, power):
    """BaBar Eq. (19); its normalization cancels in the running width."""
    return np.maximum(np.asarray(squared_energy)-16*PION_MASS_GEV**2, 0)**power/squared_energy


def cross_section(energies, parameters, four_pion_power=1.5):
    """BaBar Eq. (11), zero resonance phase and real signed background."""
    mass, width, product_eV, background, fraction = parameters
    squared_energy = np.asarray(energies)**2
    phase_space = vector_pseudoscalar_phase_space(squared_energy)
    pole_phase_space = vector_pseudoscalar_phase_space(mass**2)
    running_width = width*(fraction*four_pion_phase_space(squared_energy, four_pion_power)
                          /four_pion_phase_space(mass**2, four_pion_power)+1-fraction)
    reference_squared = BACKGROUND_REFERENCE_GEV**2
    background_coefficient = background*reference_squared/np.sqrt(
        12*np.pi*GEV_MINUS_TWO_TO_NB*vector_pseudoscalar_phase_space(reference_squared))
    resonance = np.sqrt(product_eV*1e-9*width/pole_phase_space)/(
        mass**2-squared_energy-1j*np.sqrt(squared_energy)*running_width)
    return 12*np.pi*GEV_MINUS_TWO_TO_NB*phase_space*np.abs(
        background_coefficient/squared_energy+resonance)**2


def bin_predictions(edges, parameters, four_pion_power=1.5, order=12):
    """Uniform dE averages; the acceptance and ISR response are unavailable."""
    nodes, weights = leggauss(order)
    lower, upper = np.asarray(edges).T
    energies = (lower[:, None]+upper[:, None])/2+(upper-lower)[:, None]*nodes/2
    return cross_section(energies, parameters, four_pion_power)@weights/2


def profiled_residuals(prediction, observations, statistical_errors, relative_errors=()):
    """Profile separate common multiplicative Gaussian nuisances analytically.

    The predicted cross section is (1 + sum(delta_a z_a))*prediction and
    every z_a has an independent unit Gaussian prior. No published bin
    covariance is claimed by this fully correlated approximation.
    """
    statistical_residual = (prediction-observations)/statistical_errors
    relative_errors = np.asarray(relative_errors)
    if not relative_errors.size:
        return statistical_residual, np.empty(0)
    nuisance_columns = prediction[:, None]*relative_errors/statistical_errors[:, None]
    nuisances = -np.linalg.solve(np.eye(relative_errors.size)+nuisance_columns.T@nuisance_columns,
                                nuisance_columns.T@statistical_residual)
    return np.r_[statistical_residual+nuisance_columns@nuisances, nuisances], nuisances


def fit_spectrum(edges, observations, errors, *, background, free_fraction,
                 relative_errors=(), four_pion_power=1.5):
    """Deterministic bounded fits; rank and reached bounds are reported."""
    active = [0, 1, 2]+([3] if background else [])+([4] if free_fraction else [])
    lower = np.array([KSTAR_MASS_GEV+KAON_MASS_GEV+1e-4, 0.03, 1e-6, -5, 0])
    upper = np.array([1.9, 1.5, 5000, 5, 1])
    benchmark = np.array([1.508, 0.418, 129., np.sqrt(.9), .64])
    if not background:
        benchmark[3] = 0

    def expand(values):
        parameters = benchmark.copy()
        parameters[active] = values
        return parameters

    def residual(values):
        prediction = bin_predictions(edges, expand(values), four_pion_power)
        return profiled_residuals(prediction, observations, errors, relative_errors)[0]

    fits = []
    for mass in (1.44, 1.508, 1.65):
        for width in (.2, .418, .8):
            for background_sign in ((-1, 1) if background else (1,)):
                initial = benchmark.copy()
                initial[:2] = mass, width
                initial[3] *= background_sign
                result = least_squares(residual, initial[active], bounds=(lower[active], upper[active]),
                                       x_scale="jac", ftol=1e-11, xtol=1e-11, gtol=1e-10,
                                       max_nfev=3000)
                fits.append(result)
    result = min(fits, key=lambda fit: np.dot(fit.fun, fit.fun))
    parameters = expand(result.x)
    prediction = bin_predictions(edges, parameters, four_pion_power)
    residual_values, nuisances = profiled_residuals(prediction, observations, errors, relative_errors)
    singular_values = np.linalg.svd(result.jac*np.maximum(np.abs(result.x), .01), compute_uv=False)
    rank = int(np.sum(singular_values > singular_values[0]*1e-7))
    covariance = np.linalg.pinv(result.jac.T@result.jac, rcond=1e-12)
    statistic = float(residual_values@residual_values)
    near_bounds = [PARAMETER_NAMES[index] for position, index in enumerate(active)
                   if min((result.x[position]-lower[index])/(upper[index]-lower[index]),
                          (upper[index]-result.x[position])/(upper[index]-lower[index])) < 1e-4]
    return {
        "parameters": dict(zip(PARAMETER_NAMES, map(float, parameters))),
        "conditional_local_errors": dict(zip((PARAMETER_NAMES[index] for index in active),
            [None]*len(active) if rank < len(active) or near_bounds
            else map(float, np.sqrt(np.diag(covariance))))),
        "fixed_parameters": [PARAMETER_NAMES[index] for index in range(5) if index not in active],
        "chi_squared": statistic, "degrees_of_freedom": len(observations)-len(active),
        "degrees_of_freedom_using_local_rank": len(observations)-rank,
        "chi_squared_probability": float(chi2.sf(statistic, len(observations)-len(active))),
        "probability_qualification": "Nominal parameter count; not calibrated for boundaries or unresolved directions",
        "common_relative_systematic_errors": list(relative_errors),
        "common_gaussian_nuisances": list(map(float, nuisances)),
        "four_pion_power": four_pion_power,
        "predictions_nb": prediction.tolist(), "local_jacobian_rank": rank,
        "scaled_jacobian_singular_values": singular_values.tolist(),
        "parameters_at_bounds": near_bounds,
        "optimizer_success": bool(result.success), "optimizer_message": result.message,
        "number_of_deterministic_starts": len(fits),
        "local_errors_qualification": "Conditional linearized errors; unreliable at bounds or for unresolved degeneracies",
    }


def run_fits():
    table = np.genfromtxt(ROOT/"data/vector_sector/babar_2007_isovector_kstar_k_tau_region.csv",
                         delimiter=",", names=True)
    edges = np.column_stack((table["energy_low_GeV"], table["energy_high_GeV"]))
    observed, errors = table["cross_section_nb"], table["statistical_error_nb"]
    results = {}
    for label, relatives in (("statistical", ()), ("dalitz_common", (.15,)),
                             ("dalitz_and_table1_common", (.15, .051))):
        for model, background, free_fraction in (("rho_prime", False, False),
                ("rho_prime_background", True, False), ("rho_prime_background_free_fraction", True, True)):
            results[f"{model}_{label}"] = fit_spectrum(edges, observed, errors, background=background,
                free_fraction=free_fraction, relative_errors=relatives)
    results["rho_prime_background_power_5_over_2_statistical"] = fit_spectrum(
        edges, observed, errors, background=True, free_fraction=False, four_pion_power=2.5)
    benchmark = np.array([1.508, .418, 129., np.sqrt(.9), .64])
    benchmark_prediction = bin_predictions(edges, benchmark)
    refined = bin_predictions(edges, benchmark, order=24)
    return {
        "bibtex_keys": ["BaBar:2007ceh", "ParticleDataGroup:2024cfk"],
        "observable": "BaBar Table 7 isovector K*(892)K component, 19 bins 1.40-1.78 GeV",
        "inputs": {"kaon_mass_GeV": KAON_MASS_GEV, "pion_mass_GeV": PION_MASS_GEV,
                   "kstar_mass_GeV": KSTAR_MASS_GEV, "background_reference_GeV": BACKGROUND_REFERENCE_GEV},
        "assumptions": ["Stable K* phase space from source Eq.12; no finite-K* convolution",
            "Uniform energy-bin averages; ISR luminosity and acceptance response unavailable",
            "Table7 statistical covariance unavailable; statistical errors treated as diagonal",
            "15% Dalitz and 5.1% Table1 common nuisances are correlation assumptions",
            "Real background and zero resonance phase as source Table15",
            "The Table15 benchmark takes the positive square root of its background cross section; its sign is not tabulated",
            "Fixed four-pion fraction 0.64 uses Table15 central value and is conditional, not an independent prior",
            "Only Table7 tau-overlap bins fitted; source Table15 global fit is a benchmark, not reproduced",
            "No RChT coefficient extraction; VPPP dependence cannot be inferred from this component"],
        "parameter_bounds": {"mass_GeV": [KSTAR_MASS_GEV+KAON_MASS_GEV+1e-4, 1.9],
            "width_GeV": [.03, 1.5], "ee_width_times_branching_eV": [1e-6, 5000],
            "signed_sqrt_background_nb": [-5, 5], "four_pion_fraction": [0, 1]},
        "table15_benchmark": {"parameters": dict(zip(PARAMETER_NAMES, benchmark.tolist())),
            "predictions_nb": benchmark_prediction.tolist(),
            "statistical_chi_squared": float(np.sum(((benchmark_prediction-observed)/errors)**2))},
        "bin_quadrature_refinement_max_absolute_nb": float(np.max(np.abs(refined-benchmark_prediction))),
        "bin_edges_GeV": edges.tolist(), "observed_nb": observed.tolist(),
        "statistical_errors_nb": errors.tolist(), "fits": results,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    content = json.dumps(run_fits(), indent=2, allow_nan=False)+"\n"
    if arguments.output:
        arguments.output.write_text(content)
        print(arguments.output)
    else:
        print(content, end="")
