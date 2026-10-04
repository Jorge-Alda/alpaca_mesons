"""Validated deterministic loaders for the extracted odd-parity data release.

Only the PDG omega branching fraction is connected to a predictor.  The other
loaders intentionally return provenance-preserving tables without constructing
observations.

Example:
    .venv/bin/python -m unittest scripts.implementation.test_odd_parity_data
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
from types import MappingProxyType
from typing import Any, Mapping

import numpy as np

from .correlated_fit import FitProblem, Observation, OPEConstraint, Parameter
from .odd_parity_rcht import OPEParameters, OPE_LABELS, fixed_ope_prior, ope_linear_system
from .radiative_odd_parity import (
    OmegaPiGammaInputs,
    omega_pi_gamma_branching_fraction_predictor,
)


DATA_DIRECTORY = Path(__file__).resolve().parents[2] / "data" / "odd_parity"
CANONICAL_COUPLING_ORDER = tuple(f"c{i}" for i in range(1, 8)) + tuple(f"d{i}" for i in range(1, 5))
LIKELIHOOD_SCHEMA_VERSION = "odd-parity-likelihood-v1"

FILE_SHA256 = MappingProxyType({
    "babar_3pi_born_spectrum.csv": "b30c8ddb8506196423ac9136305856347fd9cc9605d2df92464dad260e750af7",
    "babar_3pi_born_spectrum.metadata.json": "0f5d48d38aa3ead8d3952c94e684ae78a62b5debf1f5cf983a7c86c26aedd5b5",
    "babar_pi0_table11.csv": "ab85df6584eafc462a52714c4e0abbd714bff951b9322a2056b3012ccc032785",
    "babar_pi0_table11.metadata.json": "49ea47899bf3ab9c178603a07384e7d59569601290c5e2c1f938939b6a13ac96",
    "cleo_tff.csv": "d382e37a13820d9031fb66b96ccf1893a22f9c4da23e1c7d68cf343bb0b02d6b",
    "cleo_tff.metadata.json": "3a1894956729889df3573b51b5ec07dd1ecd09091d242aef7abfc69b85d618ff",
    "kloe_phi_3pi_dalitz_fit.csv": "0021247879b7fdee7acbd475e6cd28de153c055750e04641de668678a6897291",
    "kloe_phi_3pi_dalitz_fit.metadata.json": "9b287aa9ecd7fb8d0801bb64aad0910198cc9d065b039d3f9c5fc1fe25b197ae",
    "omega_pi0_gamma_pdg.csv": "33011a1a3ca94294fcb8fd22f478671cd1d7e486db647e6fa66bdf06c2eb79a9",
    "omega_pi0_gamma_pdg.metadata.json": "97762ffa5642bbd1b74a418aba92e3d2d5858b076ef67698ef6a6459863b5aba",
})


@dataclass(frozen=True)
class OddParityTable:
    rows: np.ndarray
    metadata: Mapping[str, Any]


def _freeze(value: Any) -> Any:
    if isinstance(value, dict):
        return MappingProxyType({key: _freeze(item) for key, item in value.items()})
    if isinstance(value, list):
        return tuple(_freeze(item) for item in value)
    return value


def _verified_bytes(path: Path) -> bytes:
    if path.name not in FILE_SHA256:
        raise ValueError(f"no recorded SHA256 for {path.name}")
    try:
        content = path.read_bytes()
    except OSError as error:
        raise ValueError(f"cannot read required data file {path.name}") from error
    digest = hashlib.sha256(content).hexdigest()
    if digest != FILE_SHA256[path.name]:
        raise ValueError(f"SHA256 mismatch for {path.name}")
    return content


def _metadata(
    stem: str,
    directory: Path,
    expected_units: Any,
    response_column_units: Mapping[str, str],
) -> Mapping[str, Any]:
    content = _verified_bytes(directory / f"{stem}.metadata.json")
    try:
        metadata = json.loads(content, parse_constant=lambda value: (_ for _ in ()).throw(ValueError(value)))
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as error:
        raise ValueError(f"invalid metadata JSON for {stem}") from error
    if not isinstance(metadata, dict) or metadata.get("units") != expected_units:
        raise ValueError(f"metadata units mismatch for {stem}")
    likelihood = metadata.get("likelihood")
    required = {
        "schema_version", "family", "authorization", "asymmetric_error_rule",
        "diagonal_source_ids", "nuisance_sources", "nuisance_source_ids_excluded_from_diagonal",
    }
    if not isinstance(likelihood, dict) or set(likelihood) != required:
        raise ValueError(f"invalid likelihood metadata schema for {stem}")
    if likelihood["schema_version"] != LIKELIHOOD_SCHEMA_VERSION:
        raise ValueError(f"unsupported likelihood metadata version for {stem}")
    if likelihood["family"] != "diagonal_gaussian" or likelihood["asymmetric_error_rule"] != "mean_absolute":
        raise ValueError(f"invalid likelihood convention for {stem}")
    diagonal = likelihood["diagonal_source_ids"]
    excluded = likelihood["nuisance_source_ids_excluded_from_diagonal"]
    nuisances = likelihood["nuisance_sources"]
    if not all(isinstance(value, list) for value in (diagonal, excluded, nuisances)):
        raise ValueError(f"invalid likelihood source lists for {stem}")
    nuisance_ids = [item.get("source_id") for item in nuisances if isinstance(item, dict)]
    valid_diagonal = all(isinstance(item, str) and item for item in diagonal) and len(set(diagonal)) == len(diagonal)
    valid_nuisance_ids = (
        len(nuisance_ids) == len(nuisances)
        and all(isinstance(item, str) and item for item in nuisance_ids)
        and len(set(nuisance_ids)) == len(nuisance_ids)
    )
    valid_excluded = all(isinstance(item, str) and item for item in excluded)
    if (not valid_diagonal or not valid_nuisance_ids or not valid_excluded
            or excluded != nuisance_ids or set(diagonal) & set(excluded)):
        raise ValueError(f"inconsistent diagonal and nuisance sources for {stem}")
    for item in nuisances:
        if item.get("predictor_connection") != "absent":
            raise ValueError(f"invalid unconnected nuisance metadata for {stem}")
        if item.get("response") == "multiplicative":
            if set(item) != {"source_id", "fraction", "response", "predictor_connection"}:
                raise ValueError(f"invalid multiplicative nuisance metadata for {stem}")
            fraction = item["fraction"]
            if (not isinstance(fraction, (int, float)) or isinstance(fraction, bool)
                    or not math.isfinite(float(fraction)) or float(fraction) < 0.0):
                raise ValueError(f"invalid nuisance fraction for {stem}")
        elif item.get("response") == "additive":
            if set(item) != {"source_id", "response", "response_column", "unit", "predictor_connection"}:
                raise ValueError(f"invalid additive nuisance metadata for {stem}")
            column, unit = item["response_column"], item["unit"]
            if not isinstance(column, str) or response_column_units.get(column) != unit:
                raise ValueError(f"invalid additive nuisance response column or unit for {stem}")
        else:
            raise ValueError(f"invalid nuisance response type for {stem}")
    return _freeze(metadata)


def _csv_rows(stem: str, directory: Path, expected_header: tuple[str, ...]) -> list[dict[str, str]]:
    content = _verified_bytes(directory / f"{stem}.csv")
    try:
        reader = csv.DictReader(content.decode("utf-8").splitlines())
        if tuple(reader.fieldnames or ()) != expected_header:
            raise ValueError(f"CSV schema mismatch for {stem}")
        rows = list(reader)
    except UnicodeDecodeError as error:
        raise ValueError(f"invalid UTF-8 CSV for {stem}") from error
    if not rows or any(None in row for row in rows):
        raise ValueError(f"empty or malformed CSV for {stem}")
    return rows


def _finite(text: str, field: str, *, positive: bool = False) -> float:
    try:
        value = float(text)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{field} must be finite") from error
    if not math.isfinite(value) or (positive and value <= 0.0):
        raise ValueError(f"{field} must be {'positive and ' if positive else ''}finite")
    return value


def _table(stem: str, header: tuple[str, ...], dtype, converters, units, directory, response_column_units=None) -> OddParityTable:
    directory = Path(directory)
    raw = _csv_rows(stem, directory, header)
    converted = []
    for row in raw:
        converted.append(tuple(converter(row[name], name) for name, converter in zip(header, converters)))
    rows = np.array(converted, dtype=dtype)
    rows.setflags(write=False)
    return OddParityTable(rows, _metadata(stem, directory, units, response_column_units or {}))


_text = lambda value, name: value if value != "" else (_ for _ in ()).throw(ValueError(f"{name} must not be empty"))
_number = lambda value, name: _finite(value, name)
_positive = lambda value, name: _finite(value, name, positive=True)


def load_cleo_tff(directory: Path = DATA_DIRECTORY) -> OddParityTable:
    header = ("meson", "source_table", "q2_GeV2", "q2F_GeV", "stat_GeV", "syst_GeV")
    return _table("cleo_tff", header, [(header[0], "U16"), (header[1], "U24")] + [(name, "f8") for name in header[2:]],
                  (_text, _text, _positive, _number, _positive, _positive),
                  {"q2": "GeV^2", "q2F": "GeV"}, directory)


def load_babar_pi0_table11(directory: Path = DATA_DIRECTORY) -> OddParityTable:
    header = ("q2_interval_GeV2", "q2bar_GeV2", "dsigma_stat_fb_per_GeV2", "dsigma_syst_fb_per_GeV2",
              "q2F_MeV", "q2F_err_low_MeV", "q2F_err_high_MeV")
    return _table("babar_pi0_table11", header, [(header[0], "U16")] + [(name, "f8") for name in header[1:]],
                   (_text, _positive, _positive, _positive, _positive, _positive, _positive),
                  {"q2": "GeV^2", "cross_section": "fb/GeV^2", "q2F": "MeV"}, directory)


def load_babar_3pi_born_spectrum(directory: Path = DATA_DIRECTORY) -> OddParityTable:
    header = ("sqrt_s_GeV", "Ncorr", "Ncorr_stat", "Ncorr_syst", "efficiency_percent", "efficiency_err_percent",
              "dL_nb_inv", "sigma_nb", "sigma_stat_nb", "sigma_syst_nb")
    return _table("babar_3pi_born_spectrum", header, [(name, "f8") for name in header],
                  (_positive, _number, _positive, _positive, _positive, _positive, _positive, _number, _positive, _positive),
                  {"sqrt_s": "GeV", "efficiency": "%", "efficiency_error": "%", "dL": "nb^-1", "sigma": "nb"},
                  directory, {"sigma_syst_nb": "nb"})


def _nullable_number(value: str, name: str) -> str:
    stripped = value.strip()
    if stripped:
        _finite(stripped, name)
    return stripped


def load_kloe_phi_3pi_dalitz_fit(directory: Path = DATA_DIRECTORY) -> OddParityTable:
    header = ("parameter", "value", "stat_error", "syst_error", "unit")
    table = _table("kloe_phi_3pi_dalitz_fit", header,
                   [("parameter", "U32"), ("value", "U24"), ("stat_error", "U24"),
                    ("syst_error", "U24"), ("unit", "U16")],
                   (_text, _nullable_number, _nullable_number, _nullable_number, _text),
                   "Per unit column", directory)
    if not set(table.rows["unit"]).issubset({"MeV", "dimensionless", "rad"}):
        raise ValueError("invalid KLOE row unit")
    return table


def load_omega_pi0_gamma_table(directory: Path = DATA_DIRECTORY) -> OddParityTable:
    header = ("observable", "value", "stat_or_total_error", "unit")
    table = _table("omega_pi0_gamma_pdg", header,
                   [("observable", "U48"), ("value", "f8"), ("stat_or_total_error", "f8"), ("unit", "U16")],
                   (_text, _positive, _positive, _text), {"branching_fraction": "fraction"}, directory)
    if table.rows.size != 1 or table.rows["observable"][0] != "omega_to_pi0_gamma_branching_fraction" or table.rows["unit"][0] != "fraction":
        raise ValueError("unexpected omega branching-fraction row schema or unit")
    return table


def load_omega_pi0_gamma_observation(
    inputs: OmegaPiGammaInputs,
    GammaOmegaTotal: float,
    directory: Path = DATA_DIRECTORY,
) -> Observation:
    """Connect the single PDG branching fraction to the RChT predictor."""
    if not isinstance(inputs, OmegaPiGammaInputs):
        raise TypeError("inputs must be OmegaPiGammaInputs")
    total_width = _finite(GammaOmegaTotal, "GammaOmegaTotal", positive=True)
    table = load_omega_pi0_gamma_table(directory)
    row = table.rows[0]
    config = {name: getattr(inputs, name) for name in ("e", "F", "FV", "MV", "Momega", "mpi")}
    config["GammaOmegaTotal"] = total_width
    error = float(row["stat_or_total_error"])
    return Observation(
        [float(row["value"])], [[error**2]], omega_pi_gamma_branching_fraction_predictor,
        config, unit="fraction",
    )


def build_fixed_ope_constraint(
    parameters: OPEParameters,
    truncation_fraction: float | None = None,
    rho: float = 0.0,
) -> OPEConstraint:
    """Adapt the canonical ``ope_linear_system`` and fixed prior to a fit block."""
    matrix, target = ope_linear_system(parameters)
    prior = fixed_ope_prior(parameters, truncation_fraction=truncation_fraction, rho=rho)
    return OPEConstraint(CANONICAL_COUPLING_ORDER, matrix, target, prior.covariance, OPE_LABELS)


def build_first_release_fit_problem(
    inputs: OmegaPiGammaInputs,
    GammaOmegaTotal: float,
    ope_parameters: OPEParameters,
    parameter_starts: Mapping[str, float],
    parameter_bounds: Mapping[str, tuple[float | None, float | None]],
    *,
    truncation_fraction: float | None = None,
    rho: float = 0.0,
    directory: Path = DATA_DIRECTORY,
) -> FitProblem:
    """Build, but never run, the rank-deficient first-release problem."""
    if set(parameter_starts) != set(CANONICAL_COUPLING_ORDER):
        raise ValueError("parameter_starts must contain exactly c1..c7,d1..d4")
    if set(parameter_bounds) != set(CANONICAL_COUPLING_ORDER):
        raise ValueError("parameter_bounds must contain exactly c1..c7,d1..d4")
    specs = []
    for name in CANONICAL_COUPLING_ORDER:
        bounds = parameter_bounds[name]
        if not isinstance(bounds, (tuple, list)) or len(bounds) != 2:
            raise ValueError(f"bounds for {name} must be a (lower, upper) pair")
        specs.append(Parameter(name, parameter_starts[name], bounds[0], bounds[1]))
    return FitProblem(
        tuple(specs),
        (load_omega_pi0_gamma_observation(inputs, GammaOmegaTotal, directory),),
        ope_constraints=(build_fixed_ope_constraint(ope_parameters, truncation_fraction, rho),),
        block_independence=(
            "The PDG omega branching-fraction block and fixed OPE theory block have zero cross-covariance."
        ),
    )
