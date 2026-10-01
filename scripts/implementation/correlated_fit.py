"""Generic fixed-covariance nonlinear least-squares infrastructure.

Predictor callables live only in a runtime registry.  Problems serialize
predictor identifiers and JSON configuration, never executable objects.

Example:
    .venv/bin/python -m unittest scripts.implementation.test_correlated_fit
"""

from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import json
from types import MappingProxyType
from typing import Any, Callable, Mapping

import numpy as np
from scipy import linalg
from scipy.optimize import least_squares


SCHEMA_VERSION = "correlated-fit-problem-v3"
Predictor = Callable[[Mapping[str, float], Mapping[str, Any]], np.ndarray]


def _finite_array(value, name: str, *, ndim=None) -> np.ndarray:
    try:
        candidate = np.asarray(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{name} must be a finite real array") from error
    contains_complex_object = candidate.dtype == object and any(
        isinstance(item, (complex, np.complexfloating)) for item in candidate.flat
    )
    if np.iscomplexobj(candidate) or contains_complex_object:
        raise ValueError(f"{name} must be real; complex values are not accepted")
    try:
        result = np.array(value, dtype=float, copy=True)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{name} must be a finite real array") from error
    if ndim is not None and result.ndim != ndim:
        raise ValueError(f"{name} must be {ndim}-dimensional")
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must contain only finite values")
    result.setflags(write=False)
    return result


def _covariance(value, name: str, size: int) -> np.ndarray:
    result = _finite_array(value, name, ndim=2)
    if result.shape != (size, size):
        raise ValueError(f"{name} must have shape {(size, size)}")
    if not np.array_equal(result, result.T):
        raise ValueError(f"{name} must be exactly symmetric")
    try:
        linalg.cholesky(result, lower=True, check_finite=True)
    except linalg.LinAlgError as error:
        raise ValueError(f"{name} must be symmetric positive definite") from error
    return result


def _finite_float(value, name: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{name} must be finite and real") from error
    if not np.isfinite(result):
        raise ValueError(f"{name} must be finite and real")
    return result


def _json_value(value: Any, name: str) -> Any:
    """Validate, copy, and deeply freeze a JSON value."""
    if value is None or isinstance(value, (str, bool)):
        return value
    if isinstance(value, (int, np.integer)) and not isinstance(value, bool):
        return int(value)
    if isinstance(value, (float, np.floating)):
        return _finite_float(value, name)
    if isinstance(value, Mapping):
        if not all(isinstance(key, str) for key in value):
            raise ValueError(f"{name} object keys must be strings")
        return MappingProxyType({key: _json_value(value[key], name) for key in sorted(value)})
    if isinstance(value, (list, tuple)):
        return tuple(_json_value(item, name) for item in value)
    raise ValueError(f"{name} must contain only JSON values")


def _plain_json(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {key: _plain_json(item) for key, item in value.items()}
    if isinstance(value, tuple):
        return [_plain_json(item) for item in value]
    return value


@dataclass(frozen=True)
class NuisanceEffect:
    """One standardized Gaussian nuisance response on an observation block.

    ``multiplicative`` applies ``mu * exp(log1p(fraction) * eta)`` and
    ``additive`` applies ``additive_vector * eta``.  Additive vectors use the
    observation unit when units are represented.
    """

    source_id: str
    parameter_name: str
    response: str
    fraction: float | None = None
    additive_vector: np.ndarray | None = None
    unit: str | None = None

    def __post_init__(self) -> None:
        for name in ("source_id", "parameter_name"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip():
                raise ValueError(f"{name} must be a nonempty string")
            object.__setattr__(self, name, value.strip())
        if self.response not in ("multiplicative", "additive"):
            raise ValueError("nuisance response must be 'multiplicative' or 'additive'")
        if self.unit is not None and (not isinstance(self.unit, str) or not self.unit.strip()):
            raise ValueError("nuisance unit must be a nonempty string or None")
        if self.response == "multiplicative":
            if self.additive_vector is not None:
                raise ValueError("multiplicative nuisance must not have an additive_vector")
            fraction = _finite_float(self.fraction, "multiplicative nuisance fraction")
            if fraction < 0.0:
                raise ValueError("multiplicative nuisance fraction must be nonnegative")
            if self.unit is not None:
                raise ValueError("multiplicative nuisance is dimensionless and must not specify a unit")
            object.__setattr__(self, "fraction", fraction)
        else:
            if self.fraction is not None:
                raise ValueError("additive nuisance must not have a fraction")
            if self.additive_vector is None:
                raise ValueError("additive nuisance requires an additive_vector")
            object.__setattr__(
                self, "additive_vector", _finite_array(self.additive_vector, "additive nuisance vector", ndim=1)
            )
        if self.unit is not None:
            object.__setattr__(self, "unit", self.unit.strip())


@dataclass(frozen=True)
class ParameterSpec:
    name: str
    initial: float
    lower: float | None = None
    upper: float | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or not self.name:
            raise ValueError("parameter name must be a nonempty string")
        initial = _finite_float(self.initial, f"initial value for {self.name}")
        lower = None if self.lower is None else _finite_float(self.lower, f"lower bound for {self.name}")
        upper = None if self.upper is None else _finite_float(self.upper, f"upper bound for {self.name}")
        if lower is not None and upper is not None and lower >= upper:
            raise ValueError("parameter lower bound must be smaller than upper bound")
        if (lower is not None and initial < lower) or (upper is not None and initial > upper):
            raise ValueError("parameter initial value must lie within its bounds")
        object.__setattr__(self, "initial", initial)
        object.__setattr__(self, "lower", lower)
        object.__setattr__(self, "upper", upper)


@dataclass(frozen=True)
class ObservationGroup:
    values: np.ndarray
    covariance: np.ndarray
    predictor_id: str
    config: Mapping[str, Any] = field(default_factory=dict)
    unit: str | None = None
    diagonal_source_ids: tuple[str, ...] = ()
    nuisance_effects: tuple[NuisanceEffect, ...] = ()

    def __post_init__(self) -> None:
        values = _finite_array(self.values, "observation values", ndim=1)
        if values.size == 0:
            raise ValueError("observation values must not be empty")
        if not isinstance(self.predictor_id, str) or not self.predictor_id:
            raise ValueError("predictor_id must be a nonempty string")
        config = _json_value(self.config, "predictor config")
        if not isinstance(config, Mapping):
            raise ValueError("predictor config must be a JSON object")
        unit = self.unit
        if unit is not None and (not isinstance(unit, str) or not unit.strip()):
            raise ValueError("observation unit must be a nonempty string or None")
        unit = None if unit is None else unit.strip()
        diagonal_sources = tuple(self.diagonal_source_ids)
        if not all(isinstance(item, str) and item.strip() for item in diagonal_sources):
            raise ValueError("diagonal source IDs must be nonempty strings")
        diagonal_sources = tuple(item.strip() for item in diagonal_sources)
        if len(set(diagonal_sources)) != len(diagonal_sources):
            raise ValueError("diagonal source IDs must be unique")
        effects = tuple(self.nuisance_effects)
        if not all(isinstance(item, NuisanceEffect) for item in effects):
            raise ValueError("nuisance_effects must contain NuisanceEffect objects")
        effect_sources = tuple(item.source_id for item in effects)
        if len(set(effect_sources)) != len(effect_sources):
            raise ValueError("nuisance source IDs must be unique within an observation group")
        overlap = set(diagonal_sources) & set(effect_sources)
        if overlap:
            raise ValueError(f"source_id cannot be both diagonal and nuisance: {sorted(overlap)[0]}")
        for effect in effects:
            if effect.response == "additive":
                if effect.additive_vector.shape != values.shape:
                    raise ValueError("additive nuisance vector length must match observation values")
                if unit is not None and effect.unit != unit:
                    raise ValueError("additive nuisance unit must match the observation unit")
                if unit is None and effect.unit is not None:
                    raise ValueError("additive nuisance unit cannot be represented without an observation unit")
        object.__setattr__(self, "values", values)
        object.__setattr__(self, "covariance", _covariance(self.covariance, "observation covariance", values.size))
        object.__setattr__(self, "config", config)
        object.__setattr__(self, "unit", unit)
        object.__setattr__(self, "diagonal_source_ids", diagonal_sources)
        object.__setattr__(self, "nuisance_effects", effects)


@dataclass(frozen=True)
class ParameterPrior:
    parameter_names: tuple[str, ...]
    means: np.ndarray
    covariance: np.ndarray

    def __post_init__(self) -> None:
        names = tuple(self.parameter_names)
        if not names or not all(isinstance(name, str) and name for name in names) or len(set(names)) != len(names):
            raise ValueError("prior parameter_names must be unique nonempty strings")
        means = _finite_array(self.means, "prior means", ndim=1)
        if means.size != len(names):
            raise ValueError("prior means must match parameter_names")
        object.__setattr__(self, "parameter_names", names)
        object.__setattr__(self, "means", means)
        object.__setattr__(self, "covariance", _covariance(self.covariance, "prior covariance", len(names)))


@dataclass(frozen=True)
class FixedOPEConstraint:
    """Fixed linear OPE residual ``matrix @ parameters - target``."""

    parameter_names: tuple[str, ...]
    matrix: np.ndarray
    target: np.ndarray
    covariance: np.ndarray
    labels: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        names = tuple(self.parameter_names)
        if not names or not all(isinstance(name, str) and name for name in names) or len(set(names)) != len(names):
            raise ValueError("OPE parameter_names must be unique nonempty strings")
        matrix = _finite_array(self.matrix, "OPE matrix", ndim=2)
        target = _finite_array(self.target, "OPE target", ndim=1)
        if matrix.shape != (target.size, len(names)) or target.size == 0:
            raise ValueError("OPE matrix shape must match target and parameter_names")
        labels = tuple(self.labels) if self.labels else tuple(f"ope-{index}" for index in range(target.size))
        if len(labels) != target.size or not all(isinstance(label, str) and label for label in labels):
            raise ValueError("OPE labels must match the target")
        object.__setattr__(self, "parameter_names", names)
        object.__setattr__(self, "matrix", matrix)
        object.__setattr__(self, "target", target)
        object.__setattr__(self, "covariance", _covariance(self.covariance, "OPE covariance", target.size))
        object.__setattr__(self, "labels", labels)


@dataclass(frozen=True)
class FitProblem:
    parameters: tuple[ParameterSpec, ...]
    observations: tuple[ObservationGroup, ...]
    priors: tuple[ParameterPrior, ...] = ()
    ope_constraints: tuple[FixedOPEConstraint, ...] = ()
    block_independence_assumption: str | None = None

    def __post_init__(self) -> None:
        parameters = tuple(self.parameters)
        observations = tuple(self.observations)
        priors = tuple(self.priors)
        constraints = tuple(self.ope_constraints)
        if not parameters or not all(isinstance(item, ParameterSpec) for item in parameters):
            raise ValueError("parameters must contain ParameterSpec objects")
        if not observations or not all(isinstance(item, ObservationGroup) for item in observations):
            raise ValueError("observations must contain ObservationGroup objects")
        names = tuple(item.name for item in parameters)
        if len(set(names)) != len(names):
            raise ValueError("parameter names must be unique")
        available = set(names)
        if not all(isinstance(item, ParameterPrior) for item in priors):
            raise ValueError("priors must contain ParameterPrior objects")
        if not all(isinstance(item, FixedOPEConstraint) for item in constraints):
            raise ValueError("ope_constraints must contain FixedOPEConstraint objects")
        block_count = len(observations) + len(priors) + len(constraints)
        assumption = self.block_independence_assumption
        if assumption is None:
            if block_count > 1:
                raise ValueError(
                    "block_independence_assumption is required for multiple residual covariance blocks"
                )
        elif not isinstance(assumption, str) or not assumption.strip():
            raise ValueError("block_independence_assumption must be a nonempty string or None")
        else:
            assumption = assumption.strip()
        referenced = [name for item in priors + constraints for name in item.parameter_names]
        if not set(referenced).issubset(available):
            raise ValueError("prior and OPE parameter names must occur in parameters")
        nuisance_parameter_by_source: dict[str, str] = {}
        nuisance_names: set[str] = set()
        all_diagonal_sources = {
            source_id for observation in observations for source_id in observation.diagonal_source_ids
        }
        for observation in observations:
            for effect in observation.nuisance_effects:
                if effect.source_id in all_diagonal_sources:
                    raise ValueError(f"source_id cannot be both diagonal and nuisance: {effect.source_id}")
                previous = nuisance_parameter_by_source.setdefault(effect.source_id, effect.parameter_name)
                if previous != effect.parameter_name:
                    raise ValueError("a shared nuisance source_id must use one consistent parameter name")
                nuisance_names.add(effect.parameter_name)
        if not nuisance_names.issubset(available):
            raise ValueError("nuisance parameter names must occur in parameters")
        prior_occurrences: dict[str, list[tuple[float, float]]] = {name: [] for name in nuisance_names}
        for prior in priors:
            for index, name in enumerate(prior.parameter_names):
                if name in prior_occurrences:
                    prior_occurrences[name].append((float(prior.means[index]), float(prior.covariance[index, index])))
        for name, occurrences in prior_occurrences.items():
            if len(occurrences) != 1:
                raise ValueError(f"nuisance parameter {name!r} requires exactly one prior")
            mean, variance = occurrences[0]
            if mean != 0.0 or variance != 1.0:
                raise ValueError(f"nuisance parameter {name!r} prior must be standardized N(0,1)")
        object.__setattr__(self, "parameters", parameters)
        object.__setattr__(self, "observations", observations)
        object.__setattr__(self, "priors", priors)
        object.__setattr__(self, "ope_constraints", constraints)
        object.__setattr__(self, "block_independence_assumption", assumption)

    def canonical_json(self) -> str:
        return serialize_fit_problem(self)

    def sha256(self) -> str:
        return fit_problem_digest(self)

    @property
    def ope_design_rank(self) -> int:
        """Rank of all OPE rows embedded in the canonical parameter order."""
        if not self.ope_constraints:
            return 0
        names = tuple(item.name for item in self.parameters)
        index = {name: position for position, name in enumerate(names)}
        rows = []
        for constraint in self.ope_constraints:
            embedded = np.zeros((constraint.matrix.shape[0], len(names)))
            embedded[:, [index[name] for name in constraint.parameter_names]] = constraint.matrix
            rows.append(embedded)
        return int(np.linalg.matrix_rank(np.vstack(rows)))

    @property
    def ope_rank_deficiency(self) -> int:
        """Parameter-space deficiency left by the fixed OPE rows alone."""
        return len(self.parameters) - self.ope_design_rank


@dataclass(frozen=True)
class JacobianDiagnostics:
    rank: int
    singular_values: np.ndarray
    null_vectors: np.ndarray
    condition_number: float


@dataclass(frozen=True)
class ComponentChi2:
    data: tuple[float, ...]
    priors: tuple[float, ...]
    ope: tuple[float, ...]
    total: float


@dataclass(frozen=True)
class FitDiagnostics:
    data_jacobian: JacobianDiagnostics
    full_jacobian: JacobianDiagnostics
    active_bounds: tuple[str, ...]
    chi2: ComponentChi2


@dataclass(frozen=True)
class FitResult:
    parameter_names: tuple[str, ...]
    values: np.ndarray
    success: bool
    status: int
    message: str
    nfev: int
    residuals: np.ndarray
    diagnostics: FitDiagnostics
    problem_sha256: str

    def as_mapping(self) -> Mapping[str, float]:
        return MappingProxyType(dict(zip(self.parameter_names, self.values)))


class PredictorRegistry:
    """Mutable runtime-only map from stable identifiers to predictors."""

    def __init__(self) -> None:
        self._predictors: dict[str, Predictor] = {}

    def register(self, predictor_id: str, predictor: Predictor) -> None:
        if not isinstance(predictor_id, str) or not predictor_id:
            raise ValueError("predictor_id must be a nonempty string")
        if predictor_id in self._predictors:
            raise ValueError(f"predictor {predictor_id!r} is already registered")
        if not callable(predictor):
            raise TypeError("predictor must be callable")
        self._predictors[predictor_id] = predictor

    def predict(self, predictor_id: str, parameters: Mapping[str, float], config: Mapping[str, Any]) -> np.ndarray:
        try:
            predictor = self._predictors[predictor_id]
        except KeyError as error:
            raise KeyError(f"predictor {predictor_id!r} is not registered") from error
        return _finite_array(predictor(parameters, config), f"prediction from {predictor_id}", ndim=1)


def _problem_object(problem: FitProblem) -> dict[str, Any]:
    return {
        "schema_version": SCHEMA_VERSION,
        "block_independence_assumption": problem.block_independence_assumption,
        "parameters": [
            {"name": item.name, "initial": item.initial, "lower": item.lower, "upper": item.upper}
            for item in problem.parameters
        ],
        "observations": [
            {"values": item.values.tolist(), "covariance": item.covariance.tolist(),
             "predictor_id": item.predictor_id, "config": _plain_json(item.config),
             "unit": item.unit, "diagonal_source_ids": list(item.diagonal_source_ids),
             "nuisance_effects": [
                 {"source_id": effect.source_id, "parameter_name": effect.parameter_name,
                  "response": effect.response, "fraction": effect.fraction,
                  "additive_vector": None if effect.additive_vector is None else effect.additive_vector.tolist(),
                  "unit": effect.unit}
                 for effect in item.nuisance_effects
             ]}
            for item in problem.observations
        ],
        "priors": [
            {"parameter_names": list(item.parameter_names), "means": item.means.tolist(),
             "covariance": item.covariance.tolist()}
            for item in problem.priors
        ],
        "ope_constraints": [
            {"parameter_names": list(item.parameter_names), "matrix": item.matrix.tolist(),
             "target": item.target.tolist(), "covariance": item.covariance.tolist(),
             "labels": list(item.labels)}
            for item in problem.ope_constraints
        ],
    }


def serialize_fit_problem(problem: FitProblem) -> str:
    if not isinstance(problem, FitProblem):
        raise TypeError("problem must be a FitProblem")
    return json.dumps(_problem_object(problem), sort_keys=True, separators=(",", ":"), allow_nan=False)


def fit_problem_digest(problem: FitProblem) -> str:
    return hashlib.sha256(serialize_fit_problem(problem).encode("utf-8")).hexdigest()


def apply_nuisance_effects(
    prediction: np.ndarray,
    effects: tuple[NuisanceEffect, ...],
    parameters: Mapping[str, float],
) -> np.ndarray:
    """Apply multiplicative effects first, then additive effects, deterministically."""
    result = np.array(_finite_array(prediction, "prediction", ndim=1), copy=True)
    effects = tuple(effects)
    if not all(isinstance(effect, NuisanceEffect) for effect in effects):
        raise TypeError("effects must contain NuisanceEffect objects")
    for effect in effects:
        if effect.parameter_name not in parameters:
            raise ValueError(f"missing nuisance parameter {effect.parameter_name!r}")
        eta = _finite_float(parameters[effect.parameter_name], effect.parameter_name)
        if effect.response == "multiplicative":
            with np.errstate(over="ignore", invalid="ignore"):
                result *= np.exp(np.log1p(effect.fraction) * eta)
    for effect in effects:
        if effect.response == "additive":
            eta = _finite_float(parameters[effect.parameter_name], effect.parameter_name)
            if effect.additive_vector.shape != result.shape:
                raise ValueError("additive nuisance vector length must match prediction")
            result += effect.additive_vector * eta
    if not np.all(np.isfinite(result)):
        raise ValueError("nuisance-adjusted prediction must be finite")
    result.setflags(write=False)
    return result


def _jacobian_diagnostics(jacobian: np.ndarray, parameter_count: int) -> JacobianDiagnostics:
    matrix = np.asarray(jacobian, dtype=float)
    _, singular_values, vh = linalg.svd(matrix, full_matrices=True, check_finite=True)
    tolerance = max(matrix.shape) * np.finfo(float).eps * (singular_values[0] if singular_values.size else 0.0)
    rank = int(np.sum(singular_values > tolerance))
    condition = float(np.inf if rank < parameter_count else singular_values[0] / singular_values[-1])
    singular_values = _finite_array(singular_values, "singular values", ndim=1)
    null_vectors = np.array(vh[rank:, :], copy=True)
    null_vectors.setflags(write=False)
    return JacobianDiagnostics(rank, singular_values, null_vectors, condition)


def fit_correlated(problem: FitProblem, registry: PredictorRegistry) -> FitResult:
    """Fit with fixed Gaussian covariances and ordinary quadratic loss."""
    if not isinstance(problem, FitProblem) or not isinstance(registry, PredictorRegistry):
        raise TypeError("problem and registry must be FitProblem and PredictorRegistry")
    names = tuple(item.name for item in problem.parameters)
    name_to_index = {name: index for index, name in enumerate(names)}
    x0 = np.array([item.initial for item in problem.parameters])
    lower = np.array([-np.inf if item.lower is None else item.lower for item in problem.parameters])
    upper = np.array([np.inf if item.upper is None else item.upper for item in problem.parameters])
    data_factors = [linalg.cholesky(item.covariance, lower=True) for item in problem.observations]
    prior_factors = [linalg.cholesky(item.covariance, lower=True) for item in problem.priors]
    ope_factors = [linalg.cholesky(item.covariance, lower=True) for item in problem.ope_constraints]
    data_size = sum(item.values.size for item in problem.observations)

    def components(x: np.ndarray) -> tuple[list[np.ndarray], list[np.ndarray], list[np.ndarray]]:
        parameter_map = MappingProxyType(dict(zip(names, x)))
        data_parts = []
        for item, factor in zip(problem.observations, data_factors):
            prediction = registry.predict(item.predictor_id, parameter_map, item.config)
            if prediction.shape != item.values.shape:
                raise ValueError(f"predictor {item.predictor_id!r} returned shape {prediction.shape}, expected {item.values.shape}")
            prediction = apply_nuisance_effects(prediction, item.nuisance_effects, parameter_map)
            data_parts.append(linalg.solve_triangular(factor, prediction - item.values, lower=True, check_finite=True))
        prior_parts = []
        for item, factor in zip(problem.priors, prior_factors):
            selected = x[[name_to_index[name] for name in item.parameter_names]]
            prior_parts.append(linalg.solve_triangular(factor, selected - item.means, lower=True, check_finite=True))
        ope_parts = []
        for item, factor in zip(problem.ope_constraints, ope_factors):
            selected = x[[name_to_index[name] for name in item.parameter_names]]
            ope_parts.append(linalg.solve_triangular(factor, item.matrix @ selected - item.target, lower=True, check_finite=True))
        return data_parts, prior_parts, ope_parts

    def residual(x: np.ndarray) -> np.ndarray:
        parts = components(x)
        return np.concatenate(parts[0] + parts[1] + parts[2])

    optimum = least_squares(residual, x0, bounds=(lower, upper), loss="linear")
    data_parts, prior_parts, ope_parts = components(optimum.x)
    all_parts = data_parts + prior_parts + ope_parts
    chi2 = ComponentChi2(
        tuple(float(part @ part) for part in data_parts),
        tuple(float(part @ part) for part in prior_parts),
        tuple(float(part @ part) for part in ope_parts),
        float(sum(part @ part for part in all_parts)),
    )
    active = tuple(
        f"{name}:{'lower' if mask < 0 else 'upper'}"
        for name, mask in zip(names, optimum.active_mask) if mask
    )
    values = _finite_array(optimum.x, "fit values", ndim=1)
    residuals = _finite_array(optimum.fun, "fit residuals", ndim=1)
    diagnostics = FitDiagnostics(
        _jacobian_diagnostics(optimum.jac[:data_size, :], len(names)),
        _jacobian_diagnostics(optimum.jac, len(names)),
        active,
        chi2,
    )
    return FitResult(names, values, bool(optimum.success), int(optimum.status), str(optimum.message),
                     int(optimum.nfev), residuals, diagnostics, fit_problem_digest(problem))
