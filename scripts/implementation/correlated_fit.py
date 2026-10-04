"""Fixed-covariance Gaussian least-squares fitting infrastructure.

The method is fixed-covariance whitening plus ``scipy.optimize.least_squares``
with ordinary (linear) loss.  Every residual block is whitened with the Cholesky
factor ``L`` of its covariance ``C = L L^T`` taken with ``lower=True``, the
whitened blocks are concatenated, and their sum of squares is minimized.
Covariances are inputs: each ``L`` is computed once before the minimization and
is never recomputed, so no covariance changes during the fit.  There is no
jitter, no pseudoinverse, no robust loss, and no implicit diagonalization; a
covariance that is not exactly symmetric positive definite is rejected, not
repaired.  Each ``Nuisance`` declares one standardized Gaussian ``N(0, 1)``
parameter, so that prior is implicit in the declaration and contributes to chi2
through a residual row equal to ``eta`` itself.

Invocation:
    .venv/bin/python -m unittest scripts.implementation.test_correlated_fit
"""

from __future__ import annotations

from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Any, Callable, Mapping, Sequence

import numpy as np
from scipy import linalg
from scipy.optimize import least_squares


Predictor = Callable[[Mapping[str, float], Mapping[str, Any]], np.ndarray]


def _finite_array(value, name: str, *, ndim: int | None = None) -> np.ndarray:
    """Return an owned finite real copy of ``value``, rejecting complex input."""
    if np.iscomplexobj(np.asarray(value)):
        raise ValueError(f"{name} must be real; complex values are not accepted")
    try:
        result = np.array(value, dtype=float, copy=True)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{name} must be a finite real array") from error
    if ndim is not None and result.ndim != ndim:
        raise ValueError(f"{name} must be {ndim}-dimensional")
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must contain only finite values")
    return result


def _covariance(value, name: str, size: int) -> np.ndarray:
    """Return one fixed covariance block, rejecting anything not exactly SPD."""
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


@dataclass(frozen=True)
class Parameter:
    name: str
    start: float
    lower: float | None = None
    upper: float | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or not self.name:
            raise ValueError("parameter name must be a nonempty string")
        start = _finite_float(self.start, f"start value for {self.name}")
        lower = None if self.lower is None else _finite_float(self.lower, f"lower bound for {self.name}")
        upper = None if self.upper is None else _finite_float(self.upper, f"upper bound for {self.name}")
        if lower is not None and upper is not None and lower >= upper:
            raise ValueError("parameter lower bound must be smaller than upper bound")
        if (lower is not None and start < lower) or (upper is not None and start > upper):
            raise ValueError("parameter start must lie within its bounds")
        for attribute, value in (("start", start), ("lower", lower), ("upper", upper)):
            object.__setattr__(self, attribute, value)


@dataclass(frozen=True)
class Nuisance:
    """One standardized Gaussian N(0,1) parameter acting on one observation block.

    ``fraction`` gives the multiplicative response ``mu * exp(log1p(fraction) * eta)``
    (positive for every eta); ``shift`` gives the additive response ``shift * eta``
    in the observation's ``unit``.  The N(0,1) prior is implicit in the declaration
    and is contributed to chi2 by a residual row equal to ``eta`` itself.
    """

    parameter: str
    fraction: float | None = None
    shift: Sequence[float] | None = None
    unit: str | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.parameter, str) or not self.parameter.strip():
            raise ValueError("nuisance parameter must be a nonempty string")
        object.__setattr__(self, "parameter", self.parameter.strip())
        if (self.fraction is None) == (self.shift is None):
            raise ValueError("a nuisance must declare exactly one of fraction or shift")
        if self.fraction is not None:
            fraction = _finite_float(self.fraction, "multiplicative nuisance fraction")
            if fraction < 0.0:
                raise ValueError("multiplicative nuisance fraction must be nonnegative")
            if self.unit is not None:
                raise ValueError("multiplicative nuisance is dimensionless and must not specify a unit")
            object.__setattr__(self, "fraction", fraction)
            return
        if self.unit is not None and (not isinstance(self.unit, str) or not self.unit.strip()):
            raise ValueError("nuisance unit must be a nonempty string or None")
        object.__setattr__(self, "unit", None if self.unit is None else self.unit.strip())
        object.__setattr__(self, "shift", _finite_array(self.shift, "additive nuisance shift", ndim=1))


@dataclass(frozen=True)
class Observation:
    """One data block: measured ``values`` with a fixed ``covariance``."""

    values: np.ndarray
    covariance: np.ndarray
    predict: Predictor
    config: Mapping[str, Any] = field(default_factory=dict)
    unit: str | None = None
    nuisances: tuple[Nuisance, ...] = ()

    def __post_init__(self) -> None:
        values = _finite_array(self.values, "observation values", ndim=1)
        if values.size == 0:
            raise ValueError("observation values must not be empty")
        if not callable(self.predict):
            raise TypeError("observation predictor must be callable")
        if not isinstance(self.config, Mapping):
            raise ValueError("observation config must be a mapping")
        if self.unit is not None and (not isinstance(self.unit, str) or not self.unit.strip()):
            raise ValueError("observation unit must be a nonempty string or None")
        unit = None if self.unit is None else self.unit.strip()
        nuisances = tuple(self.nuisances)
        if not all(isinstance(item, Nuisance) for item in nuisances):
            raise ValueError("observation nuisances must contain Nuisance objects")
        for nuisance in nuisances:
            if nuisance.shift is None:
                continue
            if nuisance.shift.shape != values.shape:
                raise ValueError("additive nuisance shift length must match observation values")
            if nuisance.unit != unit:
                raise ValueError("additive nuisance unit must match the observation unit")
        object.__setattr__(self, "values", values)
        object.__setattr__(self, "covariance", _covariance(self.covariance, "observation covariance", values.size))
        object.__setattr__(self, "unit", unit)
        object.__setattr__(self, "nuisances", nuisances)


@dataclass(frozen=True)
class OPEConstraint:
    """Fixed linear theory residual ``matrix @ theta - target``."""

    parameter_names: tuple[str, ...]
    matrix: np.ndarray
    target: np.ndarray
    covariance: np.ndarray
    labels: tuple[str, ...] = ()  # default: f"ope-{i}"

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
    """Parameters, correlated data blocks, and optional fixed OPE theory blocks.

    ``block_independence`` states the zero-cross-covariance assumption between
    the residual covariance blocks (observations plus OPE constraints plus the
    declared standardized nuisance priors).  It is required as a nonempty string
    when more than one such block is present; with a single block there is no
    cross-covariance question to answer and it may be left ``None``.
    """

    parameters: tuple[Parameter, ...]
    observations: tuple[Observation, ...]
    ope_constraints: tuple[OPEConstraint, ...] = ()
    block_independence: str | None = None

    def __post_init__(self) -> None:
        parameters, observations = tuple(self.parameters), tuple(self.observations)
        constraints = tuple(self.ope_constraints)
        if not parameters or not all(isinstance(item, Parameter) for item in parameters):
            raise ValueError("parameters must contain Parameter objects")
        if not observations or not all(isinstance(item, Observation) for item in observations):
            raise ValueError("observations must contain Observation objects")
        if not all(isinstance(item, OPEConstraint) for item in constraints):
            raise ValueError("ope_constraints must contain OPEConstraint objects")
        names = tuple(item.name for item in parameters)
        if len(set(names)) != len(names):
            raise ValueError("parameter names must be unique")
        declared = [nuisance.parameter for item in observations for nuisance in item.nuisances]
        referenced = set(declared) | {name for item in constraints for name in item.parameter_names}
        if not referenced.issubset(names):
            raise ValueError("OPE and nuisance parameter names must occur in parameters")
        duplicated = sorted({name for name in declared if declared.count(name) != 1})
        if duplicated:
            raise ValueError(f"nuisance parameter {duplicated[0]!r} must be declared by exactly one nuisance")
        statement = self.block_independence
        if statement is not None:
            if not isinstance(statement, str) or not statement.strip():
                raise ValueError("block_independence must be a nonempty string or None")
            statement = statement.strip()
        blocks = len(observations) + len(constraints) + len(set(declared))
        if blocks > 1 and statement is None:
            raise ValueError("block_independence is required for multiple residual covariance blocks")
        object.__setattr__(self, "parameters", parameters)
        object.__setattr__(self, "observations", observations)
        object.__setattr__(self, "ope_constraints", constraints)
        object.__setattr__(self, "block_independence", statement)

    @property
    def ope_design_rank(self) -> int:
        """Rank of all OPE rows embedded in the canonical parameter order."""
        if not self.ope_constraints:
            return 0
        index = {name: position for position, name in enumerate(item.name for item in self.parameters)}
        embedded_rows = []
        for constraint in self.ope_constraints:
            embedded = np.zeros((constraint.matrix.shape[0], len(self.parameters)))
            embedded[:, [index[name] for name in constraint.parameter_names]] = constraint.matrix
            embedded_rows.append(embedded)
        return int(np.linalg.matrix_rank(np.vstack(embedded_rows)))

    @property
    def ope_rank_deficiency(self) -> int:
        """Parameter-space deficiency left by the fixed OPE rows alone."""
        return len(self.parameters) - self.ope_design_rank


@dataclass(frozen=True)
class RankInfo:
    """Numerical rank report for a Jacobian block."""

    rank: int
    singular_values: np.ndarray
    null_directions: np.ndarray  # (n_parameters - rank, n_parameters); each ROW is a null direction
    condition_number: float  # infinite when rank < n_parameters


@dataclass(frozen=True)
class Chi2:
    data: float
    nuisances: float
    theory: float
    total: float


@dataclass(frozen=True)
class FitResult:
    parameter_names: tuple[str, ...]
    values: np.ndarray
    success: bool
    message: str
    nfev: int
    residuals: np.ndarray
    data_rank: RankInfo
    full_rank: RankInfo
    active_bounds: tuple[str, ...]
    chi2: Chi2

    def as_mapping(self) -> Mapping[str, float]:
        return MappingProxyType(dict(zip(self.parameter_names, self.values)))


def apply_nuisances(prediction: np.ndarray, nuisances: tuple[Nuisance, ...], parameters: Mapping[str, float]) -> np.ndarray:
    """Apply all multiplicative responses first, then all additive ones."""
    result = np.array(_finite_array(prediction, "prediction", ndim=1), copy=True)
    ordered = tuple(nuisances)
    for nuisance in ordered:
        if nuisance.fraction is not None:
            with np.errstate(over="ignore", invalid="ignore"):
                result *= np.exp(np.log1p(nuisance.fraction) * _standardized(nuisance, parameters))
    for nuisance in ordered:
        if nuisance.shift is not None:
            if nuisance.shift.shape != result.shape:
                raise ValueError("additive nuisance shift length must match the prediction")
            result += nuisance.shift * _standardized(nuisance, parameters)
    if not np.all(np.isfinite(result)):
        raise ValueError("nuisance-adjusted prediction must be finite")
    return result


def _standardized(nuisance: Nuisance, parameters: Mapping[str, float]) -> float:
    """Return the standardized nuisance value ``eta`` carried by ``parameters``."""
    if nuisance.parameter not in parameters:
        raise ValueError(f"missing nuisance parameter {nuisance.parameter!r}")
    return _finite_float(parameters[nuisance.parameter], nuisance.parameter)


def _rank_info(jacobian: np.ndarray, parameter_count: int) -> RankInfo:
    """Report the numerical rank, null directions, and conditioning of ``jacobian``."""
    _, singular_values, vh = linalg.svd(jacobian, full_matrices=True, check_finite=True)
    tolerance = max(jacobian.shape) * np.finfo(float).eps * (singular_values[0] if singular_values.size else 0.0)
    rank = int(np.sum(singular_values > tolerance))
    condition = float(np.inf if rank < parameter_count else singular_values[0] / singular_values[-1])
    return RankInfo(rank, singular_values, np.array(vh[rank:, :], copy=True), condition)


def fit_correlated(problem: FitProblem) -> FitResult:
    """Minimize the whitened sum of squares with fixed covariances and linear loss.

    The residual vector concatenates, in order, one ``solve_triangular(L, ...)``
    block per observation, one row equal to ``eta`` per declared nuisance
    parameter, and one ``solve_triangular(L, matrix @ theta - target)`` block per
    OPE constraint.  Every Cholesky factor is computed here, once, before the
    minimization, and reused.
    """
    if not isinstance(problem, FitProblem):
        raise TypeError("problem must be a FitProblem")
    names = tuple(item.name for item in problem.parameters)
    name_to_index = {name: index for index, name in enumerate(names)}
    x0 = np.array([item.start for item in problem.parameters])
    lower = np.array([-np.inf if item.lower is None else item.lower for item in problem.parameters])
    upper = np.array([np.inf if item.upper is None else item.upper for item in problem.parameters])
    data_factors = [linalg.cholesky(item.covariance, lower=True) for item in problem.observations]
    ope_factors = [linalg.cholesky(item.covariance, lower=True) for item in problem.ope_constraints]
    nuisance_order = tuple(dict.fromkeys(n.parameter for block in problem.observations for n in block.nuisances))
    data_size = sum(item.values.size for item in problem.observations)

    def components(x: np.ndarray) -> tuple[list[np.ndarray], list[np.ndarray], list[np.ndarray]]:
        parameter_map = dict(zip(names, x))
        data_parts = []
        for item, factor in zip(problem.observations, data_factors):
            prediction = _finite_array(item.predict(parameter_map, item.config), "observation prediction", ndim=1)
            if prediction.shape != item.values.shape:
                raise ValueError(f"predictor returned shape {prediction.shape}, expected {item.values.shape}")
            prediction = apply_nuisances(prediction, item.nuisances, parameter_map)
            data_parts.append(linalg.solve_triangular(factor, prediction - item.values, lower=True, check_finite=True))
        nuisance_parts = [np.array([x[name_to_index[name]]]) for name in nuisance_order]
        theory_parts = []
        for item, factor in zip(problem.ope_constraints, ope_factors):
            selected = x[[name_to_index[name] for name in item.parameter_names]]
            theory_parts.append(
                linalg.solve_triangular(factor, item.matrix @ selected - item.target, lower=True, check_finite=True)
            )
        return data_parts, nuisance_parts, theory_parts

    def residual(x: np.ndarray) -> np.ndarray:
        return np.concatenate([part for group in components(x) for part in group])

    optimum = least_squares(residual, x0, bounds=(lower, upper), loss="linear")
    data_parts, nuisance_parts, theory_parts = components(optimum.x)
    chi2 = Chi2(
        *[float(sum(part @ part for part in group)) for group in (data_parts, nuisance_parts, theory_parts)],
        float(optimum.fun @ optimum.fun),
    )
    active_bounds = tuple(f"{name}:{'lower' if mask < 0 else 'upper'}" for name, mask in zip(names, optimum.active_mask) if mask)
    return FitResult(
        names, np.asarray(optimum.x, dtype=float), bool(optimum.success), str(optimum.message), int(optimum.nfev),
        np.asarray(optimum.fun, dtype=float), _rank_info(optimum.jac[:data_size], len(names)),
        _rank_info(optimum.jac, len(names)), active_bounds, chi2,
    )