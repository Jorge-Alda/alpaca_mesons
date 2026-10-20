"""Numerical kernels for the non-HLS odd-parity one-multiplet RChT basis.

This module exposes flavor coefficients and short-distance residuals, not
on-shell amplitudes, widths, fitted couplings, or a conversion from unordered
``C_AB`` coefficients to an anomalous bilinear matrix ``Lambda``.

Example:
    .venv/bin/python -m unittest scripts.implementation.test_odd_parity_rcht
"""

from dataclasses import dataclass

import numpy as np
from scipy import linalg


OPE_LABELS = (
    "4*c3+c1",
    "c1-c2+c5",
    "c5-c6-short_distance",
    "d1+8*d2-short_distance",
    "d3-short_distance",
)
IDEAL_VECTOR_MIXING_LABEL = "ideal omega-phi mixing"


def _readonly_array(value, name: str, *, shape=None, ndim=None, real=False) -> np.ndarray:
    if real:
        candidate = np.asarray(value)
        if np.iscomplexobj(candidate) and np.any(np.imag(candidate) != 0.0):
            raise ValueError(f"{name} must be real")
    try:
        array = np.array(value, dtype=float if real else complex, copy=True)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{name} must be a finite {'real ' if real else ''}array") from error
    if shape is not None and array.shape != shape:
        raise ValueError(f"{name} must have shape {shape}")
    if ndim is not None and array.ndim != ndim:
        raise ValueError(f"{name} must be {ndim}-dimensional")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values")
    array.setflags(write=False)
    return array


def _hermitian_matrix(value, name: str) -> np.ndarray:
    array = _readonly_array(value, name, shape=(3, 3))
    if not np.allclose(array, array.conj().T, rtol=0.0, atol=1e-12):
        raise ValueError(f"{name} must be Hermitian")
    return array


def _source_matrix(value, name: str) -> np.ndarray:
    return _readonly_array(value, name, shape=(3, 3))


def _positive_real(value, name: str) -> float:
    if isinstance(value, (complex, np.complexfloating)) and value.imag != 0.0:
        raise ValueError(f"{name} must be real, finite, and positive")
    result = float(np.real(value))
    if not np.isfinite(result) or result <= 0.0:
        raise ValueError(f"{name} must be real, finite, and positive")
    return result


def _finite_real(value, name: str) -> float:
    if isinstance(value, (complex, np.complexfloating)) and value.imag != 0.0:
        raise ValueError(f"{name} must be finite and real")
    result = float(np.real(value))
    if not np.isfinite(result):
        raise ValueError(f"{name} must be finite and real")
    return result


def _finite_scalar(value, name: str) -> complex:
    try:
        result = complex(value)
    except (TypeError, ValueError) as error:
        raise TypeError(f"{name} must be a scalar") from error
    if not np.isfinite(result.real) or not np.isfinite(result.imag):
        raise ValueError(f"{name} must be finite")
    return result


@dataclass(frozen=True)
class ALPFlavorSpurions:
    """Physical-ALP flavor inputs; ``b0_mq`` already includes ``B0``."""

    m_a_inverse_C: np.ndarray
    c_q_AP: np.ndarray
    b0_mq: np.ndarray

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "m_a_inverse_C",
            _hermitian_matrix(self.m_a_inverse_C, "m_a_inverse_C"),
        )
        object.__setattr__(self, "c_q_AP", _hermitian_matrix(self.c_q_AP, "c_q_AP"))
        object.__setattr__(self, "b0_mq", _hermitian_matrix(self.b0_mq, "b0_mq"))


@dataclass(frozen=True)
class VJPCouplings:
    c1: float
    c2: float
    c3: float
    c4: float
    c5: float
    c6: float
    c7: float

    def __post_init__(self) -> None:
        for name in self.__dataclass_fields__:
            object.__setattr__(self, name, _finite_real(getattr(self, name), name))

    def as_array(self) -> np.ndarray:
        return _readonly_array(tuple(getattr(self, f"c{i}") for i in range(1, 8)), "vjp", ndim=1, real=True)


@dataclass(frozen=True)
class VVPCouplings:
    d1: float
    d2: float
    d3: float
    d4: float

    def __post_init__(self) -> None:
        for name in self.__dataclass_fields__:
            object.__setattr__(self, name, _finite_real(getattr(self, name), name))

    def as_array(self) -> np.ndarray:
        return _readonly_array(tuple(getattr(self, f"d{i}") for i in range(1, 5)), "vvp", ndim=1, real=True)


@dataclass(frozen=True)
class OddParityCouplings:
    vjp: VJPCouplings
    vvp: VVPCouplings

    def __post_init__(self) -> None:
        if not isinstance(self.vjp, VJPCouplings) or not isinstance(self.vvp, VVPCouplings):
            raise TypeError("vjp and vvp must be VJPCouplings and VVPCouplings")

    def as_array(self) -> np.ndarray:
        return _readonly_array(np.concatenate((self.vjp.as_array(), self.vvp.as_array())), "couplings", ndim=1, real=True)


@dataclass(frozen=True)
class OddParityFlavorKernels:
    derivative_trace: complex
    mass_trace: complex

    def __post_init__(self) -> None:
        object.__setattr__(self, "derivative_trace", _finite_scalar(self.derivative_trace, "derivative_trace"))
        object.__setattr__(self, "mass_trace", _finite_scalar(self.mass_trace, "mass_trace"))


@dataclass(frozen=True)
class OperatorCoefficients:
    """Ordered, unreduced coefficients of distinct Lorentz operators."""

    coefficients: np.ndarray

    def __post_init__(self) -> None:
        object.__setattr__(self, "coefficients", _readonly_array(self.coefficients, "coefficients", ndim=1))


@dataclass(frozen=True)
class OPEParameters:
    n_c: int
    m_v: float
    f_v: float
    f: float

    def __post_init__(self) -> None:
        if isinstance(self.n_c, bool) or not isinstance(self.n_c, (int, np.integer)) or self.n_c <= 0:
            raise ValueError("n_c must be a positive integer")
        object.__setattr__(self, "n_c", int(self.n_c))
        object.__setattr__(self, "m_v", _positive_real(self.m_v, "m_v"))
        object.__setattr__(self, "f_v", _positive_real(self.f_v, "f_v"))
        object.__setattr__(self, "f", _positive_real(self.f, "f"))


@dataclass(frozen=True)
class OPEResult:
    labels: tuple[str, ...]
    residuals: np.ndarray

    def __post_init__(self) -> None:
        labels = tuple(self.labels)
        residuals = _readonly_array(self.residuals, "residuals", ndim=1, real=True)
        if len(labels) != residuals.size or not all(isinstance(label, str) for label in labels):
            raise ValueError("labels must be strings matching the residual count")
        object.__setattr__(self, "labels", labels)
        object.__setattr__(self, "residuals", residuals)


@dataclass(frozen=True)
class OPEPrior:
    covariance: np.ndarray

    def __post_init__(self) -> None:
        covariance = _readonly_array(self.covariance, "covariance", shape=(5, 5), real=True)
        if not np.array_equal(covariance, covariance.T):
            raise ValueError("covariance must be exactly symmetric; it is not symmetrized automatically")
        try:
            linalg.cholesky(covariance, lower=True, check_finite=True)
        except linalg.LinAlgError as error:
            raise ValueError("covariance must be symmetric positive definite") from error
        object.__setattr__(self, "covariance", covariance)


def fixed_ope_prior_scales(parameters: OPEParameters) -> np.ndarray:
    """Return fixed structural scales for the five OPE residual conventions.

    The scales have the same dimensionless convention as the ``c_i,d_i``.
    They depend only on the caller-supplied reference parameters and are meant
    to be evaluated before, not during, a fit.
    """
    if not isinstance(parameters, OPEParameters):
        raise TypeError("parameters must be OPEParameters")
    common = parameters.n_c / (64.0 * np.pi**2)
    s_c = common * parameters.m_v / (np.sqrt(2.0) * parameters.f_v)
    d_g = common * parameters.m_v**2 / parameters.f_v**2
    d_f = parameters.f**2 / (8.0 * parameters.f_v**2)
    return _readonly_array(
        (s_c, s_c, s_c, d_g + 2.0 * d_f, d_g + d_f),
        "fixed OPE prior scales",
        shape=(5,),
        real=True,
    )


def fixed_ope_prior(
    parameters: OPEParameters,
    truncation_fraction: float | None = None,
    rho: float = 0.0,
) -> OPEPrior:
    """Build a fixed Gaussian OPE prior from reference parameters.

    Gaussianity, the structural scale prescription, and the baseline
    ``truncation_fraction=1/n_c`` are assumptions.  ``rho`` is a fixed
    equicorrelation and must lie in the strict SPD interval ``(-1/4, 1)``.
    The five residual definitions themselves are literature inputs.
    """
    if not isinstance(parameters, OPEParameters):
        raise TypeError("parameters must be OPEParameters")
    fraction = 1.0 / parameters.n_c if truncation_fraction is None else _positive_real(
        truncation_fraction, "truncation_fraction"
    )
    correlation = _finite_real(rho, "rho")
    if not -0.25 < correlation < 1.0:
        raise ValueError("rho must satisfy -1/4 < rho < 1 for positive definiteness")
    sigmas = fraction * fixed_ope_prior_scales(parameters)
    correlation_matrix = np.full((5, 5), correlation)
    np.fill_diagonal(correlation_matrix, 1.0)
    covariance = np.outer(sigmas, sigmas) * correlation_matrix
    return OPEPrior(covariance)


@dataclass(frozen=True)
class PhotonVectorTransform:
    epsilon: np.ndarray
    rotation: np.ndarray
    z: float
    k_sqrt: np.ndarray
    h: np.ndarray
    rho: np.ndarray
    s: np.ndarray

    def __post_init__(self) -> None:
        epsilon = _readonly_array(self.epsilon, "epsilon", ndim=1, real=True)
        n = epsilon.size
        rotation = _readonly_array(self.rotation, "rotation", shape=(n, n), real=True)
        if not np.allclose(rotation.T @ rotation, np.eye(n), rtol=0.0, atol=1e-12):
            raise ValueError("rotation must be orthogonal")
        for name, shape in (("k_sqrt", (n, n)), ("h", (n, n)), ("s", (n + 1, n + 1))):
            object.__setattr__(self, name, _readonly_array(getattr(self, name), name, shape=shape, real=True))
        object.__setattr__(self, "rho", _readonly_array(self.rho, "rho", shape=(n,), real=True))
        object.__setattr__(self, "epsilon", epsilon)
        object.__setattr__(self, "rotation", rotation)
        z = _positive_real(self.z, "z")
        expected_z = 1.0 + epsilon @ epsilon
        sqrt_z = np.sqrt(expected_z)
        expected_k = np.eye(n) - np.outer(epsilon, epsilon) / (sqrt_z * (sqrt_z + 1.0))
        expected_h = expected_k @ rotation
        expected_rho = rotation.T @ epsilon / sqrt_z
        expected_s = np.block([[np.array([[1.0 / sqrt_z]]), -expected_rho[None, :]],
                               [epsilon[:, None] / sqrt_z, expected_h]])
        consistent = (
            np.isclose(z, expected_z, rtol=1e-13, atol=1e-15)
            and np.allclose(self.k_sqrt, expected_k, rtol=1e-13, atol=1e-15)
            and np.allclose(self.h, expected_h, rtol=1e-13, atol=1e-15)
            and np.allclose(self.rho, expected_rho, rtol=1e-13, atol=1e-15)
            and np.allclose(self.s, expected_s, rtol=1e-13, atol=1e-15)
        )
        if not consistent:
            raise ValueError("z, k_sqrt, h, rho, and s must match epsilon and rotation")
        object.__setattr__(self, "z", z)


def odd_parity_flavor_kernels(source_a, source_b, spurions: ALPFlavorSpurions) -> OddParityFlavorKernels:
    """Return traces containing the explicit Alda:2026nrd physical-ALP spurions."""
    if not isinstance(spurions, ALPFlavorSpurions):
        raise TypeError("spurions must be ALPFlavorSpurions")
    ta = _source_matrix(source_a, "source_a")
    tb = _source_matrix(source_b, "source_b")
    sources = ta @ tb + tb @ ta
    c_q_AP_plus_m_a_inverse_C = spurions.c_q_AP + spurions.m_a_inverse_C
    mass = (
        c_q_AP_plus_m_a_inverse_C @ spurions.b0_mq
        + spurions.b0_mq @ c_q_AP_plus_m_a_inverse_C
    )
    return OddParityFlavorKernels(
        np.trace(sources @ spurions.m_a_inverse_C),
        np.trace(sources @ mass),
    )


def vjp_source_coefficients(kernels: OddParityFlavorKernels, e, f_a) -> OperatorCoefficients:
    if not isinstance(kernels, OddParityFlavorKernels):
        raise TypeError("kernels must be OddParityFlavorKernels")
    charge = _finite_real(e, "e")
    fa = _positive_real(f_a, "f_a")
    d, m = kernels.derivative_trace, kernels.mass_trace
    return OperatorCoefficients((4 * charge * d / fa, 4 * charge * d / fa, -8 * charge * m / fa, 0, 4 * charge * d / fa, 4 * charge * d / fa, 4 * charge * d / fa))


def weighted_vjp_coefficients(kernels, couplings: VJPCouplings, e, f_a, m_v) -> OperatorCoefficients:
    if not isinstance(couplings, VJPCouplings):
        raise TypeError("couplings must be VJPCouplings")
    source = vjp_source_coefficients(kernels, e, f_a).coefficients
    return OperatorCoefficients(source * couplings.as_array() / _positive_real(m_v, "m_v"))


def vvp_source_coefficients(kernels: OddParityFlavorKernels, f_a) -> OperatorCoefficients:
    if not isinstance(kernels, OddParityFlavorKernels):
        raise TypeError("kernels must be OddParityFlavorKernels")
    fa = _positive_real(f_a, "f_a")
    d, m = kernels.derivative_trace, kernels.mass_trace
    return OperatorCoefficients((-2 * d / fa, 4 * m / fa, -2 * d / fa, -2 * d / fa))


def weighted_vvp_coefficients(kernels, couplings: VVPCouplings, f_a) -> OperatorCoefficients:
    if not isinstance(couplings, VVPCouplings):
        raise TypeError("couplings must be VVPCouplings")
    return OperatorCoefficients(vvp_source_coefficients(kernels, f_a).coefficients * couplings.as_array())


def anomaly_seed(t_a, t_b, c_q_A, m_a_inverse_C, c_ab, n_c, identical_sources: bool) -> complex:
    if not isinstance(identical_sources, (bool, np.bool_)):
        raise TypeError("identical_sources must be an explicit boolean")
    if isinstance(n_c, bool) or not isinstance(n_c, (int, np.integer)) or n_c <= 0:
        raise ValueError("n_c must be a positive integer")
    ta, tb = _source_matrix(t_a, "t_a"), _source_matrix(t_b, "t_b")
    c_a = _hermitian_matrix(c_q_A, "c_q_A")
    m_a_inverse_C = _hermitian_matrix(m_a_inverse_C, "m_a_inverse_C")
    multiplicity = 1 if identical_sources else 2
    return _finite_scalar(c_ab, "c_ab") + multiplicity * int(n_c) * np.trace(
        (ta @ tb + tb @ ta) @ (c_a + m_a_inverse_C)
    )


def photon_anomaly_seed(q, c_q_A, m_a_inverse_C, c_gamma, n_c) -> complex:
    charge = _hermitian_matrix(q, "q")
    c_a = _hermitian_matrix(c_q_A, "c_q_A")
    m_a_inverse_C = _hermitian_matrix(m_a_inverse_C, "m_a_inverse_C")
    if isinstance(n_c, bool) or not isinstance(n_c, (int, np.integer)) or n_c <= 0:
        raise ValueError("n_c must be a positive integer")
    return _finite_scalar(c_gamma, "c_gamma") + 2 * int(n_c) * np.trace(
        charge @ charge @ (c_a + m_a_inverse_C)
    )


def ope_linear_system(parameters: OPEParameters) -> tuple[np.ndarray, np.ndarray]:
    """Return read-only ``A,b`` for ``A @ (c1..c7,d1..d4) = b``."""
    if not isinstance(parameters, OPEParameters):
        raise TypeError("parameters must be OPEParameters")
    a = np.zeros((5, 11))
    a[0, (0, 2)] = (1.0, 4.0)
    a[1, (0, 1, 4)] = (1.0, -1.0, 1.0)
    a[2, (4, 5)] = (1.0, -1.0)
    a[3, (7, 8)] = (1.0, 8.0)
    a[4, 9] = 1.0
    common = parameters.n_c / (64.0 * np.pi**2)
    b = np.array((0.0, 0.0, common * parameters.m_v / (np.sqrt(2.0) * parameters.f_v),
                  -common * parameters.m_v**2 / parameters.f_v**2 + parameters.f**2 / (4.0 * parameters.f_v**2),
                  -common * parameters.m_v**2 / parameters.f_v**2 + parameters.f**2 / (8.0 * parameters.f_v**2)))
    a.setflags(write=False)
    b.setflags(write=False)
    return a, b


def ope_residuals(couplings: OddParityCouplings, parameters: OPEParameters) -> OPEResult:
    if not isinstance(couplings, OddParityCouplings):
        raise TypeError("couplings must be OddParityCouplings")
    a, b = ope_linear_system(parameters)
    return OPEResult(OPE_LABELS, a @ couplings.as_array() - b)


def ope_prior_chi2(result: OPEResult, prior: OPEPrior) -> float:
    if not isinstance(result, OPEResult) or not isinstance(prior, OPEPrior):
        raise TypeError("result and prior must be OPEResult and OPEPrior")
    if result.residuals.shape != (5,):
        raise ValueError("the Ruiz--Femenia prior requires five residuals")
    factor = linalg.cholesky(prior.covariance, lower=True, check_finite=True)
    whitened = linalg.solve_triangular(factor, result.residuals, lower=True, check_finite=True)
    return float(whitened @ whitened)


def photon_vector_transform(epsilon, rotation) -> PhotonVectorTransform:
    eps = _readonly_array(epsilon, "epsilon", ndim=1, real=True)
    n = eps.size
    rot = _readonly_array(rotation, "rotation", shape=(n, n), real=True)
    if not np.allclose(rot.T @ rot, np.eye(n), rtol=0.0, atol=1e-12):
        raise ValueError("rotation must be orthogonal")
    z = float(1.0 + eps @ eps)
    sqrt_z = np.sqrt(z)
    k_sqrt = np.eye(n) - np.outer(eps, eps) / (sqrt_z * (sqrt_z + 1.0))
    h = k_sqrt @ rot
    rho = rot.T @ eps / sqrt_z
    s = np.block([[np.array([[1.0 / sqrt_z]]), -rho[None, :]], [eps[:, None] / sqrt_z, h]])
    return PhotonVectorTransform(eps, rot, z, k_sqrt, h, rho, s)


def transform_anomalous_bilinear(lambda_old, transform: PhotonVectorTransform) -> np.ndarray:
    """Apply ``S.T @ Lambda @ S``; no automatic ``C_AB`` conversion is made."""
    if not isinstance(transform, PhotonVectorTransform):
        raise TypeError("transform must be PhotonVectorTransform")
    matrix = _readonly_array(lambda_old, "lambda_old", shape=transform.s.shape, real=True)
    if not np.array_equal(matrix, matrix.T):
        raise ValueError("lambda_old must be exactly symmetric; it is not symmetrized automatically")
    result = transform.s.T @ matrix @ transform.s
    result.setflags(write=False)
    return result


def ideal_vector_generators() -> dict[str, np.ndarray]:
    """Return fresh generators in the repository's ideal omega--phi mixing convention."""
    def unit(row, column):
        value = np.zeros((3, 3), dtype=complex)
        value[row, column] = 1.0
        return value

    return {
        "rho+": unit(0, 1), "rho-": unit(1, 0),
        "rho0": np.diag([1.0, -1.0, 0.0]).astype(complex) / np.sqrt(2.0),
        "K*+": unit(0, 2), "K*-": unit(2, 0),
        "K*0": unit(1, 2), "K*0bar": unit(2, 1),
        "omega": np.diag([1.0, 1.0, 0.0]).astype(complex) / np.sqrt(2.0),
        "phi": np.diag([0.0, 0.0, 1.0]).astype(complex),
    }
