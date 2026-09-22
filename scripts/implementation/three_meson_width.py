"""Numerical leading-order widths for an ALP decaying to three mesons.

The coupling matrices are the hatted axial and pseudoscalar couplings in the
no-ALP-meson-mixing basis. Masses, ``f_a``, and the returned width are in GeV.

Example:
    .venv/bin/python -c "import numpy as np; from scripts.implementation import decay_width; print(decay_width(0.6, np.eye(3), np.eye(3), 'pi+', 'pi-', 'pi0', 1000.0))"
"""

from dataclasses import dataclass
from itertools import permutations

import numpy as np
from scipy.integrate import quad


F_PI = 0.093
"""Pion decay constant in GeV, in the convention of Alda:2026bar."""

# Appendix C fit values, converted from MeV^2 to GeV^2.  Keeping B_0 m_q as
# one object avoids introducing scale-dependent quark masses into the LO input.
B0_MQ = np.diag([0.006568, 0.011740, 0.235921])


@dataclass(frozen=True)
class _Meson:
    mass: float
    charge: int
    strangeness: int
    generator: np.ndarray


_SQRT2 = np.sqrt(2.0)
_ZERO = np.zeros((3, 3), dtype=complex)


def _off_diagonal(row: int, column: int) -> np.ndarray:
    generator = _ZERO.copy()
    generator[row, column] = _SQRT2
    return generator


# Phi and generator normalization follow Eqs. (III.1) and (C.6).  Observed
# masses from Appendix C Table I are used for physical phase space.
_MESONS = {
    "pi+": _Meson(0.13957039, 1, 0, _off_diagonal(0, 1)),
    "pi-": _Meson(0.13957039, -1, 0, _off_diagonal(1, 0)),
    "pi0": _Meson(
        0.1349768,
        0,
        0,
        np.diag([1.01707, -0.98257, -0.01140]).astype(complex),
    ),
    "K+": _Meson(0.493677, 1, 1, _off_diagonal(0, 2)),
    "K-": _Meson(0.493677, -1, -1, _off_diagonal(2, 0)),
    "K0": _Meson(0.497611, 0, 1, _off_diagonal(1, 2)),
    "K0bar": _Meson(0.497611, 0, -1, _off_diagonal(2, 1)),
    "eta": _Meson(
        0.547862,
        0,
        0,
        np.diag([0.79999, 0.83749, -0.81156]).astype(complex),
    ),
    "eta_prime": _Meson(
        0.95778,
        0,
        0,
        np.diag([0.57061, 0.57721, 1.15812]).astype(complex),
    ),
}

_ALIASES = {
    "pi+": "pi+",
    "pi-": "pi-",
    "pi0": "pi0",
    "pi_0": "pi0",
    "k+": "K+",
    "k-": "K-",
    "k0": "K0",
    "k0bar": "K0bar",
    "kbar0": "K0bar",
    "anti-k0": "K0bar",
    "eta": "eta",
    "eta'": "eta_prime",
    "eta_prime": "eta_prime",
    "etaprime": "eta_prime",
}


def _canonical_name(name: str) -> str:
    if not isinstance(name, str):
        raise TypeError("meson identifiers must be strings")
    try:
        return _ALIASES[name.strip().lower()]
    except KeyError as error:
        accepted = ", ".join(_MESONS)
        raise ValueError(f"unknown meson {name!r}; use one of: {accepted}") from error


def _coupling_matrix(value: np.ndarray, name: str) -> np.ndarray:
    matrix = np.asarray(value, dtype=complex)
    if matrix.shape != (3, 3):
        raise ValueError(f"{name} must have shape (3, 3)")
    if not np.all(np.isfinite(matrix)):
        raise ValueError(f"{name} must contain only finite values")
    if not np.allclose(matrix, matrix.conj().T, rtol=1e-12, atol=1e-12):
        raise ValueError(f"{name} must be Hermitian")
    return matrix


def _commutator(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    return left @ right - right @ left


def _kallen(x: float, y: float, z: float) -> float:
    return x * x + y * y + z * z - 2.0 * (x * y + x * z + y * z)


def _s23_width(s12: float, ma2: float, m1: float, m2: float, m3: float) -> float:
    product = _kallen(s12, m1 * m1, m2 * m2) * _kallen(
        ma2, s12, m3 * m3
    )
    # Roundoff can make the product slightly negative at a phase-space edge.
    return np.sqrt(max(0.0, product)) / s12


def _s23_bounds(
    s12: float, ma2: float, m1: float, m2: float, m3: float
) -> tuple[float, float]:
    center = m2 * m2 + m3 * m3 + (
        (s12 + m2 * m2 - m1 * m1) * (ma2 - s12 - m3 * m3)
    ) / (2.0 * s12)
    half_width = 0.5 * _s23_width(s12, ma2, m1, m2, m3)
    return center - half_width, center + half_width


def _trace(matrix: np.ndarray) -> complex:
    return complex(np.trace(matrix))


def _integrate_one_dimension(integrand, lower: float, upper: float) -> float:
    value, _ = quad(integrand, lower, upper, epsabs=0.0, epsrel=1e-9, limit=200)
    return float(value)


def decay_width(
    alp_mass: float,
    c_q_A: np.ndarray,
    c_q_P: np.ndarray,
    meson_1: str,
    meson_2: str,
    meson_3: str,
    f_a: float,
) -> float:
    """Return ``Gamma(a -> meson_1 meson_2 meson_3)`` in GeV.

    Supported canonical identifiers are ``pi+``, ``pi-``, ``pi0``, ``K+``,
    ``K-``, ``K0``, ``K0bar``, ``eta``, and ``eta_prime``.  ``eta'`` and a
    few common ASCII spellings are accepted as aliases.

    Charge-, ``|S| > 1``-, and threshold-forbidden channels return ``0.0``
    before the coupling matrices are inspected.  A channel with ``|S| = 1``
    is allowed only when at least one coupling has a nonzero d-s entry.
    """
    names = tuple(_canonical_name(name) for name in (meson_1, meson_2, meson_3))
    mesons = tuple(_MESONS[name] for name in names)

    if sum(meson.charge for meson in mesons) != 0:
        return 0.0
    total_strangeness = sum(meson.strangeness for meson in mesons)
    if abs(total_strangeness) > 1:
        return 0.0

    try:
        ma = float(alp_mass)
    except (TypeError, ValueError) as error:
        raise TypeError("alp_mass must be a real scalar") from error
    if not np.isfinite(ma):
        raise ValueError("alp_mass must be finite")
    if ma <= sum(meson.mass for meson in mesons):
        return 0.0

    c_a = None
    c_p = None
    if abs(total_strangeness) == 1:
        c_a = _coupling_matrix(c_q_A, "c_q_A")
        c_p = _coupling_matrix(c_q_P, "c_q_P")
        if c_a[1, 2] == 0.0 and c_p[1, 2] == 0.0:
            return 0.0

    try:
        fa = float(f_a)
    except (TypeError, ValueError) as error:
        raise TypeError("f_a must be a real scalar") from error
    if not np.isfinite(fa) or fa <= 0.0:
        raise ValueError("f_a must be finite and positive")

    if c_a is None:
        c_a = _coupling_matrix(c_q_A, "c_q_A")
        c_p = _coupling_matrix(c_q_P, "c_q_P")
    p_matrix = c_p @ B0_MQ + B0_MQ @ c_p
    kinetic_coefficient = 1.0 / (6.0 * fa * F_PI)
    breaking_coefficient = -1.0 / (6.0 * fa * F_PI)
    ma2 = ma * ma
    prefactor = 1.0 / (256.0 * np.pi**3 * ma**3)

    unique_names = set(names)
    if len(unique_names) == 1:
        meson = mesons[0]
        generator = meson.generator
        tensor_p = _trace(generator @ generator @ generator @ p_matrix)
        amplitude = 6.0 * breaking_coefficient * tensor_p
        lower = (2.0 * meson.mass) ** 2
        upper = (ma - meson.mass) ** 2
        integral = _integrate_one_dimension(
            lambda s12: _s23_width(
                s12, ma2, meson.mass, meson.mass, meson.mass
            )
            * abs(amplitude) ** 2,
            lower,
            upper,
        )
        return max(0.0, prefactor * integral / 6.0)

    if len(unique_names) == 2:
        identical_name = next(name for name in unique_names if names.count(name) == 2)
        distinct_name = next(name for name in unique_names if names.count(name) == 1)
        identical = _MESONS[identical_name]
        distinct = _MESONS[distinct_name]
        generator_i = identical.generator
        generator_j = distinct.generator
        tensor_a = _trace(
            c_a
            @ _commutator(
                generator_i, _commutator(generator_i, generator_j)
            )
        )
        tensor_p = _trace(
            (
                generator_i @ generator_i @ generator_j
                + generator_i @ generator_j @ generator_i
                + generator_j @ generator_i @ generator_i
            )
            @ p_matrix
        )

        def amplitude(s12: float) -> complex:
            return (
                0.5
                * kinetic_coefficient
                * tensor_a
                * (ma2 + 3.0 * distinct.mass**2 - 3.0 * s12)
                + 2.0 * breaking_coefficient * tensor_p
            )

        lower = (2.0 * identical.mass) ** 2
        upper = (ma - distinct.mass) ** 2
        integral = _integrate_one_dimension(
            lambda s12: _s23_width(
                s12,
                ma2,
                identical.mass,
                identical.mass,
                distinct.mass,
            )
            * abs(amplitude(s12)) ** 2,
            lower,
            upper,
        )
        return max(0.0, prefactor * integral / 2.0)

    generators = tuple(meson.generator for meson in mesons)
    masses = tuple(meson.mass for meson in mesons)
    generator_i, generator_j, generator_k = generators
    m_i, m_j, m_k = masses

    def tensor_a(first: np.ndarray, second: np.ndarray, differentiated: np.ndarray) -> complex:
        return _trace(
            c_a
            @ (
                _commutator(first, _commutator(second, differentiated))
                + _commutator(second, _commutator(first, differentiated))
            )
        )

    tensor_ij_k = tensor_a(generator_i, generator_j, generator_k)
    tensor_ik_j = tensor_a(generator_i, generator_k, generator_j)
    tensor_jk_i = tensor_a(generator_j, generator_k, generator_i)
    tensor_p = _trace(
        sum(
            (
                generators[order[0]]
                @ generators[order[1]]
                @ generators[order[2]]
                for order in permutations(range(3))
            ),
            _ZERO.copy(),
        )
        @ p_matrix
    )

    def amplitude(s12: float, s23: float) -> complex:
        s13 = ma2 + m_i**2 + m_j**2 + m_k**2 - s12 - s23
        return (
            0.5
            * kinetic_coefficient
            * (
                tensor_ij_k * (ma2 + m_k**2 - s12)
                + tensor_ik_j * (ma2 + m_j**2 - s13)
                + tensor_jk_i * (ma2 + m_i**2 - s23)
            )
            + breaking_coefficient * tensor_p
        )

    lower = (m_i + m_j) ** 2
    upper = (ma - m_k) ** 2

    def integrate_s23(s12: float) -> float:
        s23_min, s23_max = _s23_bounds(s12, ma2, m_i, m_j, m_k)
        return _integrate_one_dimension(
            lambda s23: abs(amplitude(s12, s23)) ** 2,
            s23_min,
            s23_max,
        )

    integral = _integrate_one_dimension(integrate_s23, lower, upper)
    return max(0.0, prefactor * integral)
