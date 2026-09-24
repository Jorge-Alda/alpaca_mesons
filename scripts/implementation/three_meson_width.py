"""Numerical leading-order amplitudes and widths for ALP three-meson decays.

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

EPSILON_K = 2.228e-3 * np.exp(1j * np.deg2rad(43.52))
"""Default indirect-CP-violation parameter from ParticleDataGroup:2024cfk."""

# Appendix C fit values, converted from MeV^2 to GeV^2.  Keeping B_0 m_q as
# one object avoids introducing scale-dependent quark masses into the LO input.
B0_MQ = np.diag([0.006568, 0.011740, 0.235921])


@dataclass(frozen=True)
class _Meson:
    mass: float
    charge: int
    strangeness: int | None
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

# Physical K_S and K_L are coherent flavor superpositions, including the
# configured epsilon_K, and therefore do not carry a definite strangeness.


def _kaon_mixing_coefficients(epsilon_k: complex) -> tuple[complex, complex]:
    """Return normalized ``(p, q)`` for the adopted epsilon_K convention."""
    try:
        epsilon = complex(epsilon_k)
    except (TypeError, ValueError) as error:
        raise TypeError("epsilon_k must be a complex scalar") from error
    if not np.isfinite(epsilon.real) or not np.isfinite(epsilon.imag):
        raise ValueError("epsilon_k must be finite")
    normalization = np.sqrt(2.0 * (1.0 + abs(epsilon) ** 2))
    return (1.0 + epsilon) / normalization, (1.0 - epsilon) / normalization


def _physical_neutral_kaon(name: str, epsilon_k: complex) -> _Meson:
    """Construct a final-state KS or KL generator, including epsilon_K."""
    p_mixing, q_mixing = _kaon_mixing_coefficients(epsilon_k)
    sign = 1.0 if name == "KS" else -1.0
    generator = (
        np.conjugate(p_mixing) * _MESONS["K0"].generator
        + sign * np.conjugate(q_mixing) * _MESONS["K0bar"].generator
    )
    return _Meson(0.497611, 0, None, generator)


_MESONS["KS"] = _physical_neutral_kaon("KS", EPSILON_K)
_MESONS["KL"] = _physical_neutral_kaon("KL", EPSILON_K)

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
    "ks": "KS",
    "k_s": "KS",
    "kl": "KL",
    "k_l": "KL",
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


@dataclass(frozen=True)
class ThreeMesonAmplitude:
    """Leading-order amplitude with invariants ordered like ``channel``.

    Calling the model as ``model(s12, s23)`` uses
    ``s12 = (p1 + p2)^2`` and ``s23 = (p2 + p3)^2``, where ``p1``, ``p2``,
    and ``p3`` follow the order of the three mesons in ``channel``.
    """

    alp_mass: float
    channel: tuple[str, str, str]
    masses: tuple[float, float, float]
    kinetic_coefficient: float
    breaking_coefficient: float
    kinetic_tensors: tuple[complex, ...]
    breaking_tensor: complex
    identical_pair: tuple[int, int] | None

    def __call__(self, s12: float, s23: float) -> complex:
        """Evaluate the invariant amplitude at the ordered Dalitz variables."""
        ma2 = self.alp_mass**2
        if len(set(self.channel)) == 1:
            return 6.0 * self.breaking_coefficient * self.breaking_tensor

        if self.identical_pair is not None:
            s13 = ma2 + sum(mass**2 for mass in self.masses) - s12 - s23
            pair_invariants = {(0, 1): s12, (0, 2): s13, (1, 2): s23}
            distinct_index = next(
                index for index in range(3) if index not in self.identical_pair
            )
            pair_invariant = pair_invariants[self.identical_pair]
            return (
                0.5
                * self.kinetic_coefficient
                * self.kinetic_tensors[0]
                * (
                    ma2
                    + 3.0 * self.masses[distinct_index] ** 2
                    - 3.0 * pair_invariant
                )
                + 2.0 * self.breaking_coefficient * self.breaking_tensor
            )

        m1, m2, m3 = self.masses
        s13 = ma2 + m1**2 + m2**2 + m3**2 - s12 - s23
        tensor_12_3, tensor_13_2, tensor_23_1 = self.kinetic_tensors
        return (
            0.5
            * self.kinetic_coefficient
            * (
                tensor_12_3 * (ma2 + m3**2 - s12)
                + tensor_13_2 * (ma2 + m2**2 - s13)
                + tensor_23_1 * (ma2 + m1**2 - s23)
            )
            + self.breaking_coefficient * self.breaking_tensor
        )

    def _at_identical_pair_invariant(self, invariant: float) -> complex:
        """Evaluate a repeated-particle amplitude using its sole invariant."""
        if self.identical_pair == (0, 1):
            return self(invariant, 0.0)
        if self.identical_pair == (1, 2):
            return self(0.0, invariant)
        if self.identical_pair == (0, 2):
            invariant_sum = self.alp_mass**2 + sum(
                mass**2 for mass in self.masses
            )
            return self(0.0, invariant_sum - invariant)
        raise ValueError("channel does not contain exactly two identical mesons")


def _channel_mesons(
    names: tuple[str, str, str], epsilon_k: complex
) -> tuple[_Meson, _Meson, _Meson]:
    return tuple(
        _physical_neutral_kaon(name, epsilon_k)
        if name in {"KS", "KL"}
        else _MESONS[name]
        for name in names
    )


def _build_amplitude(
    alp_mass: float,
    f_a: float,
    names: tuple[str, str, str],
    mesons: tuple[_Meson, _Meson, _Meson],
    c_a: np.ndarray,
    c_p: np.ndarray,
) -> ThreeMesonAmplitude:
    """Build an amplitude from validated numerical and channel inputs."""
    p_matrix = c_p @ B0_MQ + B0_MQ @ c_p
    kinetic_coefficient = 1.0 / (6.0 * f_a * F_PI)
    breaking_coefficient = -1.0 / (6.0 * f_a * F_PI)
    generators = tuple(meson.generator for meson in mesons)
    masses = tuple(meson.mass for meson in mesons)
    unique_names = set(names)

    if len(unique_names) == 1:
        tensor_p = _trace(generators[0] @ generators[0] @ generators[0] @ p_matrix)
        return ThreeMesonAmplitude(
            alp_mass,
            names,
            masses,
            kinetic_coefficient,
            breaking_coefficient,
            (),
            tensor_p,
            None,
        )

    if len(unique_names) == 2:
        identical_name = next(name for name in unique_names if names.count(name) == 2)
        identical_indices = tuple(
            index for index, name in enumerate(names) if name == identical_name
        )
        distinct_index = next(
            index for index in range(3) if index not in identical_indices
        )
        generator_i = generators[identical_indices[0]]
        generator_j = generators[distinct_index]
        tensor_a = _trace(
            c_a
            @ _commutator(generator_i, _commutator(generator_i, generator_j))
        )
        tensor_p = _trace(
            (
                generator_i @ generator_i @ generator_j
                + generator_i @ generator_j @ generator_i
                + generator_j @ generator_i @ generator_i
            )
            @ p_matrix
        )
        return ThreeMesonAmplitude(
            alp_mass,
            names,
            masses,
            kinetic_coefficient,
            breaking_coefficient,
            (tensor_a,),
            tensor_p,
            identical_indices,
        )

    def tensor_a(
        first: np.ndarray, second: np.ndarray, differentiated: np.ndarray
    ) -> complex:
        return _trace(
            c_a
            @ (
                _commutator(first, _commutator(second, differentiated))
                + _commutator(second, _commutator(first, differentiated))
            )
        )

    generator_1, generator_2, generator_3 = generators
    kinetic_tensors = (
        tensor_a(generator_1, generator_2, generator_3),
        tensor_a(generator_1, generator_3, generator_2),
        tensor_a(generator_2, generator_3, generator_1),
    )
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
    return ThreeMesonAmplitude(
        alp_mass,
        names,
        masses,
        kinetic_coefficient,
        breaking_coefficient,
        kinetic_tensors,
        tensor_p,
        None,
    )


def build_amplitude(
    alp_mass: float,
    c_q_A: np.ndarray,
    c_q_P: np.ndarray,
    meson_1: str,
    meson_2: str,
    meson_3: str,
    f_a: float,
    *,
    epsilon_k: complex = EPSILON_K,
) -> ThreeMesonAmplitude:
    """Return the LO amplitude model for an ordered three-meson channel."""
    try:
        ma = float(alp_mass)
        fa = float(f_a)
    except (TypeError, ValueError) as error:
        raise TypeError("alp_mass and f_a must be real scalars") from error
    if not np.isfinite(ma) or ma <= 0.0:
        raise ValueError("alp_mass must be finite and positive")
    if not np.isfinite(fa) or fa <= 0.0:
        raise ValueError("f_a must be finite and positive")

    names = tuple(_canonical_name(name) for name in (meson_1, meson_2, meson_3))
    mesons = _channel_mesons(names, epsilon_k)
    c_a = _coupling_matrix(c_q_A, "c_q_A")
    c_p = _coupling_matrix(c_q_P, "c_q_P")
    return _build_amplitude(ma, fa, names, mesons, c_a, c_p)


def decay_width(
    alp_mass: float,
    c_q_A: np.ndarray,
    c_q_P: np.ndarray,
    meson_1: str,
    meson_2: str,
    meson_3: str,
    f_a: float,
    *,
    epsilon_k: complex = EPSILON_K,
) -> float:
    """Return ``Gamma(a -> meson_1 meson_2 meson_3)`` in GeV.

    Supported canonical identifiers are ``pi+``, ``pi-``, ``pi0``, ``K+``,
    ``K-``, ``K0``, ``K0bar``, ``KS``, ``KL``, ``eta``, and ``eta_prime``.
    The physical neutral states use ``KS=p*K0+q*K0bar`` and
    ``KL=p*K0-q*K0bar``, with ``p`` and ``q`` constructed from ``epsilon_k``.
    The default is the phenomenological PDG value; pass ``epsilon_k=0`` for
    the CP-conserving limit. ``eta'`` and common ASCII spellings are aliases.

    Definite-flavor channels forbidden by charge, ``|S| > 1``, or threshold
    return ``0.0`` before the coupling matrices are inspected. A definite
    flavor channel with ``|S| = 1`` is allowed only when at least one coupling
    has a nonzero d-s entry. Channels containing ``KS`` or ``KL`` have no
    definite strangeness; their allowed components and interference are
    selected directly by the flavor traces.
    """
    names = tuple(_canonical_name(name) for name in (meson_1, meson_2, meson_3))
    has_physical_neutral_kaon = any(name in {"KS", "KL"} for name in names)
    mesons = _channel_mesons(names, epsilon_k)

    if sum(meson.charge for meson in mesons) != 0:
        return 0.0
    total_strangeness = None
    if not has_physical_neutral_kaon:
        total_strangeness = sum(int(meson.strangeness) for meson in mesons)
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
    if total_strangeness is not None and abs(total_strangeness) == 1:
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
    amplitude = _build_amplitude(ma, fa, names, mesons, c_a, c_p)
    ma2 = ma * ma
    prefactor = 1.0 / (256.0 * np.pi**3 * ma**3)

    unique_names = set(names)
    if len(unique_names) == 1:
        meson = mesons[0]
        lower = (2.0 * meson.mass) ** 2
        upper = (ma - meson.mass) ** 2
        integral = _integrate_one_dimension(
            lambda s12: _s23_width(
                s12, ma2, meson.mass, meson.mass, meson.mass
            )
            * abs(amplitude(0.0, 0.0)) ** 2,
            lower,
            upper,
        )
        return max(0.0, prefactor * integral / 6.0)

    if len(unique_names) == 2:
        identical_name = next(name for name in unique_names if names.count(name) == 2)
        distinct_name = next(name for name in unique_names if names.count(name) == 1)
        identical = mesons[names.index(identical_name)]
        distinct = mesons[names.index(distinct_name)]

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
            * abs(amplitude._at_identical_pair_invariant(s12)) ** 2,
            lower,
            upper,
        )
        return max(0.0, prefactor * integral / 2.0)

    masses = tuple(meson.mass for meson in mesons)
    m_i, m_j, m_k = masses

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
