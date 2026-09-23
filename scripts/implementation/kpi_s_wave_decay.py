"""Charge-resolved ``a -> K Kbar pi`` amplitudes with K-pi S waves.

This module implements the independent-band matching discussed in
``memory/kpi-s-wave-chiral-matching.md``.  Despite the short name, the
implemented final state is ``K Kbar pi`` (not a net-strangeness ``K pi pi``
state).  Masses, ``f_a``, Mandelstam invariants, and returned widths use GeV.

The optional ``fsi`` callable has signature ``fsi(s, pair) -> complex``, where
``pair`` is a canonical two-string tuple such as ``("K+", "pi0")``.  It must
be normalized consistently with the low-energy matching.  ``fsi=None`` means
unity and exactly restores the full leading-order contact amplitude.

Example:
    .venv/bin/python -c "import numpy as np; from scripts.implementation.kpi_s_wave_decay import s_wave_decay_width; print(s_wave_decay_width(1.5, np.diag([1., 0., 0.]), np.zeros((3, 3)), ('K+', 'K-', 'pi0'), 1000.))"

This is a charge-resolved, diagonal-FSI independent-isobar approximation.  It
does not solve the coupled charge-basis Omnes problem, include dynamical
``m_d-m_u`` corrections to K-pi scattering, electromagnetic corrections, or
crossed-channel three-body rescattering.
"""

from dataclasses import dataclass
from itertools import permutations
from typing import Callable

import numpy as np
from scipy.integrate import quad

# Physical generators, masses, f_pi, and B_0 m_q use the fit conventions and
# inputs of Alda:2026bar at the precise point where they enter this module.
from .three_meson_width import B0_MQ, F_PI, _MESONS, _canonical_name


ScalarFSI = Callable[[float, tuple[str, str]], complex]


@dataclass(frozen=True)
class AffineKernel:
    """A single-band production kernel ``P(s) = slope*s + intercept``."""

    slope: complex
    intercept: complex
    pair: tuple[str, str]

    def __call__(self, invariant: float) -> complex:
        return self.slope * invariant + self.intercept


@dataclass(frozen=True)
class SWaveModel:
    """Kinematics and the two matched K-pi production kernels for one mode."""

    alp_mass: float
    f_a: float
    channel: tuple[str, str, str]
    masses: tuple[float, float, float]
    s_kernel: AffineKernel
    t_kernel: AffineKernel

    def contact_amplitude(self, s: float, t: float) -> complex:
        """Return the undressed leading-order contact amplitude."""
        return self.s_kernel(s) + self.t_kernel(t)

    def amplitude(self, s: float, t: float, fsi: ScalarFSI | None = None) -> complex:
        """Return ``P_s(s) F_s(s) + P_t(t) F_t(t)``."""
        if fsi is None:
            return self.contact_amplitude(s, t)
        # Omnes:1958hv motivates multiplying each matched kernel by a
        # normalized analytic FSI function; see implementation.omnes.
        return self.s_kernel(s) * complex(fsi(s, self.s_kernel.pair)) + self.t_kernel(
            t
        ) * complex(fsi(t, self.t_kernel.pair))


# The ordering is (kaon, pion, antikaon), so s=(p_K+p_pi)^2 and
# t=(p_pi+p_Kbar)^2.  These are all charge-conserving K Kbar pi flavor modes.
_CHANNELS = {
    frozenset(("K+", "pi0", "K-")): ("K+", "pi0", "K-"),
    frozenset(("K0", "pi0", "K0bar")): ("K0", "pi0", "K0bar"),
    frozenset(("K+", "pi-", "K0bar")): ("K+", "pi-", "K0bar"),
    frozenset(("K0", "pi+", "K-")): ("K0", "pi+", "K-"),
}


def supported_channels() -> tuple[tuple[str, str, str], ...]:
    """Return canonical supported ``(K, pi, Kbar)`` channels."""
    return tuple(_CHANNELS.values())


def kpi_thresholds() -> dict[tuple[str, str], float]:
    """Return physical K-pi thresholds in GeV for all implemented bands."""
    pairs = {
        (channel[0], channel[1]) for channel in _CHANNELS.values()
    } | {(channel[2], channel[1]) for channel in _CHANNELS.values()}
    # Alda:2026bar supplies the charged/neutral physical masses used here.
    return {pair: _MESONS[pair[0]].mass + _MESONS[pair[1]].mass for pair in pairs}


def _trace(matrix: np.ndarray) -> complex:
    return complex(np.trace(matrix))


def _commutator(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    return left @ right - right @ left


def _diagonal_coupling(value: np.ndarray, name: str) -> np.ndarray:
    matrix = np.asarray(value, dtype=complex)
    if matrix.shape != (3, 3):
        raise ValueError(f"{name} must have shape (3, 3)")
    if not np.all(np.isfinite(matrix)):
        raise ValueError(f"{name} must contain only finite values")
    if not np.allclose(matrix, matrix.conj().T, rtol=1e-12, atol=1e-12):
        raise ValueError(f"{name} must be Hermitian")
    if not np.allclose(matrix, np.diag(np.diag(matrix)), rtol=0.0, atol=1e-14):
        raise ValueError(f"{name} must be diagonal for the neutral ALP source")
    return matrix


def _canonical_channel(channel: tuple[str, str, str]) -> tuple[str, str, str]:
    if len(channel) != 3:
        raise ValueError("channel must contain exactly three meson identifiers")
    names = tuple(_canonical_name(name) for name in channel)
    try:
        return _CHANNELS[frozenset(names)]
    except KeyError as error:
        choices = ", ".join(" ".join(item) for item in _CHANNELS.values())
        raise ValueError(f"unsupported K Kbar pi channel; use one of: {choices}") from error


def _kinetic_tensor(
    c_a: np.ndarray,
    first: np.ndarray,
    second: np.ndarray,
    differentiated: np.ndarray,
) -> complex:
    return _trace(
        c_a
        @ (
            _commutator(first, _commutator(second, differentiated))
            + _commutator(second, _commutator(first, differentiated))
        )
    )


def _kernels_from_inputs(
    alp_mass: float,
    f_a: float,
    names: tuple[str, str, str],
    masses: tuple[float, float, float],
    generators: tuple[np.ndarray, np.ndarray, np.ndarray],
    c_a: np.ndarray,
    c_p: np.ndarray,
    b0_mq: np.ndarray,
) -> tuple[AffineKernel, AffineKernel]:
    """Construct kernels; explicit inputs also permit algebraic unit tests."""
    g_k, g_pi, g_kbar = generators
    m_k, m_pi, m_kbar = masses
    tensor_s = _kinetic_tensor(c_a, g_k, g_pi, g_kbar)
    tensor_u = _kinetic_tensor(c_a, g_k, g_kbar, g_pi)
    tensor_t = _kinetic_tensor(c_a, g_pi, g_kbar, g_k)
    symmetry_breaking = c_p @ b0_mq + b0_mq @ c_p
    ordered_sum = sum(
        (order[0] @ order[1] @ order[2] for order in permutations(generators)),
        np.zeros((3, 3), dtype=complex),
    )
    tensor_p = _trace(ordered_sum @ symmetry_breaking)

    # Alda:2026bar: expansion of the chiral-invariant LO a Phi^3 interaction.
    # Eliminating u leaves an affine function of the two K-pi invariants.
    normalization = 1.0 / (12.0 * f_a * F_PI)
    slope_s = normalization * (tensor_u - tensor_s)
    slope_t = normalization * (tensor_u - tensor_t)
    constant = normalization * (
        tensor_s * (alp_mass**2 + m_kbar**2)
        - tensor_u * (m_k**2 + m_kbar**2)
        + tensor_t * (alp_mass**2 + m_k**2)
        - 2.0 * tensor_p
    )

    # Assumption: split the invariant-independent contact term equally between
    # the two bands.  This is the unique split invariant under band exchange;
    # only their sum is fixed for non-charge-conjugate bands.
    return (
        AffineKernel(slope_s, 0.5 * constant, (names[0], names[1])),
        AffineKernel(slope_t, 0.5 * constant, (names[2], names[1])),
    )


def build_s_wave_model(
    alp_mass: float,
    c_q_A: np.ndarray,
    c_q_P: np.ndarray,
    channel: tuple[str, str, str],
    f_a: float,
) -> SWaveModel:
    """Build a charge-resolved two-band model for one physical flavor mode.

    Couplings are diagonal hatted couplings in the no-mixing basis.  Physical
    neutral generators and nondegenerate ``B_0 m_q`` values retain source-level
    isospin breaking.  The scattering FSI itself is supplied separately.
    """
    try:
        ma = float(alp_mass)
        fa = float(f_a)
    except (TypeError, ValueError) as error:
        raise TypeError("alp_mass and f_a must be real scalars") from error
    if not np.isfinite(ma) or ma <= 0.0:
        raise ValueError("alp_mass must be finite and positive")
    if not np.isfinite(fa) or fa <= 0.0:
        raise ValueError("f_a must be finite and positive")
    names = _canonical_channel(channel)
    c_a = _diagonal_coupling(c_q_A, "c_q_A")
    c_p = _diagonal_coupling(c_q_P, "c_q_P")
    mesons = tuple(_MESONS[name] for name in names)
    masses = tuple(meson.mass for meson in mesons)
    generators = tuple(meson.generator for meson in mesons)
    kernels = _kernels_from_inputs(
        ma, fa, names, masses, generators, c_a, c_p, B0_MQ
    )
    return SWaveModel(ma, fa, names, masses, *kernels)


def _kallen(x: float, y: float, z: float) -> float:
    return x * x + y * y + z * z - 2.0 * (x * y + x * z + y * z)


def _t_bounds(s: float, ma2: float, m_k: float, m_pi: float, m_kbar: float) -> tuple[float, float]:
    center = m_pi**2 + m_kbar**2 + (
        (s + m_pi**2 - m_k**2) * (ma2 - s - m_kbar**2)
    ) / (2.0 * s)
    product = _kallen(s, m_k**2, m_pi**2) * _kallen(ma2, s, m_kbar**2)
    half_width = 0.5 * np.sqrt(max(0.0, product)) / s
    return center - half_width, center + half_width


def s_wave_decay_width(
    alp_mass: float,
    c_q_A: np.ndarray,
    c_q_P: np.ndarray,
    channel: tuple[str, str, str],
    f_a: float,
    fsi: ScalarFSI | None = None,
) -> float:
    """Return the independent-band ``a -> K Kbar pi`` width in GeV."""
    model = build_s_wave_model(alp_mass, c_q_A, c_q_P, channel, f_a)
    m_k, m_pi, m_kbar = model.masses
    if model.alp_mass <= sum(model.masses):
        return 0.0
    ma2 = model.alp_mass**2
    s_min = (m_k + m_pi) ** 2
    s_max = (model.alp_mass - m_kbar) ** 2

    def integrate_t(s: float) -> float:
        t_min, t_max = _t_bounds(s, ma2, m_k, m_pi, m_kbar)
        value, _ = quad(
            lambda t: abs(model.amplitude(s, t, fsi)) ** 2,
            t_min,
            t_max,
            epsabs=0.0,
            epsrel=1e-8,
            limit=200,
        )
        return float(value)

    integral, _ = quad(
        integrate_t, s_min, s_max, epsabs=0.0, epsrel=1e-8, limit=200
    )
    return max(0.0, integral / (256.0 * np.pi**3 * model.alp_mass**3))
