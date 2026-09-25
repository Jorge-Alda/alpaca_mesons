"""Canonical-vector pole and tensor-to-Proca matching contributions.

The ordered arguments and Dalitz variables match :func:`build_amplitude` in
``three_meson_width``.  The pole and local matching terms are exposed as
separate builders so callers can use the standalone Proca pole or assemble the
representation-equivalent pole-plus-contact result explicitly.

``c_q_P`` is validated for interface compatibility but does not enter this
leading ``G_V`` exchange contribution.

Example:
    .venv/bin/python -c "import numpy as np; from scripts.implementation import build_vector_exchange_amplitude; m=build_vector_exchange_amplitude(0.8,np.diag([1.,0.,0.]),np.eye(3),'pi+','pi-','pi0',1000.); print(m(0.2,0.2), m.exchange_names)"
"""

from dataclasses import dataclass
from typing import Mapping

import numpy as np

from three_meson_width import (
    EPSILON_K,
    F_PI,
    _canonical_name,
    _channel_mesons,
    _commutator,
    _coupling_matrix,
    _trace,
)


@dataclass(frozen=True)
class VectorPoleParameters:
    """Constant-width Breit--Wigner inputs in GeV."""

    mass: float
    width: float


# PDG 2024 listing central values.  Charged/neutral category values are kept
# separate as documented in memory/extracts/ParticleDataGroup:2024cfk/.
DEFAULT_VECTOR_POLES: dict[str, VectorPoleParameters] = {
    "rho+": VectorPoleParameters(0.77511, 0.1491),
    "rho-": VectorPoleParameters(0.77511, 0.1491),
    "rho0": VectorPoleParameters(0.77526, 0.1474),
    "K*+": VectorPoleParameters(0.89167, 0.0514),
    "K*-": VectorPoleParameters(0.89167, 0.0514),
    "K*0": VectorPoleParameters(0.89555, 0.0473),
    "K*0bar": VectorPoleParameters(0.89555, 0.0473),
    "omega": VectorPoleParameters(0.78266, 0.00868),
    "phi": VectorPoleParameters(1.019461, 0.004249),
}


def _matrix_unit(row: int, column: int) -> np.ndarray:
    result = np.zeros((3, 3), dtype=complex)
    result[row, column] = 1.0
    return result


# Ideal omega--phi mixing is an explicit model assumption.  Unlike the Alda
# pseudoscalar generators imported above, the physical-vector matrix has unit
# off-diagonal entries.
_VECTOR_GENERATORS = {
    "rho+": _matrix_unit(0, 1),
    "rho-": _matrix_unit(1, 0),
    "rho0": np.diag([1.0, -1.0, 0.0]).astype(complex) / np.sqrt(2.0),
    "K*+": _matrix_unit(0, 2),
    "K*-": _matrix_unit(2, 0),
    "K*0": _matrix_unit(1, 2),
    "K*0bar": _matrix_unit(2, 1),
    "omega": np.diag([1.0, 1.0, 0.0]).astype(complex) / np.sqrt(2.0),
    "phi": np.diag([0.0, 0.0, 1.0]).astype(complex),
}

_DIAGONAL_VECTORS = frozenset(("rho0", "omega", "phi"))


def _roundoff_zero(value: complex, scale: float) -> bool:
    """Test a computed trace against its scale, without a fixed absolute cut."""
    return abs(value) <= 64.0 * np.finfo(float).eps * scale


def _is_diagonal(matrix: np.ndarray) -> bool:
    off_diagonal = matrix - np.diag(np.diag(matrix))
    scale = np.linalg.norm(matrix)
    return np.linalg.norm(off_diagonal) <= 64.0 * np.finfo(float).eps * scale


def _is_proportional_to_identity(matrix: np.ndarray) -> bool:
    """Return whether a validated matrix is identity-like at roundoff level."""
    singlet = np.trace(matrix) * np.eye(3, dtype=complex) / 3.0
    scale = np.linalg.norm(matrix)
    return np.linalg.norm(matrix - singlet) <= 64.0 * np.finfo(float).eps * scale


def _validated_poles(
    overrides: Mapping[str, VectorPoleParameters | tuple[float, float] | Mapping[str, float]] | None,
) -> dict[str, VectorPoleParameters]:
    poles = dict(DEFAULT_VECTOR_POLES)
    if overrides is None:
        return poles
    if not isinstance(overrides, Mapping):
        raise TypeError("pole_parameters must be a mapping keyed by exchange name")
    for name, value in overrides.items():
        if name not in poles:
            raise ValueError(f"unknown vector pole {name!r}")
        if isinstance(value, VectorPoleParameters):
            pole = value
        elif isinstance(value, Mapping):
            try:
                pole = VectorPoleParameters(float(value["mass"]), float(value["width"]))
            except KeyError as error:
                raise ValueError("pole mappings require mass and width") from error
        else:
            try:
                mass, width = value
                pole = VectorPoleParameters(float(mass), float(width))
            except (TypeError, ValueError) as error:
                raise TypeError("each pole override must provide (mass, width)") from error
        if not np.isfinite(pole.mass) or pole.mass <= 0.0:
            raise ValueError(f"{name} mass must be finite and positive")
        if not np.isfinite(pole.width) or pole.width < 0.0:
            raise ValueError(f"{name} width must be finite and nonnegative")
        poles[name] = pole
    return poles


@dataclass(frozen=True)
class VectorExchangeContribution:
    """One labeled spectator graph retained by the flavor traces."""

    exchange_name: str
    spectator_index: int
    pair_indices: tuple[int, int]
    vpp_trace: complex
    alp_trace: complex
    flavor_factor: complex
    pole: VectorPoleParameters


@dataclass(frozen=True)
class VectorExchangeAmplitude:
    """Callable vector-pole amplitude with ordered Dalitz variables.

    Calling ``model(s12, s23)`` means exactly ``s12=(p1+p2)^2`` and
    ``s23=(p2+p3)^2``, with particle labels following ``channel``.
    """

    alp_mass: float
    channel: tuple[str, str, str]
    masses: tuple[float, float, float]
    f_a: float
    f_pi: float
    g_v: float
    contributions: tuple[VectorExchangeContribution, ...]

    @property
    def exchange_names(self) -> tuple[str, ...]:
        """Exchange name for each retained labeled spectator graph."""
        return tuple(item.exchange_name for item in self.contributions)

    def __call__(self, s12: float, s23: float) -> complex:
        """Evaluate the sum of canonical-vector pole graphs."""
        pair_invariants = _ordered_pair_invariants(
            s12, s23, self.alp_mass, self.masses
        )
        result = 0.0j
        for item in self.contributions:
            i, j = item.pair_indices
            k = item.spectator_index
            s_ij = pair_invariants[(i, j)]
            s_ik = pair_invariants[tuple(sorted((i, k)))]
            s_jk = pair_invariants[tuple(sorted((j, k)))]
            mi2, mj2, mk2 = self.masses[i] ** 2, self.masses[j] ** 2, self.masses[k] ** 2

            # J_a.J_PP derived in derive_a3p_vector_exchange_amplitude.py.
            numerator = 0.25 * s_ij * (
                (s_ij + mj2 - mi2) * (s_ik - mi2 - mk2)
                - (s_ij + mi2 - mj2) * (s_jk - mj2 - mk2)
            )
            mass = item.pole.mass
            denominator = mass**2 - s_ij - 1j * mass * item.pole.width
            coefficient = self.g_v**2 / (self.f_pi**3 * self.f_a * mass**2)
            result += coefficient * item.flavor_factor * numerator / denominator
        return complex(result)


def _ordered_pair_invariants(
    s12: float,
    s23: float,
    alp_mass: float,
    masses: tuple[float, float, float],
) -> dict[tuple[int, int], float]:
    """Validate ordered Dalitz inputs and return all three pair invariants."""
    try:
        s12_value = float(s12)
        s23_value = float(s23)
    except (TypeError, ValueError) as error:
        raise TypeError("s12 and s23 must be real scalars") from error
    if not np.isfinite(s12_value) or not np.isfinite(s23_value):
        raise ValueError("s12 and s23 must be finite")

    invariant_sum = alp_mass**2 + sum(mass**2 for mass in masses)
    return {
        (0, 1): s12_value,
        (1, 2): s23_value,
        (0, 2): invariant_sum - s12_value - s23_value,
    }


@dataclass(frozen=True)
class VectorMatchingAmplitude:
    """Callable local tensor-to-Proca matching amplitude.

    Calling ``model(s12, s23)`` means exactly ``s12=(p1+p2)^2`` and
    ``s23=(p2+p3)^2``, with particle labels following ``channel``.  Pole
    widths are retained only as contribution metadata and never enter the
    local amplitude.
    """

    alp_mass: float
    channel: tuple[str, str, str]
    masses: tuple[float, float, float]
    f_a: float
    f_pi: float
    g_v: float
    contributions: tuple[VectorExchangeContribution, ...]

    @property
    def exchange_names(self) -> tuple[str, ...]:
        """Exchange name for each retained labeled matching contribution."""
        return tuple(item.exchange_name for item in self.contributions)

    def __call__(self, s12: float, s23: float) -> complex:
        """Evaluate the sum of local representation-matching terms."""
        pair_invariants = _ordered_pair_invariants(
            s12, s23, self.alp_mass, self.masses
        )
        ma2 = self.alp_mass**2
        result = 0.0j
        for item in self.contributions:
            i, j = item.pair_indices
            k = item.spectator_index
            s_ik = pair_invariants[tuple(sorted((i, k)))]
            s_jk = pair_invariants[tuple(sorted((j, k)))]
            mi2 = self.masses[i] ** 2
            mj2 = self.masses[j] ** 2
            mk2 = self.masses[k] ** 2

            # L_(k;ij), derived in derive_a3p_vector_exchange_amplitude.py.
            local_kinematic = 0.25 * (
                (s_ik - mi2 - mk2) * (ma2 + mj2 - s_ik)
                - (s_jk - mj2 - mk2) * (ma2 + mi2 - s_jk)
            )
            mass = item.pole.mass
            coefficient = self.g_v**2 / (self.f_pi**3 * self.f_a * mass**2)
            result += coefficient * item.flavor_factor * local_kinematic
        return complex(result)


def build_vector_exchange_amplitude(
    alp_mass: float,
    c_q_A: np.ndarray,
    c_q_P: np.ndarray,
    meson_1: str,
    meson_2: str,
    meson_3: str,
    f_a: float,
    *,
    epsilon_k: complex = EPSILON_K,
    g_v: float | None = None,
    pole_parameters: Mapping[
        str, VectorPoleParameters | tuple[float, float] | Mapping[str, float]
    ]
    | None = None,
    f_pi: float = F_PI,
) -> VectorExchangeAmplitude:
    """Build the leading ``G_V`` canonical-vector pole amplitude.

    Both hatted coupling arrays are validated as finite Hermitian matrices.
    ``c_q_P`` is intentionally absent from the result: it is accepted only to
    preserve the main positional interface of ``three_meson_width``.

    The default ``G_V=f_pi/sqrt(2)`` is the single-resonance short-distance
    assumption.  Pole masses and widths may be overridden by exchange name.
    Widths are constant, which is an explicit line-shape model assumption.
    """
    try:
        ma = float(alp_mass)
        fa = float(f_a)
        fpi = float(f_pi)
        gv = fpi / np.sqrt(2.0) if g_v is None else float(g_v)
    except (TypeError, ValueError) as error:
        raise TypeError("alp_mass, f_a, f_pi, and g_v must be real scalars") from error
    for value, name in ((ma, "alp_mass"), (fa, "f_a"), (fpi, "f_pi")):
        if not np.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be finite and positive")
    if not np.isfinite(gv) or gv < 0.0:
        raise ValueError("g_v must be finite and nonnegative")

    names = tuple(_canonical_name(name) for name in (meson_1, meson_2, meson_3))
    mesons = _channel_mesons(names, epsilon_k)
    c_a = _coupling_matrix(c_q_A, "c_q_A")
    _coupling_matrix(c_q_P, "c_q_P")  # validated, absent at leading G_V order
    poles = _validated_poles(pole_parameters)
    generators = tuple(meson.generator for meson in mesons)
    diagonal_c_a = _is_diagonal(c_a)

    retained: list[VectorExchangeContribution] = []
    # Every ALP trace contains [Lambda_k,c_q_A], so a singlet spurion makes the
    # complete topology vanish before any exchange-specific work is needed.
    if _is_proportional_to_identity(c_a):
        return VectorExchangeAmplitude(
            alp_mass=ma,
            channel=names,
            masses=tuple(meson.mass for meson in mesons),
            f_a=fa,
            f_pi=fpi,
            g_v=gv,
            contributions=(),
        )

    for spectator_index in range(3):
        pair_indices = tuple(index for index in range(3) if index != spectator_index)
        i, j = pair_indices  # naturally ordered i<j; no compressed multiplicities
        pair_commutator = _commutator(generators[i], generators[j])
        for exchange_name, vector_generator in _VECTOR_GENERATORS.items():
            # A diagonal vector is forbidden at the ALP vertex for diagonal
            # c_q_A, independently of the pair, so skip it before any trace.
            if diagonal_c_a and exchange_name in _DIAGONAL_VECTORS:
                continue
            # For the remaining candidates, reject a VPP zero before evaluating
            # the ALP trace or any kinematics.
            conjugate_generator = vector_generator.conj().T
            vpp_trace = _trace(conjugate_generator @ pair_commutator)
            vpp_scale = np.linalg.norm(conjugate_generator) * np.linalg.norm(pair_commutator)
            if _roundoff_zero(vpp_trace, vpp_scale):
                continue
            spectator_commutator = _commutator(generators[spectator_index], c_a)
            alp_trace = _trace(vector_generator @ spectator_commutator)
            alp_scale = np.linalg.norm(vector_generator) * np.linalg.norm(spectator_commutator)
            if _roundoff_zero(alp_trace, alp_scale):
                continue
            retained.append(
                VectorExchangeContribution(
                    exchange_name=exchange_name,
                    spectator_index=spectator_index,
                    pair_indices=(i, j),
                    vpp_trace=vpp_trace,
                    alp_trace=alp_trace,
                    flavor_factor=vpp_trace * alp_trace,
                    pole=poles[exchange_name],
                )
            )

    return VectorExchangeAmplitude(
        alp_mass=ma,
        channel=names,
        masses=tuple(meson.mass for meson in mesons),
        f_a=fa,
        f_pi=fpi,
        g_v=gv,
        contributions=tuple(retained),
    )


def build_vector_matching_amplitude(
    alp_mass: float,
    c_q_A: np.ndarray,
    c_q_P: np.ndarray,
    meson_1: str,
    meson_2: str,
    meson_3: str,
    f_a: float,
    *,
    epsilon_k: complex = EPSILON_K,
    g_v: float | None = None,
    pole_parameters: Mapping[
        str, VectorPoleParameters | tuple[float, float] | Mapping[str, float]
    ]
    | None = None,
    f_pi: float = F_PI,
) -> VectorMatchingAmplitude:
    """Build the local term matching antisymmetric tensors to Proca fields.

    The positional and optional inputs intentionally match
    :func:`build_vector_exchange_amplitude`.  The same flavor classification,
    traces, and vector masses are reused.  Widths are accepted and validated
    for interface consistency but do not enter this local contribution.
    ``c_q_P`` is likewise validated but absent at leading ``G_V`` order.
    """
    pole_model = build_vector_exchange_amplitude(
        alp_mass,
        c_q_A,
        c_q_P,
        meson_1,
        meson_2,
        meson_3,
        f_a,
        epsilon_k=epsilon_k,
        g_v=g_v,
        pole_parameters=pole_parameters,
        f_pi=f_pi,
    )
    return VectorMatchingAmplitude(
        alp_mass=pole_model.alp_mass,
        channel=pole_model.channel,
        masses=pole_model.masses,
        f_a=pole_model.f_a,
        f_pi=pole_model.f_pi,
        g_v=pole_model.g_v,
        contributions=pole_model.contributions,
    )
