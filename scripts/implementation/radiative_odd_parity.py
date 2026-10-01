"""On-shell omega -> pi0 gamma predictor in odd-parity RChT.

This transcribes Ruiz-Femenia:2003jdx, Eqs. (27)--(29).  It is not an HLS
model and does not generalize the omega flavor factors to other channels.

Invocation:
    .venv/bin/python -m unittest scripts.implementation.test_radiative_odd_parity
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Mapping

import numpy as np

from .odd_parity_rcht import OddParityCouplings, VJPCouplings, VVPCouplings


OMEGA_PI_GAMMA_PREDICTOR_ID = "rcht.omega_pi0_gamma.width"
OMEGA_PI_GAMMA_BRANCHING_FRACTION_PREDICTOR_ID = "rcht.omega_pi0_gamma.branching_fraction"
Z_OMEGA = 1.0
_PARAMETER_NAMES = tuple(f"c{i}" for i in range(1, 8)) + tuple(f"d{i}" for i in range(1, 5))
_CONFIG_NAMES = ("e", "F", "FV", "MV", "Momega", "mpi")


def _finite_real(value: Any, name: str) -> float:
    if isinstance(value, (complex, np.complexfloating)) and value.imag != 0.0:
        raise ValueError(f"{name} must be finite and real")
    try:
        result = float(np.real(value))
    except (TypeError, ValueError) as error:
        raise ValueError(f"{name} must be finite and real") from error
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite and real")
    return result


def _positive_real(value: Any, name: str) -> float:
    result = _finite_real(value, name)
    if result <= 0.0:
        raise ValueError(f"{name} must be positive")
    return result


@dataclass(frozen=True)
class OmegaPiGammaInputs:
    """Fixed Eq. (27) inputs in GeV; ``MV`` is not ``Momega``.

    ``MV`` normalizes the resonance operators and internal rho exchange;
    ``Momega`` is the external on-shell mass.  The canonical external tensor
    residue is fixed to ``Z_OMEGA=1`` and is intentionally not configurable.
    """

    e: float
    F: float
    FV: float
    MV: float
    Momega: float
    mpi: float

    def __post_init__(self) -> None:
        for name in ("e", "F", "FV", "MV", "Momega"):
            object.__setattr__(self, name, _positive_real(getattr(self, name), name))
        mpi = _finite_real(self.mpi, "mpi")
        if mpi < 0.0:
            raise ValueError("mpi must be nonnegative")
        if mpi >= self.Momega:
            raise ValueError("mpi must be smaller than Momega for an open decay")
        object.__setattr__(self, "mpi", mpi)

    @classmethod
    def from_config(cls, config: Mapping[str, Any]) -> "OmegaPiGammaInputs":
        if not isinstance(config, Mapping):
            raise TypeError("config must be a mapping")
        missing = tuple(name for name in _CONFIG_NAMES if name not in config)
        if missing:
            raise ValueError(f"config is missing required keys: {', '.join(missing)}")
        unexpected = tuple(sorted(set(config) - set(_CONFIG_NAMES)))
        if unexpected:
            raise ValueError(f"config has unexpected keys: {', '.join(unexpected)}")
        return cls(**{name: config[name] for name in _CONFIG_NAMES})


def _couplings_from_mapping(parameters: Mapping[str, Any]) -> OddParityCouplings:
    if not isinstance(parameters, Mapping):
        raise TypeError("parameters must be a mapping")
    missing = tuple(name for name in _PARAMETER_NAMES if name not in parameters)
    if missing:
        raise ValueError(f"parameters are missing required keys: {', '.join(missing)}")
    values = {name: _finite_real(parameters[name], name) for name in _PARAMETER_NAMES}
    return OddParityCouplings(
        VJPCouplings(*(values[f"c{i}"] for i in range(1, 8))),
        VVPCouplings(*(values[f"d{i}"] for i in range(1, 5))),
    )


def omega_pi_gamma_direct_coupling(
    couplings: OddParityCouplings, inputs: OmegaPiGammaInputs
) -> float:
    """Return the direct Eq. (27) scalar coefficient in GeV^-1."""
    if not isinstance(couplings, OddParityCouplings) or not isinstance(inputs, OmegaPiGammaInputs):
        raise TypeError("couplings and inputs must be OddParityCouplings and OmegaPiGammaInputs")
    c = couplings.vjp
    mass_combination = (
        (c.c2 - c.c1 + c.c5 - 2.0 * c.c6) * inputs.Momega**2
        + (c.c1 + c.c2 + 8.0 * c.c3 - c.c5) * inputs.mpi**2
    )
    result = 2.0 * math.sqrt(2.0) * inputs.e * mass_combination / (
        inputs.Momega * inputs.MV * inputs.F
    )
    return _finite_real(result, "direct omega-pi-gamma coupling")


def omega_pi_gamma_rho_coupling(
    couplings: OddParityCouplings, inputs: OmegaPiGammaInputs
) -> float:
    """Return the signed rho-exchange Eq. (27) coefficient in GeV^-1."""
    if not isinstance(couplings, OddParityCouplings) or not isinstance(inputs, OmegaPiGammaInputs):
        raise TypeError("couplings and inputs must be OddParityCouplings and OmegaPiGammaInputs")
    d = couplings.vvp
    mass_combination = d.d3 * inputs.Momega**2 + (
        d.d1 + 8.0 * d.d2 - d.d3
    ) * inputs.mpi**2
    result = -4.0 * inputs.e * inputs.FV * mass_combination / (
        inputs.MV**2 * inputs.Momega * inputs.F
    )
    return _finite_real(result, "rho-exchange omega-pi-gamma coupling")


def G_omega_pi_gamma(couplings: OddParityCouplings, inputs: OmegaPiGammaInputs) -> float:
    """Return the full uneliminated Eq. (27) coefficient in GeV^-1."""
    result = omega_pi_gamma_direct_coupling(couplings, inputs) + omega_pi_gamma_rho_coupling(
        couplings, inputs
    )
    return _finite_real(result, "omega-pi-gamma coupling")


def constrained_G_omega_pi_gamma(inputs: OmegaPiGammaInputs, n_c: int = 3) -> float:
    """Return the OPE-constrained published Eq. (28) coefficient in GeV^-1."""
    if not isinstance(inputs, OmegaPiGammaInputs):
        raise TypeError("inputs must be OmegaPiGammaInputs")
    if isinstance(n_c, bool) or not isinstance(n_c, (int, np.integer)) or n_c <= 0:
        raise ValueError("n_c must be a positive integer")
    result = inputs.e / inputs.FV * (
        int(n_c) * inputs.Momega / (8.0 * math.pi**2 * inputs.F)
        - inputs.F * inputs.Momega / (2.0 * inputs.MV**2)
        * (1.0 + inputs.mpi**2 / inputs.Momega**2)
    )
    return _finite_real(result, "constrained omega-pi-gamma coupling")


def vector_to_pseudoscalar_gamma_width(coupling: complex, Mvector: float, mpseudoscalar: float) -> float:
    """Return ``|G|^2 (MV^2-MP^2)^3/(96*pi*MV^3)`` in GeV."""
    try:
        scalar = complex(coupling)
    except (TypeError, ValueError) as error:
        raise ValueError("coupling must be a finite scalar") from error
    if not math.isfinite(scalar.real) or not math.isfinite(scalar.imag):
        raise ValueError("coupling must be a finite scalar")
    vector_mass = _positive_real(Mvector, "Mvector")
    pseudoscalar_mass = _finite_real(mpseudoscalar, "mpseudoscalar")
    if pseudoscalar_mass < 0.0 or pseudoscalar_mass >= vector_mass:
        raise ValueError("mpseudoscalar must be nonnegative and smaller than Mvector")
    try:
        result = abs(scalar) ** 2 * (vector_mass**2 - pseudoscalar_mass**2) ** 3 / (
            96.0 * math.pi * vector_mass**3
        )
    except OverflowError as error:
        raise ValueError("width must be finite") from error
    return _finite_real(result, "width")


def omega_pi_gamma_width(couplings: OddParityCouplings, inputs: OmegaPiGammaInputs) -> float:
    """Return the uneliminated Eq. (27) on-shell width in GeV."""
    return vector_to_pseudoscalar_gamma_width(
        G_omega_pi_gamma(couplings, inputs), inputs.Momega, inputs.mpi
    )


def omega_pi_gamma_width_predictor(
    parameters: Mapping[str, float], config: Mapping[str, Any]
) -> np.ndarray:
    """Correlated-fit predictor returning one uneliminated width in GeV."""
    inputs = OmegaPiGammaInputs.from_config(config)
    couplings = _couplings_from_mapping(parameters)
    prediction = np.array([omega_pi_gamma_width(couplings, inputs)], dtype=float)
    if prediction.shape != (1,) or not np.all(np.isfinite(prediction)):
        raise ValueError("omega-pi-gamma predictor must return one finite value")
    prediction.setflags(write=False)
    return prediction


def omega_pi_gamma_branching_fraction_predictor(
    parameters: Mapping[str, float], config: Mapping[str, Any]
) -> np.ndarray:
    """Return the uneliminated partial width divided by ``GammaOmegaTotal``."""
    if not isinstance(config, Mapping):
        raise TypeError("config must be a mapping")
    if "GammaOmegaTotal" not in config:
        raise ValueError("config is missing required key: GammaOmegaTotal")
    total_width = _positive_real(config["GammaOmegaTotal"], "GammaOmegaTotal")
    width_config = {name: config[name] for name in config if name != "GammaOmegaTotal"}
    prediction = omega_pi_gamma_width_predictor(parameters, width_config) / total_width
    if prediction.shape != (1,) or not np.all(np.isfinite(prediction)):
        raise ValueError("omega-pi-gamma branching predictor must return one finite value")
    prediction.setflags(write=False)
    return prediction
