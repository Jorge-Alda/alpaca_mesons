"""Numerical once-subtracted elastic Omnes functions.

Invariant masses and subtraction points are expressed in GeV^2, while phase
callbacks return radians.  A phase callback and its high-energy continuation
are explicit inputs: this module does not embed a particular K-pi phase fit.

Example:
    omega = OmnesFunction(
        phase=lambda s: phase_fit(s),
        threshold=(0.493677 + 0.1349768) ** 2,
        phase_input_max=1.3**2,
        subtraction_point=0.0,
        tail="constant",
        asymptotic_phase=np.pi,
    )
    value = omega(0.8)
"""

from dataclasses import dataclass
from typing import Callable, Literal, Mapping

import numpy as np
from scipy.integrate import quad


PhaseFunction = Callable[[float], float]


def _real_finite(value: float, name: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise TypeError(f"{name} must be a real scalar") from error
    if not np.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


@dataclass(frozen=True)
class OmnesFunction:
    """A normalized, once-subtracted scalar Omnes function.

    ``phase`` is used from ``threshold`` through ``phase_input_max``.  With
    ``tail='constant'``, the phase above that endpoint is
    ``asymptotic_phase``; continuity at the matching point is required.  With
    ``tail='truncate'``, the phase is set to zero above the endpoint, which is
    an explicit cutoff model rather than a physical asymptotic prediction.

    Only real values of ``s`` are accepted.  On the upper rim of the elastic
    cut, the returned phase is the supplied scattering phase.  Exact
    evaluation at the integration endpoints is excluded because it requires a
    separately specified one-sided limiting prescription.
    """

    phase: PhaseFunction
    threshold: float
    phase_input_max: float
    subtraction_point: float = 0.0
    tail: Literal["constant", "truncate"] = "constant"
    asymptotic_phase: float | None = None
    epsabs: float = 1e-10
    epsrel: float = 1e-9
    limit: int = 300

    def __post_init__(self) -> None:
        threshold = _real_finite(self.threshold, "threshold")
        endpoint = _real_finite(self.phase_input_max, "phase_input_max")
        subtraction = _real_finite(self.subtraction_point, "subtraction_point")
        if not subtraction < threshold < endpoint:
            raise ValueError(
                "require subtraction_point < threshold < phase_input_max"
            )
        if self.tail not in ("constant", "truncate"):
            raise ValueError("tail must be 'constant' or 'truncate'")
        if self.epsabs <= 0.0 or self.epsrel <= 0.0 or self.limit < 1:
            raise ValueError("quadrature tolerances and limit must be positive")

        endpoint_phase = self._phase(endpoint)
        if self.tail == "constant":
            if self.asymptotic_phase is None:
                raise ValueError(
                    "asymptotic_phase is required for a constant phase tail"
                )
            asymptotic = _real_finite(self.asymptotic_phase, "asymptotic_phase")
            if not np.isclose(endpoint_phase, asymptotic, rtol=1e-8, atol=1e-10):
                raise ValueError(
                    "constant tail must match phase(phase_input_max) continuously"
                )
        elif self.asymptotic_phase is not None:
            raise ValueError("asymptotic_phase is only valid with tail='constant'")

    def _phase(self, invariant: float) -> float:
        return _real_finite(self.phase(float(invariant)), "phase(s)")

    def __call__(
        self, invariant: float, pair: tuple[str, str] | None = None
    ) -> complex:
        """Evaluate ``Omega(s; s0)`` on the physical upper rim of the cut.

        ``pair`` is accepted and ignored so a single object can be passed
        directly to the ``fsi`` argument of ``s_wave_decay_width``.
        """
        del pair
        s = _real_finite(invariant, "invariant")
        s0 = float(self.subtraction_point)
        threshold = float(self.threshold)
        endpoint = float(self.phase_input_max)
        if s == s0:
            return 1.0 + 0.0j
        if s == threshold or s == endpoint:
            raise ValueError(
                "evaluate an explicit one-sided limit at an integration endpoint"
            )

        if self.tail == "constant":
            delta_infinity = float(self.asymptotic_phase)
            # Omnes:1958hv: the once-subtracted elastic solution is the
            # exponential of the phase dispersion integral.  The constant
            # phase from threshold to infinity is integrated analytically.
            real_exponent = (delta_infinity / np.pi) * np.log(
                (threshold - s0) / abs(threshold - s)
            )
            imaginary_exponent = delta_infinity if s > threshold else 0.0

            def numerator(x: float) -> float:
                return (self._phase(x) - delta_infinity) / (x - s0)

        else:
            real_exponent = 0.0
            imaginary_exponent = 0.0

            def numerator(x: float) -> float:
                return self._phase(x) / (x - s0)

        if threshold < s < endpoint:
            principal_value, _ = quad(
                numerator,
                threshold,
                endpoint,
                weight="cauchy",
                wvar=s,
                epsabs=self.epsabs,
                epsrel=self.epsrel,
                limit=self.limit,
            )
            # Watson:1952ji and Omnes:1958hv: on the elastic upper rim, the
            # Omnes phase equals the elastic scattering phase.
            if self.tail == "constant":
                imaginary_exponent += self._phase(s) - delta_infinity
            else:
                imaginary_exponent = self._phase(s)
        else:
            principal_value, _ = quad(
                lambda x: numerator(x) / (x - s),
                threshold,
                endpoint,
                epsabs=self.epsabs,
                epsrel=self.epsrel,
                limit=self.limit,
            )

        real_exponent += (s - s0) * principal_value / np.pi
        return complex(np.exp(real_exponent + 1j * imaginary_exponent))


class PairOmnesFSI:
    """Dispatch pair-dependent scalar Omnes functions to decay bands.

    This remains a diagonal independent-band approximation.  It is not the
    charge-basis Muskhelishvili--Omnes matrix required for exact coupled-channel
    unitarity.  Jamin:2001zq discusses the matrix treatment when inelastic
    channels are retained.
    """

    def __init__(self, functions: Mapping[tuple[str, str], OmnesFunction]):
        if not functions:
            raise ValueError("at least one pair-dependent Omnes function is required")
        normalized: dict[frozenset[str], OmnesFunction] = {}
        for pair, function in functions.items():
            if len(pair) != 2 or pair[0] == pair[1]:
                raise ValueError("each pair key must contain two distinct mesons")
            if not isinstance(function, OmnesFunction):
                raise TypeError("pair values must be OmnesFunction instances")
            key = frozenset(pair)
            if key in normalized:
                raise ValueError(f"duplicate unordered pair {pair!r}")
            normalized[key] = function
        self._functions = normalized

    def __call__(self, invariant: float, pair: tuple[str, str]) -> complex:
        try:
            function = self._functions[frozenset(pair)]
        except KeyError as error:
            raise KeyError(f"no Omnes function configured for pair {pair!r}") from error
        return function(invariant)
