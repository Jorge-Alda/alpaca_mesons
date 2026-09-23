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
from scipy.interpolate import PchipInterpolator


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
    evaluation at an integration endpoint is obtained from its continuous
    two-sided numerical limit.
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
        if s == threshold:
            step = 1e-10 * max(1.0, abs(s))
            left = self(s - step)
            right = self(s + step)
            return complex(0.5 * (left.real + right.real), 0.0)
        if s == endpoint:
            step = 1e-8 * max(1.0, abs(s))
            return 0.5 * (self(s - step) + self(s + step))

        # Omnes:1958hv: the once-subtracted elastic solution is the
        # exponential of the phase dispersion integral.  Use
        # x=threshold+u^2 so the square-root threshold behavior is smooth.
        upper_u = np.sqrt(endpoint - threshold)

        def transformed_numerator(u: float) -> float:
            x = threshold + u * u
            return 2.0 * u * self._phase(x) / (x - s0)

        real_exponent = 0.0
        imaginary_exponent = 0.0
        if threshold < s < endpoint:
            singular_u = np.sqrt(s - threshold)

            def cauchy_numerator(u: float) -> float:
                return transformed_numerator(u) / (u + singular_u)

            numerator_at_s = cauchy_numerator(singular_u)
            derivative_step = min(singular_u, upper_u - singular_u) * 1e-4
            derivative = (
                cauchy_numerator(singular_u + derivative_step)
                - cauchy_numerator(singular_u - derivative_step)
            ) / (2.0 * derivative_step)

            def regular_integrand(u: float) -> float:
                if u == singular_u:
                    return derivative
                return (cauchy_numerator(u) - numerator_at_s) / (
                    u - singular_u
                )

            principal_value, _ = quad(
                regular_integrand,
                0.0,
                upper_u,
                epsabs=self.epsabs,
                epsrel=self.epsrel,
                limit=self.limit,
            )
            principal_value += numerator_at_s * np.log(
                (upper_u - singular_u) / singular_u
            )
            # Watson:1952ji and Omnes:1958hv: on the elastic upper rim, the
            # Omnes phase equals the elastic scattering phase.
            imaginary_exponent = self._phase(s)
        else:
            def regular_integrand(u: float) -> float:
                denominator = threshold + u * u - s
                if u == 0.0 and denominator == 0.0:
                    probe = 1e-8 * upper_u
                    return transformed_numerator(probe) / (probe * probe)
                return transformed_numerator(u) / denominator

            principal_value, _ = quad(
                regular_integrand,
                0.0,
                upper_u,
                epsabs=self.epsabs,
                epsrel=self.epsrel,
                limit=self.limit,
            )

        real_exponent += (s - s0) * principal_value / np.pi
        if self.tail == "constant":
            delta_infinity = float(self.asymptotic_phase)
            # The constant phase above the finite input endpoint is integrated
            # analytically, including its upper-rim phase when s is above it.
            real_exponent += (delta_infinity / np.pi) * np.log(
                (endpoint - s0) / abs(endpoint - s)
            )
            if s > endpoint:
                imaginary_exponent = delta_infinity
        return complex(np.exp(real_exponent + 1j * imaginary_exponent))


class TabulatedOmnes:
    """Adaptive PCHIP interpolation of an expensive ``OmnesFunction``.

    The logarithm of the magnitude is interpolated in a signed square-root
    threshold coordinate, while the inexpensive input phase is retained
    exactly.  Each refinement step checks two omitted interior points per
    interval against direct dispersion-integral evaluations.  Construction is
    expensive once; subsequent calls require no dispersion quadrature.
    """

    def __init__(
        self,
        function: OmnesFunction,
        minimum: float,
        maximum: float,
        relative_tolerance: float = 1e-6,
        initial_points: int = 17,
        maximum_points: int = 4097,
    ) -> None:
        if not isinstance(function, OmnesFunction):
            raise TypeError("function must be an OmnesFunction")
        self.function = function
        self.minimum = _real_finite(minimum, "minimum")
        self.maximum = _real_finite(maximum, "maximum")
        self.relative_tolerance = _real_finite(
            relative_tolerance, "relative_tolerance"
        )
        if not self.minimum < self.maximum:
            raise ValueError("minimum must be smaller than maximum")
        if self.relative_tolerance <= 0.0:
            raise ValueError("relative_tolerance must be positive")
        if initial_points < 5 or maximum_points < initial_points:
            raise ValueError("require 5 <= initial_points <= maximum_points")
        self.maximum_points = int(maximum_points)

        coordinate_minimum = self._coordinate(self.minimum)
        coordinate_maximum = self._coordinate(self.maximum)
        coordinate_grid = np.linspace(
            coordinate_minimum, coordinate_maximum, int(initial_points)
        )
        grid = self._inverse_coordinate(coordinate_grid)
        # The signed square-root coordinate resolves the threshold cusp.
        # Include exactly known normalization and phase-matching landmarks.
        landmarks = (
            function.subtraction_point,
            function.threshold,
            function.phase_input_max,
        )
        grid = np.unique(
            np.concatenate(
                [grid, [point for point in landmarks if self.minimum < point < self.maximum]]
            )
        )
        values = np.asarray([function(float(x)) for x in grid], dtype=complex)
        self.max_validation_error = np.inf

        while True:
            coordinates = self._coordinate(grid)
            magnitude_spline = self._build_spline(coordinates, values)
            left = coordinates[:-1]
            widths = np.diff(coordinates)
            candidate_coordinates = np.sort(
                np.concatenate((left + widths / 3.0, left + 2.0 * widths / 3.0))
            )
            candidates = self._inverse_coordinate(candidate_coordinates)
            direct = np.asarray(
                [function(float(x)) for x in candidates], dtype=complex
            )
            approximate = self._evaluate_spline(
                candidates, candidate_coordinates, magnitude_spline
            )
            errors = np.abs(approximate - direct) / np.abs(direct)
            self.max_validation_error = float(np.max(errors, initial=0.0))
            failing = errors > self.relative_tolerance
            if not np.any(failing):
                break
            additions = candidates[failing]
            if grid.size + additions.size > self.maximum_points:
                raise RuntimeError(
                    "Omnes interpolation did not reach the requested tolerance "
                    f"before maximum_points={self.maximum_points}; "
                    f"last maximum relative error={self.max_validation_error:.3e}"
                )
            grid = np.concatenate((grid, additions))
            values = np.concatenate((values, direct[failing]))
            previous_size = grid.size - additions.size
            grid, unique_indices = np.unique(grid, return_index=True)
            values = values[unique_indices]
            if grid.size == previous_size:
                raise RuntimeError(
                    "Omnes interpolation refinement stagnated at floating-point "
                    f"resolution; last maximum relative error="
                    f"{self.max_validation_error:.3e}"
                )

        self.grid = grid
        self.values = values
        self._log_magnitude = self._build_spline(
            self._coordinate(grid), values
        )

    def _coordinate(self, invariant: np.ndarray | float) -> np.ndarray:
        displacement = np.asarray(invariant) - self.function.threshold
        return np.sign(displacement) * np.sqrt(np.abs(displacement))

    def _inverse_coordinate(self, coordinate: np.ndarray | float) -> np.ndarray:
        value = np.asarray(coordinate)
        return self.function.threshold + np.sign(value) * value**2

    @staticmethod
    def _build_spline(
        grid: np.ndarray, values: np.ndarray
    ) -> tuple[PchipInterpolator | None, PchipInterpolator | None]:
        log_magnitude = np.log(np.abs(values))
        lower = grid <= 0.0
        upper = grid >= 0.0
        lower_spline = (
            PchipInterpolator(grid[lower], log_magnitude[lower], extrapolate=False)
            if np.count_nonzero(lower) >= 2
            else None
        )
        upper_spline = (
            PchipInterpolator(grid[upper], log_magnitude[upper], extrapolate=False)
            if np.count_nonzero(upper) >= 2
            else None
        )
        return lower_spline, upper_spline

    @staticmethod
    def _log_magnitude_value(
        coordinate: np.ndarray | float,
        splines: tuple[PchipInterpolator | None, PchipInterpolator | None],
    ) -> np.ndarray:
        values = np.atleast_1d(np.asarray(coordinate, dtype=float))
        result = np.empty_like(values)
        lower = values <= 0.0
        lower_spline, upper_spline = splines
        if np.any(lower):
            spline = lower_spline if lower_spline is not None else upper_spline
            result[lower] = spline(values[lower])
        if np.any(~lower):
            spline = upper_spline if upper_spline is not None else lower_spline
            result[~lower] = spline(values[~lower])
        return result.reshape(np.shape(coordinate))

    def _physical_phase(self, invariant: np.ndarray | float) -> np.ndarray:
        values = np.atleast_1d(np.asarray(invariant, dtype=float))
        phases = np.zeros_like(values)
        elastic = (values > self.function.threshold) & (
            values <= self.function.phase_input_max
        )
        phases[elastic] = [self.function._phase(float(s)) for s in values[elastic]]
        if self.function.tail == "constant":
            phases[values > self.function.phase_input_max] = float(
                self.function.asymptotic_phase
            )
        return phases.reshape(np.shape(invariant))

    def _evaluate_spline(
        self,
        invariant: np.ndarray | float,
        coordinate: np.ndarray | float,
        magnitude_spline: tuple[PchipInterpolator | None, PchipInterpolator | None],
    ) -> np.ndarray:
        return np.exp(
            self._log_magnitude_value(coordinate, magnitude_spline)
            + 1j * self._physical_phase(invariant)
        )

    @property
    def grid_size(self) -> int:
        return int(self.grid.size)

    def __call__(
        self, invariant: float, pair: tuple[str, str] | None = None
    ) -> complex:
        del pair
        s = _real_finite(invariant, "invariant")
        if s < self.minimum or s > self.maximum:
            raise ValueError(
                f"invariant must lie in [{self.minimum}, {self.maximum}] GeV^2"
            )
        if s == self.function.subtraction_point:
            return 1.0 + 0.0j
        return complex(
            self._evaluate_spline(
                s, self._coordinate(s), self._log_magnitude
            )
        )


class PairOmnesFSI:
    """Dispatch pair-dependent scalar Omnes functions to decay bands.

    This remains a diagonal independent-band approximation.  It is not the
    charge-basis Muskhelishvili--Omnes matrix required for exact coupled-channel
    unitarity.  Jamin:2001zq discusses the matrix treatment when inelastic
    channels are retained.
    """

    def __init__(
        self,
        functions: Mapping[tuple[str, str], OmnesFunction | TabulatedOmnes],
    ):
        if not functions:
            raise ValueError("at least one pair-dependent Omnes function is required")
        normalized: dict[frozenset[str], OmnesFunction] = {}
        for pair, function in functions.items():
            if len(pair) != 2 or pair[0] == pair[1]:
                raise ValueError("each pair key must contain two distinct mesons")
            if not isinstance(function, (OmnesFunction, TabulatedOmnes)):
                raise TypeError(
                    "pair values must be OmnesFunction or TabulatedOmnes instances"
                )
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
