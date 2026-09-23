"""Literature phase inputs for elastic scalar K-pi Omnes functions.

The implemented Pelaez--Rodas parametrizations are isospin-symmetric central
fits.  The I=1/2 expression is exposed only through the K-eta threshold, where
the cited analysis switches to an inelastic parametrization.  No uncertainties
or parameter correlations are propagated here.
"""

from dataclasses import dataclass
from typing import Literal

import numpy as np


@dataclass(frozen=True)
class PelaezRodas2016ElasticPhase:
    """Central elastic ``K pi`` S-wave phase in radians.

    Parameters use GeV.  ``fit='CFD'`` selects the dispersively constrained
    fit and is the default; ``fit='UFD'`` selects the unconstrained fit to data.
    """

    isospin: Literal["1/2", "3/2"] = "1/2"
    fit: Literal["CFD", "UFD"] = "CFD"

    # Pelaez:2016tgi uses these isospin-limit masses in its phase fits.
    pion_mass: float = 0.13957
    kaon_mass: float = 0.496
    eta_mass: float = 0.547

    def __post_init__(self) -> None:
        if self.isospin not in ("1/2", "3/2"):
            raise ValueError("isospin must be '1/2' or '3/2'")
        if self.fit not in ("CFD", "UFD"):
            raise ValueError("fit must be 'CFD' or 'UFD'")

    @property
    def threshold(self) -> float:
        return (self.kaon_mass + self.pion_mass) ** 2

    @property
    def maximum(self) -> float:
        if self.isospin == "1/2":
            # Pelaez:2016tgi switches the I=1/2 wave to an inelastic form here.
            return (self.kaon_mass + self.eta_mass) ** 2
        # Pelaez:2016tgi reports no evidence for I=3/2 inelasticity to ~1.8 GeV.
        return 1.8**2

    @staticmethod
    def _momentum(s: float, first_mass: float, second_mass: float) -> float:
        product = (s - (first_mass + second_mass) ** 2) * (
            s - (first_mass - second_mass) ** 2
        )
        return np.sqrt(max(0.0, product)) / (2.0 * np.sqrt(s))

    def __call__(self, invariant: float) -> float:
        try:
            s = float(invariant)
        except (TypeError, ValueError) as error:
            raise TypeError("invariant must be a real scalar") from error
        if not np.isfinite(s):
            raise ValueError("invariant must be finite")
        tolerance = 64.0 * np.finfo(float).eps * max(1.0, self.maximum)
        if self.threshold - tolerance <= s < self.threshold:
            s = self.threshold
        if self.maximum < s <= self.maximum + tolerance:
            s = self.maximum
        if s < self.threshold or s > self.maximum:
            raise ValueError(
                f"phase is valid only for {self.threshold} <= s <= {self.maximum} GeV^2"
            )
        if s == self.threshold:
            return 0.0

        m_k = self.kaon_mass
        m_pi = self.pion_mass
        q = self._momentum(s, m_k, m_pi)
        delta_kpi = m_k**2 - m_pi**2
        y = ((s - delta_kpi) / (s + delta_kpi)) ** 2

        if self.isospin == "1/2":
            # Pelaez:2016tgi gives the conformal fit and parameters;
            # Pelaez:2020gnd gives the Adler-zero expression consistent with
            # the quoted s_A ~= 0.236 GeV^2.
            alpha = 1.15
            conformal_endpoint = 1.1**2
            b0, b1 = (0.411, 0.162) if self.fit == "CFD" else (0.411, 0.181)
            sigma = m_k**2 + m_pi**2
            adler = (
                sigma + 2.0 * np.sqrt(4.0 * delta_kpi**2 + m_k**2 * m_pi**2)
            ) / 5.0
            coefficients = (b0, b1)
            adler_factor = s - adler
        else:
            # Pelaez:2016tgi, Eqs. (conformal variable, cot delta_0^{3/2})
            # and Table of S^{3/2} parameters.
            alpha = 1.4
            conformal_endpoint = 1.84**2
            coefficients = (
                (2.27, 3.94, 3.36)
                if self.fit == "CFD"
                else (2.25, 4.21, 2.45)
            )
            adler = m_k**2 + m_pi**2
            adler_factor = adler - s

        y0 = (
            (conformal_endpoint - delta_kpi)
            / (conformal_endpoint + delta_kpi)
        ) ** 2
        conformal_root = np.sqrt(max(0.0, y0 - y))
        omega = (np.sqrt(y) - alpha * conformal_root) / (
            np.sqrt(y) + alpha * conformal_root
        )
        polynomial = sum(value * omega**power for power, value in enumerate(coefficients))
        cotangent = np.sqrt(s) * polynomial / (2.0 * q * adler_factor)
        return float(np.arctan(1.0 / cotangent))
