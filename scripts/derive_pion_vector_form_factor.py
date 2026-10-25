#!/usr/bin/env python3
"""Calculate a central three-subtracted pion vector form factor.

Run with `.venv/bin/python scripts/derive_pion_vector_form_factor.py`.
All masses are GeV; invariant masses squared are GeV^2. The prescription
is GomezDumm:2013sib, exact-isospin fit I, with a declared phase continuation
to pi at M_tau^2 and integration to infinity. This is a central model,
without electromagnetic corrections or statistical replicas.
"""
from __future__ import annotations

import argparse
import cmath
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad

ROOT = Path(__file__).resolve().parents[1]


class PionVectorFormFactor:
    """Three-subtracted Omnes form on the upper rim of the physical cut."""

    def __init__(self, phase_join=0.98, phase_end=1.77686**2, tolerance=2e-10):
        self.rho_mass = 0.8430
        self.pion_decay_constant = 0.0901
        self.alpha1 = 1.87
        self.alpha2 = 4.29
        # Physical charged masses used as the common isospin-limit masses.
        self.pion_mass = 0.13957039
        self.kaon_mass = 0.493677
        self.threshold = 4 * self.pion_mass**2
        self.phase_join = phase_join
        self.phase_end = phase_end
        self.tolerance = tolerance

    def loop(self, invariant_mass_squared, meson_mass):
        """Return A_P(s,M_rho^2), with Im A_P=-pi sigma_P^3 above threshold."""
        sigma = cmath.sqrt(complex(1 - 4 * meson_mass**2 / invariant_mass_squared))
        if invariant_mass_squared > 4 * meson_mass**2:
            logarithm = np.log((1 + sigma.real) / (1 - sigma.real)) - 1j * np.pi
        elif sigma == 0:
            logarithm = 0j
        else:
            logarithm = cmath.log((sigma + 1) / (sigma - 1))
        return (np.log(meson_mass**2 / self.rho_mass**2)
                + 8 * meson_mass**2 / invariant_mass_squared - 5 / 3
                + sigma**3 * logarithm)

    def elastic_phase(self, invariant_mass_squared):
        """Continuous phase of the resummed one-loop form, in radians."""
        denominator = (self.rho_mass**2 * (1 + invariant_mass_squared
                       / (96 * np.pi**2 * self.pion_decay_constant**2)
                       * (self.loop(invariant_mass_squared, self.pion_mass)
                          + self.loop(invariant_mass_squared, self.kaon_mass) / 2))
                       - invariant_mass_squared)
        return float(np.angle(self.rho_mass**2 / denominator))

    def phase(self, invariant_mass_squared):
        """Elastic phase followed by linear continuation to pi."""
        if invariant_mass_squared <= self.threshold:
            return 0.0
        if invariant_mass_squared <= self.phase_join:
            return self.elastic_phase(invariant_mass_squared)
        if invariant_mass_squared < self.phase_end:
            join_phase = self.elastic_phase(self.phase_join)
            return join_phase + (np.pi - join_phase) * (
                invariant_mass_squared - self.phase_join) / (self.phase_end - self.phase_join)
        return np.pi

    def evaluate(self, invariant_mass_squared):
        """Return F(s+i0) using an explicit Cauchy principal value and cut phase."""
        if invariant_mass_squared == 0:
            return 1 + 0j
        def spectral_weight(value):
            return self.phase(value) / value**3
        options = dict(epsabs=self.tolerance, epsrel=self.tolerance, limit=300)
        if self.threshold < invariant_mass_squared < self.phase_end:
            integral = quad(spectral_weight, self.threshold, self.phase_end,
                            weight='cauchy', wvar=invariant_mass_squared, **options)[0]
        else:
            integral = quad(lambda value: spectral_weight(value)
                            / (value - invariant_mass_squared), self.threshold,
                            self.phase_end, points=[self.phase_join], **options)[0]
        integral += quad(lambda value: np.pi / (value**3 * (value - invariant_mass_squared)),
                         self.phase_end, np.inf, **options)[0]
        real_exponent = (self.alpha1 * invariant_mass_squared
                         + self.alpha2 * invariant_mass_squared**2 / 2
                         + invariant_mass_squared**3 * integral / np.pi)
        return cmath.exp(real_exponent + 1j * self.phase(invariant_mass_squared))


def main():
    """Write the finite complex central grid and deterministic numerical diagnostics."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--phase-join', type=float, default=0.98)
    parser.add_argument('--phase-end', type=float, default=1.77686**2)
    parser.add_argument('--maximum-energy', type=float, default=1.3)
    parser.add_argument('--points', type=int, default=401)
    parser.add_argument('--output', type=Path,
                        default=ROOT / 'data/odd_parity/pion_vector_form_factor_central.json')
    arguments = parser.parse_args()
    form_factor = PionVectorFormFactor(arguments.phase_join, arguments.phase_end)
    invariant_masses_squared = np.linspace(form_factor.threshold,
                                         arguments.maximum_energy**2, arguments.points)
    values = np.array([form_factor.evaluate(value) for value in invariant_masses_squared])
    selected_points = [0.0, form_factor.threshold, 0.1, 0.3, 0.6, 0.98, 1.0, 1.5]
    refined = PionVectorFormFactor(arguments.phase_join, arguments.phase_end, 2e-12)
    differences = [abs(form_factor.evaluate(value) - refined.evaluate(value))
                   for value in selected_points]
    phase_residuals = np.angle(values * np.exp(-1j * np.array([
        form_factor.phase(value) for value in invariant_masses_squared])))
    result = {
        'bibtex_keys': ['GomezDumm:2013sib'],
        'convention': 'F_pi^V(0)=1; s+i0 upper rim; dimensionless F_pi^V',
        'units': {'x': 'GeV^2', 'masses': 'GeV', 'phase': 'radian'},
        'prescription': 'three-subtracted Omnes, exact-isospin fit I central parameters',
        'parameters': {'M_rho_GeV': form_factor.rho_mass,
                       'F_pi_GeV': form_factor.pion_decay_constant,
                       'alpha1_GeV^-2': form_factor.alpha1,
                       'alpha2_GeV^-4': form_factor.alpha2,
                       'm_pi_GeV': form_factor.pion_mass,
                       'm_K_GeV': form_factor.kaon_mass,
                       's1_GeV^2': form_factor.phase_join,
                       's2_GeV^2': form_factor.phase_end},
        'assumptions': [
            'Charged-pion and charged-kaon masses used as common isospin-limit masses.',
            'Linear phase from s1 to pi at s2; constant pi to infinity.',
            'No electromagnetic corrections, rho-omega mixing, or excited-rho matching.',
            'Central model only; fit covariance not available; no statistical replicas.',
            'Grid above 1.5 GeV^2 extrapolates the low-energy representation.'],
        'variation_options': {'s1_GeV^2': [0.95, 1.1],
                              'description': 'Use CLI phase-join/phase-end for discrete phase variations.'},
        'x': invariant_masses_squared.tolist(),
        'FpiV_real': values.real.tolist(), 'FpiV_imag': values.imag.tolist(),
        'diagnostics': {'max_tolerance_change_absolute': max(differences),
                        'max_phase_residual_radian': float(np.max(abs(phase_residuals))),
                        'normalization': [1.0, 0.0],
                        'selected_values': [{'x_GeV^2': value,
                                             'real': form_factor.evaluate(value).real,
                                             'imag': form_factor.evaluate(value).imag}
                                            for value in selected_points]}}
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(result, indent=2) + '\n')
    print(arguments.output.relative_to(ROOT))


if __name__ == '__main__':
    main()
