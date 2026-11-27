#!/usr/bin/env python3
"""Check whether the eta'-like pole is the heaviest meson-like pole in the only-up model.

Run:   .venv/bin/python scripts/check_up_only_eta_prime_heaviest.py
Tests: .venv/bin/python -m unittest discover -s scripts -p test_check_up_only_eta_prime_heaviest.py

Model and inputs: memory/alp/meson-mixing/finite-mixing/spectrum-and-poles/quadratic-mixing/only-up-finite-mixing-inputs.md
(n_u=1, n_d=n_s=0, m_u = kappa_u f_a/sqrt(2), radial field frozen, LO, finite ALP-meson mixing).
The poles come from the existing solver scripts/implementation/calculate_alp_meson_masses_finite.py
through scripts/implementation/only_up_decays.spectrum. Inputs are those of only_up_inputs.py
(m_pi+-, f_pi, m_u, m_d, m_s, B_0 = m_pi+-^2/(m_u+m_d), m_1^2 fixed from m_eta' = 0.95778 GeV).

Units: masses in GeV; kappa_u and epsilon = f_pi kappa_u/(sqrt(2) m_u) dimensionless.

Identification. The solver labels the six neutral poles (pi0, eta, eta_prime, K0_1, K0_2, a) by the
one-to-one assignment that maximizes the squared overlaps with the unmixed fields. The "eta'-like"
pole is the one labelled eta_prime. The nine meson-like poles are pi0, eta, eta_prime, K0_1, K0_2,
pi+, pi-, K+, K-; the ALP-like pole 'a' is excluded. The radial state is not part of the solver.
At strong mixing the labels are conventional; the check "heaviest" is label independent for the
comparison of the largest of the nine poles with the eta'-like one only if they coincide.

Reference m_sigma: only the lambda-independent term sqrt(B_0 m_u) epsilon of
m_sigma^2 = B_0 m_u epsilon^2 + 2 lambda f_a^2 is given, because lambda is not a sourced input.
It is a lower bound for m_sigma (lambda >= 0 assumed).
"""
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent / "implementation"))
from only_up_inputs import ETA_PRIME_MASS, UP_MASS, B_ZERO, PION_DECAY_CONSTANT, only_up_inputs  # noqa: E402
from only_up_decays import spectrum  # noqa: E402

MESON_LABELS = ("pi0", "eta", "eta_prime", "K0_1", "K0_2", "pi+", "pi-", "K+", "K-")


def kappa_from_epsilon(epsilon):
    """kappa_u = sqrt(2) m_u epsilon / f_pi."""
    return math.sqrt(2) * UP_MASS * epsilon / PION_DECAY_CONSTANT


def epsilon_from_kappa(kappa_u):
    """epsilon = f_pi kappa_u / (sqrt(2) m_u)."""
    return PION_DECAY_CONSTANT * kappa_u / (math.sqrt(2) * UP_MASS)


def default_epsilon_grid(points=25):
    """Logarithmic grid in epsilon from 0.1 to 150."""
    return np.logspace(math.log10(0.1), math.log10(150.), points)


def scan(epsilon_values):
    """Return a list of dicts, one per grid point, with the poles and the checks."""
    rows = []
    for epsilon in epsilon_values:
        kappa_u = kappa_from_epsilon(epsilon)
        masses = spectrum(only_up_inputs(kappa_u))
        meson_masses = {label: masses[label] for label in MESON_LABELS}
        heaviest = max(meson_masses, key=meson_masses.get)
        eta_prime_mass = masses["eta_prime"]
        rows.append({
            "kappa_u": kappa_u, "epsilon": float(epsilon),
            "eta_prime_mass": eta_prime_mass,
            "relative_shift": eta_prime_mass / ETA_PRIME_MASS - 1,
            "heaviest_label": heaviest,
            "eta_prime_is_heaviest": heaviest == "eta_prime",
            "alp_mass": masses["a"],
            "sigma_mass_lambda_independent": math.sqrt(B_ZERO * UP_MASS) * epsilon,
        })
    return rows


def main():
    print("Only-up LO poles (GeV). eta'-like = pole labelled eta_prime by the solver.")
    print("m_sigma(l=0) = sqrt(B_0 m_u) epsilon: lambda-independent part only (lambda not sourced).")
    print(f"{'kappa_u':>10} {'epsilon':>9} {'m_eta_p':>10} {'shift':>10} {'heaviest':>10} {'ok':>3} {'m_a':>10} {'m_sig(l=0)':>11}")
    for row in scan(default_epsilon_grid()):
        print(f"{row['kappa_u']:10.4e} {row['epsilon']:9.4f} {row['eta_prime_mass']:10.6f} "
              f"{row['relative_shift']:10.3e} {row['heaviest_label']:>10} {'Y' if row['eta_prime_is_heaviest'] else 'N':>3} "
              f"{row['alp_mass']:10.5g} {row['sigma_mass_lambda_independent']:11.5g}")


if __name__ == "__main__":
    main()
