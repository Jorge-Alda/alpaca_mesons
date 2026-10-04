"""Diagnose whether the heavier isoscalar states can be excluded from the BaBar 3pi data.

Scoping calculation for the odd-parity RChT coupling fit.  It answers one
question: if we keep the single-multiplet baseline fixed at the level of
rho, omega and phi, is there any mass window in the tabulated
e+e- -> pi+ pi- pi0 spectrum that is not already dominated by the
omega'(1350) and omega''(1660) states that BaBar had to introduce?

Method.  Inside the coherent sum of BaBar Eq. (14) the magnitude of the
contribution of a resonance V is

    w_V(m) = Gamma_V m_V^(3/2) sqrt(B(V->e+e-) B(V->3pi)) / |D_V(m)|,
    D_V(m) = m_V^2 - m^2 - i m Gamma_V,

which is independent of the common factor F_rhopi(m) and of the phases.  The
relative sizes of w_V therefore answer the separability question without
needing F_rhopi or the interference phases, neither of which is tabulated.
This is NOT a prediction of the cross section and is not a fit.

Conventions and units.  Natural units; masses and widths in MeV inside the
propagator; resonance parameters and fitted branching-fraction products as
printed in BaBar:2004ytv, Summary Eqs. (19)-(20).  Natural units for the
covariance-free scoping arithmetic.  Uncertainties are the quoted statistical
and systematic split of the table, used only to compare systematic against
statistical size.

Inputs.
  data/odd_parity/babar_3pi_born_spectrum.csv, checksum-validated by
  scripts/extract_odd_parity_tables.py; see
  memory/extracts/BaBar:2004ytv/mass-binning-and-fitted-products.md.
  Resonance parameters from BaBar:2004ytv [@BaBar:2004ytv].

Invocation:
    .venv/bin/python scripts/diagnose_three_pion_resonance_separation.py
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np


DATA_FILE = Path(__file__).resolve().parents[1] / "data" / "odd_parity" / "babar_3pi_born_spectrum.csv"

# name -> (m_V [MeV], Gamma_V [MeV], dm_V [MeV], dGamma_V [MeV], B*B, d(B*B)/B*B)
# BaBar:2004ytv Summary Eqs. (19)-(20); omega and phi masses/widths from PDG.
RESONANCES: dict[str, tuple[float, float, float, float, float, float]] = {
    "omega": (782.7, 8.68, 0.0, 0.0, 6.70e-5, 0.0),
    "phi": (1019.5, 4.25, 0.0, 0.0, 4.30e-5, 0.0),
    "omega'": (1350.0, 450.0, 20.0, 70.0, 0.82e-6, 0.12),
    "omega''": (1660.0, 230.0, 10.0, 30.0, 1.30e-6, 0.08),
}
HEAVY = ("omega'", "omega''")
PROBE_MASSES = (1.0625, 1.1625, 1.3125, 1.6375, 1.9875, 2.9875)
CLEAN_THRESHOLD = 0.10


def load_table(path: Path = DATA_FILE) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return sorted (sqrt_s [GeV], sigma [nb], sigma_stat [nb], sigma_syst [nb])."""
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"no rows in {path.name}")
    columns = [np.array([float(row[name]) for row in rows]) for name in
               ("sqrt_s_GeV", "sigma_nb", "sigma_stat_nb", "sigma_syst_nb")]
    order = np.argsort(columns[0])
    return tuple(column[order] for column in columns)  # type: ignore[return-value]


def weight(m_V: float, gamma: float, branching_product: float, m: np.ndarray) -> np.ndarray:
    """Magnitude of one resonance term of the BaBar Eq. (14) coherent sum."""
    propagator = m_V**2 - m**2 - 1j * m * gamma
    return gamma * m_V**1.5 * np.sqrt(branching_product) / np.abs(propagator)


def weights(masses: np.ndarray,
            parameters: dict[str, tuple[float, float, float]] | None = None) -> dict[str, np.ndarray]:
    """Weight of every resonance at every mass.

    ``parameters`` maps a name to ``(m_V, Gamma_V, branching product)`` and
    defaults to the tabulated values, discarding the uncertainty columns.
    """
    if parameters is None:
        parameters = {name: (value[0], value[1], value[4]) for name, value in RESONANCES.items()}
    return {name: weight(m_V, gamma, product, masses)
            for name, (m_V, gamma, product) in parameters.items()}


def lorentzian_floor(name: str) -> float:
    """``min_{m>0} |D(m)|^-2`` relative to the peak value; equals ``(Gamma/m_V)^2``."""
    m_V, gamma = RESONANCES[name][0], RESONANCES[name][1]
    return (gamma / m_V) ** 2


def separation_in_widths(name: str) -> float:
    """Distance of a heavier state above the phi mass, in units of its own width."""
    m_phi, gamma_phi = RESONANCES["phi"][0], RESONANCES["phi"][1]
    del gamma_phi
    return (RESONANCES[name][0] - m_phi) / RESONANCES[name][1]


def worst_case_shift(masses: np.ndarray, base_fraction: np.ndarray, index: int) -> float:
    """Largest shift of the combined heavy-state weight fraction from one parameter."""
    baseline = sum(weights(masses).values())
    shifts = []
    for name in HEAVY:
        m_V, gamma, d_m, d_gamma, product, d_product = RESONANCES[name]
        for perturbed in ((m_V + d_m, gamma, product), (m_V, gamma + d_gamma, product),
                          (m_V, gamma, product * (1.0 + d_product))):
            parameters = {key: (value[0], value[1], value[4]) for key, value in RESONANCES.items()}
            parameters[name] = perturbed
            trial = weights(masses, parameters)
            heavy_fraction = sum(trial[key] for key in HEAVY) / sum(trial.values())
            shifts.append(abs(heavy_fraction[index] - base_fraction[index]))
    return max(shifts)


def main() -> None:
    sqrt_s, sigma, sigma_stat, sigma_syst = load_table()
    masses = sqrt_s * 1e3

    base = weights(masses)
    total = sum(base.values())
    fraction = {name: base[name] / total for name in RESONANCES}
    heavy_fraction = sum(fraction[name] for name in HEAVY)

    print("=" * 72)
    print("BaBar 3pi: can the one-multiplet baseline be isolated in any mass window?")
    print("=" * 72)
    print(f"bins: {sqrt_s.size}   sqrt_s range: {sqrt_s[0]:.4f} - {sqrt_s[-1]:.4f} GeV"
          f"   bin width: {np.diff(sqrt_s)[0] * 1e3:.1f} MeV")
    print(f"omega at {RESONANCES['omega'][0]:.1f} MeV inside table: "
          f"{sqrt_s[0] * 1e3 <= RESONANCES['omega'][0]}")
    print(f"phi   at {RESONANCES['phi'][0]:.1f} MeV inside table: "
          f"{sqrt_s[0] * 1e3 <= RESONANCES['phi'][0]}")

    print("\n-- 1. Separation of the heavier states from the phi, in their own widths")
    for name in HEAVY:
        m_V = RESONANCES[name][0]
        print(f"  {name:8s} m-m_phi = {m_V - RESONANCES['phi'][0]:6.1f} MeV"
              f" = {separation_in_widths(name):5.2f} Gamma_V")

    print("\n-- 2. Lorentzian floor: min_{m>0} |D|^-2 relative to the peak value")
    for name in RESONANCES:
        m_V, gamma = RESONANCES[name][0], RESONANCES[name][1]
        print(f"  {name:8s} Gamma/m_V = {gamma / m_V:6.4f}   floor = {lorentzian_floor(name):9.5f}"
              f" ({100 * lorentzian_floor(name):6.2f}% of peak)")

    print("\n-- 3. Weight fraction inside the BaBar coherent sum")
    header = "  m[GeV]  sigma[nb] |" + "".join(f"{name:>9s}" for name in RESONANCES)
    print(header)
    print("-" * len(header))
    for target in PROBE_MASSES:
        index = int(np.argmin(np.abs(sqrt_s - target)))
        cells = "".join(f"{100 * fraction[name][index]:8.1f}%" for name in RESONANCES)
        print(f" {sqrt_s[index]:7.4f} {sigma[index]:9.2f} |{cells}")

    print("\n-- 4. Systematic from fixing the heavier states at their fitted values")
    print("   m[GeV]  stat%  syst% | heavy frac | worst-case shift")
    for target in PROBE_MASSES:
        index = int(np.argmin(np.abs(sqrt_s - target)))
        shift = 100 * worst_case_shift(masses, heavy_fraction, index)
        print(f" {sqrt_s[index]:7.4f} {100 * sigma_stat[index] / sigma[index]:6.1f}"
              f" {100 * sigma_syst[index] / sigma[index]:6.1f} |"
              f" {100 * heavy_fraction[index]:8.1f}%  |{shift:8.2f}%")

    clean = np.flatnonzero(heavy_fraction < CLEAN_THRESHOLD)
    print("\n-- 5. Verdict")
    print(f"  minimum heavy-state weight fraction anywhere in the table: "
          f"{100 * heavy_fraction.min():.1f}% at {sqrt_s[int(heavy_fraction.argmin())]:.4f} GeV")
    print(f"  maximum phi-tail weight fraction anywhere in the table: "
          f"{100 * fraction['phi'].max():.1f}%")
    if clean.size == 0:
        print(f"  bins with heavy-state weight below {100 * CLEAN_THRESHOLD:.0f}%: NONE")
    else:
        print(f"  bins below threshold: {clean.size}, up to {sqrt_s[clean[-1]]:.4f} GeV")

    assert heavy_fraction.min() > CLEAN_THRESHOLD, "a clean window exists; verdict text is stale"
    print("\n  No mass window isolates the one-multiplet baseline: the heavier states own")
    print("  at least {:.0f}% of the coherent sum everywhere, and their Lorentzian floors".format(
        100 * heavy_fraction.min()))
    print("  keep them from decaying away. A restricted-window fit would be dominated by the")
    print("  parameters supplied for omega' and omega'', not by the c_i and d_i being fitted.")


if __name__ == "__main__":
    main()