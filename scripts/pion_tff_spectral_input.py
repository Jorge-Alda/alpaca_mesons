"""Read explicitly supplied pion-TFF spectral inputs and calculate moments.

No physical input is bundled. Run ``.venv/bin/python
scripts/pion_tff_spectral_input.py`` to report the missing input contract.
The source formula is Hoferichter:2018kwz, Eq. (3.14). The supplied f1
uses that source's positive pion-TFF phase convention. The output changes
the external pion phase to the repository convention F_pi(0,0)<0; it does
not rescale physical input when the decay-constant convention changes.

Invariant grids have units GeV^2, F_pi^V is dimensionless, f1 has units
GeV^-3, and the resulting double-spectral density has units GeV^-1.
Composite Simpson quadrature is an explicit finite-grid approximation.
This module adds neither a tree rho pole nor an effective/asymptotic term.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import simpson


INPUT_CONTRACT = {
    "status": "physical_input_required",
    "metadata_fields": {
        "bibtex_keys": "nonempty list of source keys already in bibliography.bib",
        "source_convention": "Hoferichter:2018kwz-positive-pion-TFF",
        "input_decay_constant_GeV": "declared source normalization, not silently rescaled",
        "charged_pion_mass_GeV": "mass used in q_pi and input amplitudes",
        "isoscalar_threshold_GeV2": "declared lowest included isoscalar threshold",
        "spectral_file": "path to npz, relative to this metadata file",
        "provenance": "source URLs, fit version, normalization, and extraction method",
        "uncertainty_description": "correlated fit replicas/covariance and propagation prescription",
        "isovector_description": "F_pi^V model, phase input, polynomial and high-energy continuation",
        "isoscalar_description": "gamma*->3pi f1 normalization, subtractions, crossed cuts and fitted data",
        "matching_description": "spectral cutoffs; explicitly exclude duplicated tree rho strength",
    },
    "npz_arrays": {
        "isovector_GeV2": "increasing 1D x grid from 4 M_pi^2 through s_v",
        "isoscalar_GeV2": "increasing 1D y grid from declared threshold through s_s",
        "pion_vector_form_factor": "complex dimensionless array with shape (len(x),)",
        "three_pion_partial_wave_GeVminus3": "complex f1(x,y), shape (len(x), len(y))",
    },
    "limitations": [
        "A total e+e- -> 3pi cross section does not supply the complex f1(x,y).",
        "An elastic Omnes function alone does not fix the pion vector form factor or the three-pion amplitude.",
        "No moments or odd-coefficient constraint can be inferred until physical amplitudes and correlated uncertainties are supplied.",
    ],
}


def spectral_density(
    isovector_GeV2: np.ndarray,
    pion_vector_form_factor: np.ndarray,
    three_pion_partial_wave_GeVminus3: np.ndarray,
    charged_pion_mass_GeV: float,
) -> np.ndarray:
    """Return canonical signed rho(x,y) from source-normalized amplitudes.

    The minus sign converts the *complete* positive source scalar to the
    negative repository scalar. The electromagnetic factors must already
    follow Eq. (3.14); no additional e or e^2 factor is inserted.
    """
    isovector = np.asarray(isovector_GeV2, dtype=float)
    pion_momentum = np.sqrt(np.maximum(isovector / 4 - charged_pion_mass_GeV**2, 0))
    return -pion_momentum[:, None] ** 3 / (12 * np.pi * np.sqrt(isovector[:, None])) * np.imag(
        np.conjugate(pion_vector_form_factor)[:, None]
        * three_pion_partial_wave_GeVminus3
    )


def spectral_moments(
    isovector_GeV2: np.ndarray,
    isoscalar_GeV2: np.ndarray,
    density_GeVminus1: np.ndarray,
) -> dict[str, float]:
    """Integrate the four correlated moments of the Bose-symmetrized TFF.

    C_single is lim Q^2 F_disp(-Q^2,0) for this finite spectral rectangle.
    Both placements contribute, including when the integration axes differ.
    """
    isovector = np.asarray(isovector_GeV2)[:, None]
    isoscalar = np.asarray(isoscalar_GeV2)[None, :]

    def integrate(weight: np.ndarray) -> float:
        return float(simpson(simpson(density_GeVminus1 * weight, x=isoscalar_GeV2, axis=1), x=isovector_GeV2) / np.pi**2)

    return {
        "F_disp_0_0_GeVminus1": integrate(2 / (isovector * isoscalar)),
        "F_disp_derivative_1_GeVminus3": integrate(1 / (isovector**2 * isoscalar) + 1 / (isovector * isoscalar**2)),
        "F_disp_derivative_12_GeVminus5": integrate(2 / (isovector**2 * isoscalar**2)),
        "C_disp_single_GeV": integrate(1 / isovector + 1 / isoscalar),
    }


def load_input(metadata_path: Path) -> tuple[dict, dict[str, np.ndarray]]:
    """Require declared conventions and complex amplitudes over a full grid.

    Metadata declarations are preserved for audit; this does not certify
    source extraction, the adequacy of a fit, or its uncertainty propagation.
    """
    metadata = json.loads(metadata_path.read_text())
    missing = [name for name in INPUT_CONTRACT["metadata_fields"] if not metadata.get(name)]
    if missing:
        raise ValueError("missing metadata: " + ", ".join(missing))
    if metadata["source_convention"] != "Hoferichter:2018kwz-positive-pion-TFF":
        raise ValueError("source_convention must be Hoferichter:2018kwz-positive-pion-TFF")
    if not isinstance(metadata["bibtex_keys"], list):
        raise ValueError("bibtex_keys must be a nonempty list")
    bibliography = (Path(__file__).resolve().parents[1] / "bibliography.bib").read_text()
    for key in metadata["bibtex_keys"]:
        if "{" + key + "," not in bibliography:
            raise ValueError(f"source {key!r} is absent from bibliography.bib")
    with np.load(metadata_path.parent / metadata["spectral_file"], allow_pickle=False) as archive:
        arrays = {name: archive[name] for name in INPUT_CONTRACT["npz_arrays"]}
    isovector = arrays["isovector_GeV2"]
    isoscalar = arrays["isoscalar_GeV2"]
    for name, grid in (("isovector", isovector), ("isoscalar", isoscalar)):
        if grid.ndim != 1 or len(grid) < 3 or not np.all(np.isfinite(grid)) or not np.all(np.diff(grid) > 0):
            raise ValueError(f"{name} grid must be finite, increasing, and have at least three points")
    mass = float(metadata["charged_pion_mass_GeV"])
    if not np.isclose(isovector[0], 4 * mass**2, rtol=1e-10, atol=1e-12):
        raise ValueError("isovector grid must begin at 4 M_pi^2")
    if not np.isclose(isoscalar[0], metadata["isoscalar_threshold_GeV2"], rtol=1e-10, atol=1e-12):
        raise ValueError("isoscalar grid must begin at the declared threshold")
    for name, shape in (("pion_vector_form_factor", (len(isovector),)), ("three_pion_partial_wave_GeVminus3", (len(isovector), len(isoscalar)))):
        values = arrays[name]
        if values.shape != shape or not np.all(np.isfinite(values)):
            raise ValueError(f"{name} must have finite values and shape {shape}")
    return metadata, arrays


def main() -> int:
    """Print a missing-input contract or moments, without writing a fit result."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("metadata", nargs="?", type=Path)
    arguments = parser.parse_args()
    if arguments.metadata is None:
        print(json.dumps(INPUT_CONTRACT, indent=2))
        return 2
    metadata, arrays = load_input(arguments.metadata)
    density = spectral_density(arrays["isovector_GeV2"], arrays["pion_vector_form_factor"], arrays["three_pion_partial_wave_GeVminus3"], metadata["charged_pion_mass_GeV"])
    print(json.dumps({
        "status": "central_grid_diagnostic_only",
        "source_metadata": metadata,
        "quadrature": "composite Simpson; grid convergence and correlated uncertainty propagation remain required",
        "moments": spectral_moments(arrays["isovector_GeV2"], arrays["isoscalar_GeV2"], density),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
