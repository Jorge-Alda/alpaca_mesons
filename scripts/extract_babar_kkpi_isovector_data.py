#!/usr/bin/env python3
"""Transcribe BaBar:2007ceh Table 7 below 1.78 GeV.

Run from the repository root:
    .venv/bin/python scripts/extract_babar_kkpi_isovector_data.py

The tabulated Dalitz-fit isovector K*(892)K component is retained in nb.
This script supplies data only and performs no CVC conversion or fit.
"""

import csv
import json
from pathlib import Path


# Consecutive 20 MeV bins beginning at 1.40 GeV, Table 7, columns 1–4.
CROSS_SECTIONS_NB = (
    ("0.32", "0.67"), ("0.68", "0.40"), ("1.40", "0.43"),
    ("1.27", "0.39"), ("2.33", "0.50"), ("2.96", "0.58"),
    ("2.83", "0.52"), ("3.08", "0.55"), ("2.05", "0.42"),
    ("1.67", "0.33"), ("3.46", "0.68"), ("1.52", "0.38"),
    ("1.83", "0.37"), ("2.82", "0.50"), ("2.10", "0.36"),
    ("1.28", "0.30"), ("1.53", "0.35"), ("1.60", "0.32"),
    ("1.44", "0.32"),
)


def write_data(output_directory: Path) -> tuple[Path, Path]:
    """Write the published bins and explicit uncertainty qualifications."""
    output_directory.mkdir(parents=True, exist_ok=True)
    csv_path = output_directory / "babar_2007_isovector_kstar_k_tau_region.csv"
    with csv_path.open("w", newline="") as output_file:
        writer = csv.writer(output_file, lineterminator="\n")
        writer.writerow(("energy_low_GeV", "energy_high_GeV", "cross_section_nb",
                         "statistical_error_nb"))
        for bin_number, (cross_section, error) in enumerate(CROSS_SECTIONS_NB):
            writer.writerow((f"{(140 + 2 * bin_number) / 100:.2f}",
                             f"{(142 + 2 * bin_number) / 100:.2f}",
                             cross_section, error))
    metadata = {
        "bibtex_key": "BaBar:2007ceh",
        "arxiv_id": "0710.4451v3",
        "source_url": "https://arxiv.org/html/0710.4451v3",
        "source_table": 7,
        "source_section": "VI.1",
        "observable": "isovector e+e- -> K*(892)K cross section from Dalitz fit",
        "energy_unit": "GeV", "cross_section_unit": "nb",
        "number_of_bins": len(CROSS_SECTIONS_NB),
        "selection": "All published bins contained in [1.40, 1.78] GeV",
        "bin_interpretation": "Published bin estimates, not point values or bin integrals",
        "statistical_errors": "Table 7 quoted errors; statistical covariance unavailable",
        "dalitz_fit_systematic_relative": 0.15,
        "dalitz_fit_systematic_source": "Paragraph following Tables 6 and 7",
        "dalitz_fit_systematic_relation_to_table_1": "Additional to the Table 1 systematic uncertainties",
        "dalitz_fit_systematic_correlation": {
            "source_statement": "No bin correlation matrix or correlation prescription supplied",
            "optional_model": "One common multiplicative nuisance with standard deviation 0.15",
            "optional_model_status": "Assumption, equivalent to full bin correlation; not a measured covariance",
        },
        "table_1_systematics": {
            "total_relative": 0.051,
            "components_relative": {
                "ISR luminosity": 0.030,
                "qqbar background": 0.011,
                "ISR background": 0.005,
                "track reconstruction efficiency": 0.010,
                "PID efficiency": 0.020,
                "KS reconstruction efficiency": 0.011,
                "KS to pi+pi- branching fraction": 0.007,
                "MC model": 0.030,
            },
            "correlations": "Unavailable; no correlation prescription introduced",
            "treatment": "Recorded separately; not included in the CSV statistical errors",
            "published_corrections": "The +4.2% total correction is already applied; do not apply it again",
        },
        "kinematic_qualifications": [
            "K*(892)K nominal threshold lies below the first 1.40 GeV bin; finite K* width belongs to the amplitude model",
            "The 1.76-1.78 GeV bin straddles the physical tau endpoint; do not truncate or treat it as a full tau-accessible bin",
            "No numerical tau mass or CVC/isospin normalization has been introduced",
        ],
        "fit_qualifications": [
            "Dalitz-model-dependent K*(892)K component, not the inclusive KKpi spectrum",
            "No RChT coefficient extraction or higher-resonance subtraction performed",
        ],
    }
    json_path = csv_path.with_suffix(".json")
    json_path.write_text(json.dumps(metadata, indent=2) + "\n")
    return csv_path, json_path


if __name__ == "__main__":
    repository_root = Path(__file__).resolve().parents[1]
    for written_path in write_data(repository_root / "data" / "vector_sector"):
        print(written_path.relative_to(repository_root))
