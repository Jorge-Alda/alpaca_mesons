"""Tests for deterministic odd-parity data loading and first-release assembly."""

import hashlib
import json
from pathlib import Path
import shutil
import tempfile
from types import MappingProxyType
import unittest
from unittest.mock import patch

import numpy as np

from scripts.implementation.correlated_fit import PredictorRegistry
from scripts.implementation.odd_parity_data import (
    CANONICAL_COUPLING_ORDER,
    DATA_DIRECTORY,
    FILE_SHA256,
    build_first_release_fit_problem,
    build_fixed_ope_constraint,
    load_babar_3pi_born_spectrum,
    load_babar_pi0_table11,
    load_cleo_tff,
    load_kloe_phi_3pi_dalitz_fit,
    load_omega_pi0_gamma_observation,
    load_omega_pi0_gamma_table,
)
from scripts.implementation.odd_parity_rcht import OPEParameters, fixed_ope_prior, ope_linear_system
from scripts.implementation.radiative_odd_parity import (
    OMEGA_PI_GAMMA_BRANCHING_FRACTION_PREDICTOR_ID,
    OmegaPiGammaInputs,
    omega_pi_gamma_branching_fraction_predictor,
)


class OddParityDataTests(unittest.TestCase):
    def setUp(self):
        self.inputs = OmegaPiGammaInputs(0.3028, 0.092, 0.145, 0.770, 0.783, 0.135)
        self.ope = OPEParameters(3, 0.770, 0.145, 0.092)

    def test_all_table_counts_immutability_and_asymmetric_preservation(self):
        tables = (
            (load_cleo_tff(), 63),
            (load_babar_pi0_table11(), 17),
            (load_babar_3pi_born_spectrum(), 78),
            (load_kloe_phi_3pi_dalitz_fit(), 10),
            (load_omega_pi0_gamma_table(), 1),
        )
        for table, count in tables:
            self.assertEqual(table.rows.size, count)
            self.assertFalse(table.rows.flags.writeable)
            self.assertEqual(table.metadata["likelihood"]["schema_version"], "odd-parity-likelihood-v1")
            self.assertEqual(table.metadata["likelihood"]["asymmetric_error_rule"], "mean_absolute")
        cleo = tables[0][0]
        self.assertAlmostEqual(cleo.rows[0]["q2F_GeV"], 0.121)
        self.assertAlmostEqual(cleo.rows[0]["stat_GeV"], 0.008)
        babar = tables[1][0]
        self.assertEqual(babar.rows[-2]["q2F_err_low_MeV"], 40.0)
        self.assertEqual(babar.rows[-2]["q2F_err_high_MeV"], 33.0)
        self.assertEqual((babar.rows[-2]["q2F_err_low_MeV"] + babar.rows[-2]["q2F_err_high_MeV"]) / 2.0, 36.5)
        spectrum = tables[2][0]
        self.assertEqual(tuple(spectrum.rows.dtype.names), (
            "sqrt_s_GeV", "Ncorr", "Ncorr_stat", "Ncorr_syst", "efficiency_percent",
            "efficiency_err_percent", "dL_nb_inv", "sigma_nb", "sigma_stat_nb", "sigma_syst_nb"))
        self.assertEqual(tuple(spectrum.rows[0][name] for name in ("Ncorr", "Ncorr_stat", "Ncorr_syst")), (40.0, 13.0, 3.0))
        spectrum_likelihood = spectrum.metadata["likelihood"]
        self.assertEqual(spectrum_likelihood["diagonal_source_ids"], ("babar-2004-statistical",))
        self.assertEqual(
            spectrum_likelihood["nuisance_sources"][0],
            {"source_id": "babar-2004-mass-bin-systematic", "response": "additive",
             "response_column": "sigma_syst_nb", "unit": "nb", "predictor_connection": "absent"},
        )
        self.assertEqual(spectrum_likelihood["nuisance_sources"][1]["response"], "multiplicative")
        kloe = tables[3][0]
        self.assertEqual(kloe.rows[0]["value"], "")
        self.assertEqual(kloe.rows[1]["value"], "775.8")
        self.assertEqual(kloe.rows[2]["value"], "")
        nuisance = babar.metadata["likelihood"]["nuisance_sources"]
        self.assertEqual(tuple(item["fraction"] for item in nuisance), (0.03, 0.023))
        self.assertTrue(all(item["predictor_connection"] == "absent" for item in nuisance))

    def test_omega_observation_is_connected_with_diagonal_variance(self):
        total = 0.00849
        observation = load_omega_pi0_gamma_observation(self.inputs, total)
        np.testing.assert_array_equal(observation.values, [0.0835])
        np.testing.assert_allclose(observation.covariance, [[0.0027**2]])
        self.assertEqual(observation.unit, "fraction")
        self.assertEqual(observation.diagonal_source_ids, ("pdg-2024-omega-pi0-gamma-total",))
        self.assertEqual(observation.predictor_id, OMEGA_PI_GAMMA_BRANCHING_FRACTION_PREDICTOR_ID)
        self.assertEqual(observation.config["GammaOmegaTotal"], total)
        registry = PredictorRegistry()
        registry.register(OMEGA_PI_GAMMA_BRANCHING_FRACTION_PREDICTOR_ID,
                          omega_pi_gamma_branching_fraction_predictor)
        parameters = {name: 0.01 for name in CANONICAL_COUPLING_ORDER}
        self.assertEqual(registry.predict(observation.predictor_id, parameters, observation.config).shape, (1,))
        with self.assertRaisesRegex(ValueError, "positive"):
            load_omega_pi0_gamma_observation(self.inputs, 0.0)

    def test_ope_adapter_and_first_release_rank_deficiency(self):
        constraint = build_fixed_ope_constraint(self.ope)
        matrix, target = ope_linear_system(self.ope)
        np.testing.assert_array_equal(constraint.matrix, matrix)
        np.testing.assert_array_equal(constraint.target, target)
        np.testing.assert_array_equal(constraint.covariance, fixed_ope_prior(self.ope).covariance)
        self.assertEqual(constraint.parameter_names, CANONICAL_COUPLING_ORDER)

        starts = {name: 0.0 for name in CANONICAL_COUPLING_ORDER}
        bounds = {name: (-1.0, 1.0) for name in CANONICAL_COUPLING_ORDER}
        problem = build_first_release_fit_problem(self.inputs, 0.00849, self.ope, starts, bounds)
        self.assertEqual(len(problem.observations), 1)
        self.assertEqual(len(problem.ope_constraints), 1)
        self.assertIn("zero cross-covariance", problem.block_independence_assumption)
        self.assertEqual(problem.ope_design_rank, 5)
        self.assertEqual(problem.ope_rank_deficiency, 6)
        with self.assertRaisesRegex(ValueError, "exactly"):
            build_first_release_fit_problem(self.inputs, 0.00849, self.ope,
                                            {key: value for key, value in starts.items() if key != "c7"}, bounds)

    def test_hash_schema_nonfinite_and_unit_rejections(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            for suffix in ("csv", "metadata.json"):
                shutil.copy2(DATA_DIRECTORY / f"omega_pi0_gamma_pdg.{suffix}", directory)
            csv_path = directory / "omega_pi0_gamma_pdg.csv"
            csv_path.write_text(csv_path.read_text().replace("0.0835", "0.0836"))
            with self.assertRaisesRegex(ValueError, "SHA256 mismatch"):
                load_omega_pi0_gamma_table(directory)

            def accept_modified(path):
                hashes = dict(FILE_SHA256)
                hashes[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
                return MappingProxyType(hashes)

            csv_path.write_text("observable,value,stat_or_total_error,unit\n"
                                "omega_to_pi0_gamma_branching_fraction,nan,0.0027,fraction\n")
            with patch("scripts.implementation.odd_parity_data.FILE_SHA256", accept_modified(csv_path)):
                with self.assertRaisesRegex(ValueError, "finite"):
                    load_omega_pi0_gamma_table(directory)

            csv_path.write_text((DATA_DIRECTORY / "omega_pi0_gamma_pdg.csv").read_text())
            metadata_path = directory / "omega_pi0_gamma_pdg.metadata.json"
            metadata = json.loads(metadata_path.read_text())
            metadata["units"]["branching_fraction"] = "percent"
            metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
            with patch("scripts.implementation.odd_parity_data.FILE_SHA256", accept_modified(metadata_path)):
                with self.assertRaisesRegex(ValueError, "units mismatch"):
                    load_omega_pi0_gamma_table(directory)

    def test_additive_nuisance_metadata_rejections(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            for suffix in ("csv", "metadata.json"):
                shutil.copy2(DATA_DIRECTORY / f"babar_3pi_born_spectrum.{suffix}", directory)
            metadata_path = directory / "babar_3pi_born_spectrum.metadata.json"

            def accept_metadata():
                hashes = dict(FILE_SHA256)
                hashes[metadata_path.name] = hashlib.sha256(metadata_path.read_bytes()).hexdigest()
                return MappingProxyType(hashes)

            metadata = json.loads(metadata_path.read_text())
            metadata["likelihood"]["nuisance_sources"][0]["response_column"] = "sigma_nb"
            metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
            with patch("scripts.implementation.odd_parity_data.FILE_SHA256", accept_metadata()):
                with self.assertRaisesRegex(ValueError, "response column or unit"):
                    load_babar_3pi_born_spectrum(directory)

            metadata = json.loads((DATA_DIRECTORY / "babar_3pi_born_spectrum.metadata.json").read_text())
            metadata["likelihood"]["nuisance_sources"][0]["response"] = "multiplicative"
            metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
            with patch("scripts.implementation.odd_parity_data.FILE_SHA256", accept_metadata()):
                with self.assertRaisesRegex(ValueError, "multiplicative nuisance"):
                    load_babar_3pi_born_spectrum(directory)


if __name__ == "__main__":
    unittest.main()
