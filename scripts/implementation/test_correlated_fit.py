"""Tests for fixed-covariance Gaussian least squares.

Invocation:
    .venv/bin/python -m unittest scripts.implementation.test_correlated_fit
"""

import math
import unittest
import warnings
from dataclasses import FrozenInstanceError

import numpy as np

from scripts.implementation.correlated_fit import (
    FitProblem,
    Nuisance,
    Observation,
    OPEConstraint,
    Parameter,
    apply_nuisances,
    fit_correlated,
)


def affine(parameters, config):
    x = np.asarray(config["x"], dtype=float)
    return parameters["intercept"] + parameters["slope"] * x


def single(parameter):
    def predict(parameters, config):
        return np.array([parameters[parameter]])

    return predict


def sum_of_two(parameters, config):
    return np.array([parameters["intercept"] + parameters["slope"]])


class CorrelatedFitTests(unittest.TestCase):
    def test_correlated_whitening_and_component_chi2(self):
        covariance = np.array([[1.0, 0.2, 0.0], [0.2, 2.0, 0.3], [0.0, 0.3, 1.5]])
        config = {"x": [0.0, 1.0, 2.0]}
        problem = FitProblem(
            (Parameter("intercept", 0.0), Parameter("slope", 0.0)),
            (Observation([1.0, 3.0, 6.0], covariance, affine, config),),
        )
        result = fit_correlated(problem)
        self.assertTrue(result.success)
        expected = np.linalg.solve(
            np.linalg.cholesky(covariance), affine(result.as_mapping(), config) - np.array([1.0, 3.0, 6.0])
        )
        self.assertGreater(np.linalg.norm(expected), 0.1)
        np.testing.assert_allclose(result.residuals, expected, rtol=2e-8, atol=2e-10)
        self.assertAlmostEqual(result.chi2.data, expected @ expected, places=12)
        self.assertEqual(result.chi2.nuisances, 0.0)
        self.assertEqual(result.chi2.theory, 0.0)
        self.assertAlmostEqual(result.chi2.total, result.chi2.data, places=12)
        self.assertEqual(result.data_rank.rank, 2)
        self.assertEqual(result.full_rank.rank, 2)
        self.assertEqual(result.data_rank.null_directions.shape, (0, 2))

    def test_data_rank_deficiency_null_directions_and_chi2(self):
        problem = FitProblem(
            (Parameter("intercept", 0.0), Parameter("slope", 0.0)),
            (Observation([3.0], [[4.0]], sum_of_two),),
            ope_constraints=(OPEConstraint(("intercept", "slope"), [[1.0, -1.0]], [0.0], [[9.0]], ("toy OPE",)),),
            block_independence="data and fixed OPE theory blocks are independent",
        )
        result = fit_correlated(problem)
        self.assertEqual(result.data_rank.rank, 1)
        self.assertEqual(result.data_rank.null_directions.shape, (1, 2))
        data_jacobian = np.array([[0.5, 0.5]])
        for direction in result.data_rank.null_directions:
            np.testing.assert_allclose(data_jacobian @ direction, 0.0, atol=2e-8)
        self.assertEqual(result.full_rank.rank, 2)
        self.assertEqual(result.full_rank.null_directions.shape, (0, 2))
        self.assertAlmostEqual(
            result.chi2.total, result.chi2.data + result.chi2.nuisances + result.chi2.theory, places=14
        )
        np.testing.assert_allclose(result.residuals @ result.residuals, result.chi2.total)
        self.assertEqual(problem.ope_design_rank, 1)
        self.assertEqual(problem.ope_rank_deficiency, 1)

    def test_known_singular_values_and_condition_number(self):
        problem = FitProblem(
            (Parameter("x", 0.0), Parameter("y", 0.0)),
            (Observation([1.0, 6.0], np.eye(2), lambda p, config: np.array([p["x"], 3.0 * p["y"]])),),
        )
        result = fit_correlated(problem)
        np.testing.assert_allclose(result.data_rank.singular_values, [3.0, 1.0], rtol=2e-8)
        self.assertAlmostEqual(result.data_rank.condition_number, 3.0, places=7)
        np.testing.assert_allclose(result.as_mapping()["y"], 2.0, atol=1e-9)

    def test_active_upper_bound(self):
        problem = FitProblem(
            (Parameter("x", 0.0, lower=-1.0, upper=1.0),),
            (Observation([3.0], [[1.0]], single("x")),),
        )
        result = fit_correlated(problem)
        np.testing.assert_allclose(result.values, [1.0], atol=1e-9)
        self.assertEqual(result.active_bounds, ("x:upper",))
        self.assertAlmostEqual(result.chi2.total, 4.0, places=10)

    def test_ope_block_residual_matches_explicit_cholesky_solve(self):
        constraint = OPEConstraint(("intercept", "slope"), [[1.0, -1.0]], [0.25], [[9.0]])
        problem = FitProblem(
            (Parameter("intercept", 0.0), Parameter("slope", 0.0)),
            (Observation([0.0], [[1.0]], lambda p, config: np.array([0.0])),),
            ope_constraints=(constraint,),
            block_independence="the data and OPE theory blocks have zero cross-covariance",
        )
        result = fit_correlated(problem)
        fitted = result.as_mapping()
        residual = constraint.matrix @ np.array([fitted["intercept"], fitted["slope"]]) - constraint.target
        expected = np.linalg.solve(np.linalg.cholesky(constraint.covariance), residual)
        np.testing.assert_allclose(result.residuals[-1:], expected, rtol=1e-8, atol=1e-12)
        self.assertAlmostEqual(result.chi2.theory, expected @ expected, places=10)
        self.assertAlmostEqual(result.chi2.data, 0.0, places=12)
        self.assertEqual(constraint.labels, ("ope-0",))

    def test_frozen_inputs_and_validation(self):
        parameter = Parameter("x", 0.0)
        with self.assertRaises(FrozenInstanceError):
            parameter.start = 2.0
        with self.assertRaisesRegex(ValueError, "finite"):
            Observation([np.nan], [[1.0]], single("x"))
        with self.assertRaisesRegex(ValueError, "finite"):
            Observation([np.inf], [[1.0]], single("x"))
        with self.assertRaisesRegex(ValueError, "symmetric"):
            Observation([1.0, 2.0], [[1.0, 0.2], [0.1, 1.0]], single("x"))
        with self.assertRaisesRegex(ValueError, "positive definite"):
            Observation([1.0, 2.0], [[1.0, 1.0], [1.0, 1.0]], single("x"))
        with self.assertRaisesRegex(ValueError, "shape"):
            Observation([1.0, 2.0], [[1.0]], single("x"))
        with self.assertRaisesRegex(ValueError, "within"):
            Parameter("x", 2.0, upper=1.0)
        with self.assertRaisesRegex(ValueError, "smaller than upper bound"):
            Parameter("x", 0.0, lower=1.0, upper=1.0)
        with self.assertRaisesRegex(ValueError, "must not be empty"):
            Observation([], np.empty((0, 0)), single("x"))
        with self.assertRaisesRegex(ValueError, "unique"):
            FitProblem(
                (Parameter("x", 0.0), Parameter("x", 0.0)),
                (Observation([1.0], [[1.0]], single("x")),),
            )

    def test_complex_input_is_rejected_without_warnings(self):
        with warnings.catch_warnings(record=True) as caught:
            with self.assertRaisesRegex(ValueError, "complex values"):
                Observation(np.array([1.0 + 1.0j]), [[1.0]], single("x"))
            with self.assertRaisesRegex(ValueError, "complex values"):
                Observation([1.0], np.array([[1.0 + 1.0j]]), single("x"))
            with self.assertRaisesRegex(ValueError, "complex values"):
                OPEConstraint(("x",), np.array([[1.0 + 1.0j]]), [0.0], [[1.0]])
            with self.assertRaisesRegex(ValueError, "complex values"):
                Nuisance("eta", shift=np.array([1.0 + 1.0j]), unit="nb")
            with self.assertRaisesRegex(ValueError, "complex values"):
                apply_nuisances(np.array([1.0 + 1.0j]), (), {})
            with self.assertRaisesRegex(ValueError, "finite and real"):
                Parameter("x", 1.0 + 1.0j)
        self.assertEqual(caught, [])

        problem = FitProblem(
            (Parameter("x", 0.0),),
            (Observation([0.0], [[1.0]], lambda p, config: np.array([p["x"] + 1.0j])),),
        )
        with warnings.catch_warnings(record=True) as caught:
            with self.assertRaisesRegex(ValueError, "complex values"):
                fit_correlated(problem)
        self.assertEqual(caught, [])

    def test_predictor_output_validation(self):
        nonfinite = FitProblem(
            (Parameter("x", 0.0),), (Observation([0.0], [[1.0]], lambda p, config: np.array([np.inf])),),
        )
        with self.assertRaisesRegex(ValueError, "finite"):
            fit_correlated(nonfinite)
        wrong_shape = FitProblem(
            (Parameter("x", 0.0),),
            (Observation([0.0, 0.0], np.eye(2), lambda p, config: np.array([p["x"]])),),
        )
        with self.assertRaisesRegex(ValueError, "shape"):
            fit_correlated(wrong_shape)
        with self.assertRaises(TypeError):
            fit_correlated((Parameter("x", 0.0), (Observation([0.0], [[1.0]], single("x")),)))
        with self.assertRaises(TypeError):
            Observation([1.0], [[1.0]], "not callable")

    def test_multiple_blocks_require_independence_declaration(self):
        parameters = (Parameter("x", 0.0),)
        observations = (
            Observation([0.0], [[1.0]], single("x")),
            Observation([0.0], [[1.0]], single("x")),
        )
        with self.assertRaisesRegex(ValueError, "block_independence"):
            FitProblem(parameters, observations)
        with self.assertRaisesRegex(ValueError, "nonempty"):
            FitProblem(parameters, observations, block_independence="   ")
        with self.assertRaisesRegex(ValueError, "block_independence"):
            FitProblem(
                parameters + (Parameter("eta", 0.0),),
                (Observation([0.0], [[1.0]], single("x"), nuisances=(Nuisance("eta", fraction=0.1),)),),
            )
        self.assertIsNone(FitProblem(parameters, observations[:1]).block_independence)
        # a single block raises no cross-covariance question, so a statement is tolerated
        self.assertEqual(
            FitProblem(
                parameters, observations[:1], block_independence="a single block needs no assumption"
            ).block_independence,
            "a single block needs no assumption",
        )

    def test_ope_and_nuisance_names_must_be_declared(self):
        valid = OPEConstraint(("intercept", "slope"), [[1.0, -1.0]], [0.0], [[9.0]])
        self.assertEqual(
            FitProblem(
                (Parameter("intercept", 0.0), Parameter("slope", 0.0)),
                (Observation([3.0], [[4.0]], sum_of_two),),
                ope_constraints=(valid,),
                block_independence="the data and OPE theory blocks are independent",
            ).ope_design_rank,
            1,
        )
        with self.assertRaisesRegex(ValueError, "occur in parameters"):
            FitProblem(
                (Parameter("intercept", 0.0), Parameter("slope", 0.0)),
                (Observation([3.0], [[4.0]], sum_of_two),),
                ope_constraints=(OPEConstraint(("intercept", "absent"), [[1.0, 1.0]], [0.0], [[9.0]]),),
                block_independence="the data and OPE theory blocks are independent",
            )

    def test_standardized_multiplicative_nuisance_response_and_manual_chi2(self):
        nuisance = Nuisance("eta_norm", fraction=0.1)
        base = np.array([10.0, 20.0])
        np.testing.assert_array_equal(apply_nuisances(base, (nuisance,), {"eta_norm": 0.0}), base)
        np.testing.assert_allclose(apply_nuisances(base, (nuisance,), {"eta_norm": 1.0}), 1.1 * base, rtol=2e-15)
        self.assertTrue(np.all(apply_nuisances(base, (nuisance,), {"eta_norm": -100.0}) > 0.0))

        problem = FitProblem(
            (Parameter("eta_norm", 0.0),),
            (Observation([12.0], [[4.0]], lambda p, config: np.array([10.0]), nuisances=(nuisance,)),),
            block_independence="the data and standardized nuisance prior blocks are independent",
        )
        result = fit_correlated(problem)
        eta = result.as_mapping()["eta_norm"]
        prediction = 10.0 * math.exp(math.log1p(0.1) * eta)
        self.assertAlmostEqual(result.chi2.total, (prediction - 12.0) ** 2 / 4.0 + eta**2, places=12)
        self.assertAlmostEqual(result.chi2.nuisances, eta**2, places=12)
        self.assertAlmostEqual(result.chi2.data, (prediction - 12.0) ** 2 / 4.0, places=12)
        self.assertEqual(result.residuals.size, 2)

    def test_nuisance_linear_limit_and_additive_response(self):
        eta = 0.7
        errors = []
        for fraction in (1e-2, 5e-3):
            exact = apply_nuisances(np.array([3.0]), (Nuisance("eta", fraction=fraction),), {"eta": eta})[0]
            errors.append(abs(exact - 3.0 * (1.0 + fraction * eta)))
        self.assertAlmostEqual(errors[0] / errors[1], 4.0, delta=0.05)
        # (1 + fraction) ** eta = 1 + fraction * eta + eta (eta - 1) / 2 * fraction**2 + O(fraction**3)
        fraction = 5e-3
        exact = apply_nuisances(np.array([3.0]), (Nuisance("eta", fraction=fraction),), {"eta": eta})[0]
        self.assertAlmostEqual(
            exact - 3.0 * (1.0 + fraction * eta), 1.5 * (eta**2 - eta) * fraction**2, delta=1e-7
        )
        additive = Nuisance("eta_offset", shift=[1.0, -2.0], unit="nb")
        np.testing.assert_array_equal(
            apply_nuisances(np.array([4.0, 5.0]), (additive,), {"eta_offset": 0.0}), [4.0, 5.0]
        )
        np.testing.assert_allclose(
            apply_nuisances(np.array([4.0, 5.0]), (additive,), {"eta_offset": 1.0}), [5.0, 3.0]
        )
        mixed = (additive, Nuisance("eta_norm", fraction=0.1))
        np.testing.assert_allclose(
            apply_nuisances(np.array([4.0, 5.0]), mixed, {"eta_offset": 1.0, "eta_norm": 1.0}),
            1.1 * np.array([4.0, 5.0]) + np.array([1.0, -2.0]),
            rtol=2e-15,
        )

    def test_nuisance_metadata_rejections(self):
        additive = Nuisance("eta", shift=[1.0], unit="nb")
        with self.assertRaisesRegex(ValueError, "length"):
            Observation([1.0, 2.0], np.eye(2), single("x"), unit="nb", nuisances=(additive,))
        with self.assertRaisesRegex(ValueError, "unit"):
            Observation([1.0], [[1.0]], single("x"), unit="GeV", nuisances=(additive,))
        with self.assertRaisesRegex(ValueError, "unit"):
            Observation([1.0], [[1.0]], single("x"), nuisances=(additive,))
        with self.assertRaisesRegex(ValueError, "nonnegative"):
            Nuisance("eta", fraction=-1.0)
        with self.assertRaisesRegex(ValueError, "unit"):
            Nuisance("eta", fraction=0.1, unit="nb")
        with self.assertRaisesRegex(ValueError, "exactly one of"):
            Nuisance("eta")
        with self.assertRaisesRegex(ValueError, "exactly one of"):
            Nuisance("eta", fraction=0.1, shift=[1.0])

        observation = Observation([1.0], [[1.0]], single("x"), unit="nb", nuisances=(additive,))
        with self.assertRaisesRegex(ValueError, "occur in parameters"):
            FitProblem(
                (Parameter("x", 0.0),),
                (observation,),
                block_independence="the data and nuisance blocks are independent",
            )
        duplicated = Observation(
            [1.0], [[1.0]], single("x"), unit="nb",
            nuisances=(Nuisance("eta", fraction=0.1), Nuisance("eta", shift=[1.0], unit="nb")),
        )
        with self.assertRaisesRegex(ValueError, "exactly one nuisance"):
            FitProblem(
                (Parameter("x", 0.0), Parameter("eta", 0.0)),
                (duplicated,),
                block_independence="the data and nuisance blocks are independent",
            )


if __name__ == "__main__":
    unittest.main()