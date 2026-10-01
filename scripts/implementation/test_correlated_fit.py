"""Tests for generic fixed-covariance nonlinear least squares.

Invocation:
    .venv/bin/python -m unittest scripts.implementation.test_correlated_fit
"""

import json
import math
import unittest
import warnings
from dataclasses import FrozenInstanceError

import numpy as np

from scripts.implementation.correlated_fit import (
    SCHEMA_VERSION,
    FitProblem,
    FixedOPEConstraint,
    NuisanceEffect,
    ObservationGroup,
    ParameterPrior,
    ParameterSpec,
    PredictorRegistry,
    apply_nuisance_effects,
    fit_correlated,
    fit_problem_digest,
    serialize_fit_problem,
)


def affine(parameters, config):
    x = np.asarray(config["x"], dtype=float)
    return parameters["intercept"] + parameters["slope"] * x


class CorrelatedFitTests(unittest.TestCase):
    def setUp(self):
        self.registry = PredictorRegistry()
        self.registry.register("toy.affine", affine)

    def test_correlated_whitening_and_component_chi2(self):
        covariance = np.array(
            [[1.0, 0.2, 0.0], [0.2, 2.0, 0.3], [0.0, 0.3, 1.5]]
        )
        problem = FitProblem(
            parameters=(ParameterSpec("intercept", 0.0), ParameterSpec("slope", 0.0)),
            observations=(ObservationGroup(
                [1.0, 3.0, 6.0],
                covariance,
                "toy.affine",
                {"x": [0.0, 1.0, 2.0]},
            ),),
        )
        result = fit_correlated(problem, self.registry)
        self.assertTrue(result.success)
        prediction = affine(result.as_mapping(), {"x": [0.0, 1.0, 2.0]})
        expected = np.linalg.solve(np.linalg.cholesky(covariance), prediction - [1.0, 3.0, 6.0])
        self.assertGreater(np.linalg.norm(expected), 0.1)
        np.testing.assert_allclose(result.residuals, expected, rtol=2e-8, atol=2e-10)
        np.testing.assert_allclose(result.diagnostics.chi2.data, (expected @ expected,))
        np.testing.assert_allclose(result.diagnostics.chi2.total, expected @ expected)
        self.assertEqual(result.diagnostics.data_jacobian.rank, 2)
        self.assertEqual(result.diagnostics.full_jacobian.rank, 2)
        self.assertEqual(result.problem_sha256, fit_problem_digest(problem))

    def test_priors_ope_rank_null_vectors_and_chi2_accounting(self):
        registry = PredictorRegistry()
        registry.register("toy.sum", lambda p, config: np.array([p["intercept"] + p["slope"]]))
        problem = FitProblem(
            parameters=(ParameterSpec("intercept", 0.0), ParameterSpec("slope", 0.0)),
            observations=(ObservationGroup([3.0], [[4.0]], "toy.sum"),),
            priors=(ParameterPrior(("intercept",), [1.0], [[1.0]]),),
            ope_constraints=(FixedOPEConstraint(
                ("intercept", "slope"), [[1.0, -1.0]], [0.0], [[9.0]], ("toy OPE",)
            ),),
            block_independence_assumption="data, parameter prior, and OPE blocks are independent",
        )
        result = fit_correlated(problem, registry)
        self.assertEqual(result.diagnostics.data_jacobian.rank, 1)
        self.assertEqual(result.diagnostics.data_jacobian.null_vectors.shape, (1, 2))
        data_jacobian = np.array([[0.5, 0.5]])
        for null_vector in result.diagnostics.data_jacobian.null_vectors:
            np.testing.assert_allclose(data_jacobian @ null_vector, 0.0, atol=2e-8)
        self.assertEqual(result.diagnostics.full_jacobian.rank, 2)
        self.assertEqual(result.diagnostics.full_jacobian.null_vectors.shape, (0, 2))
        pieces = result.diagnostics.chi2
        self.assertAlmostEqual(pieces.total, sum(pieces.data + pieces.priors + pieces.ope), places=14)
        np.testing.assert_allclose(result.residuals @ result.residuals, pieces.total)

    def test_known_jacobian_condition_number(self):
        registry = PredictorRegistry()
        registry.register("toy.scaled", lambda p, config: np.array([p["x"], 3.0 * p["y"]]))
        problem = FitProblem(
            (ParameterSpec("x", 0.0), ParameterSpec("y", 0.0)),
            (ObservationGroup([1.0, 6.0], np.eye(2), "toy.scaled"),),
        )
        result = fit_correlated(problem, registry)
        np.testing.assert_allclose(result.diagnostics.data_jacobian.singular_values, [3.0, 1.0], rtol=2e-8)
        self.assertAlmostEqual(result.diagnostics.data_jacobian.condition_number, 3.0, places=7)

    def test_active_bound(self):
        registry = PredictorRegistry()
        registry.register("toy.value", lambda p, config: np.array([p["x"]]))
        problem = FitProblem(
            (ParameterSpec("x", 0.0, lower=-1.0, upper=1.0),),
            (ObservationGroup([3.0], [[1.0]], "toy.value"),),
        )
        result = fit_correlated(problem, registry)
        np.testing.assert_allclose(result.values, [1.0], atol=1e-9)
        self.assertEqual(result.diagnostics.active_bounds, ("x:upper",))

    def test_canonical_serialization_contains_ids_not_callables(self):
        assumption = "the two toy observation blocks are independent"
        first = FitProblem(
            (ParameterSpec("intercept", 0.0), ParameterSpec("slope", 1.0)),
            (ObservationGroup([1.0], [[2.0]], "toy.affine", {"z": 2, "a": [1, True]}),
             ObservationGroup([2.0], [[3.0]], "toy.affine", {"x": [1.0]})),
            block_independence_assumption=assumption,
        )
        second = FitProblem(
            first.parameters,
            (ObservationGroup([1.0], [[2.0]], "toy.affine", {"a": [1, True], "z": 2}),
             ObservationGroup([2.0], [[3.0]], "toy.affine", {"x": [1.0]})),
            block_independence_assumption=assumption,
        )
        serialized = serialize_fit_problem(first)
        self.assertEqual(serialized, serialize_fit_problem(second))
        self.assertEqual(first.canonical_json(), serialized)
        self.assertEqual(first.sha256(), fit_problem_digest(second))
        decoded = json.loads(serialized)
        self.assertEqual(decoded["schema_version"], SCHEMA_VERSION)
        self.assertEqual(decoded["block_independence_assumption"], assumption)
        self.assertEqual(decoded["observations"][0]["predictor_id"], "toy.affine")
        self.assertNotIn("callable", serialized)

    def test_frozen_inputs_and_validation(self):
        parameter = ParameterSpec("x", 0.0)
        with self.assertRaises(FrozenInstanceError):
            parameter.initial = 2.0
        with self.assertRaisesRegex(ValueError, "finite"):
            ObservationGroup([np.nan], [[1.0]], "toy")
        with self.assertRaisesRegex(ValueError, "symmetric"):
            ObservationGroup([1.0, 2.0], [[1.0, 0.2], [0.1, 1.0]], "toy")
        with self.assertRaisesRegex(ValueError, "positive definite"):
            ObservationGroup([1.0, 2.0], [[1.0, 1.0], [1.0, 1.0]], "toy")
        with self.assertRaisesRegex(ValueError, "JSON"):
            ObservationGroup([1.0], [[1.0]], "toy", {"bad": object()})
        with self.assertRaisesRegex(ValueError, "within"):
            ParameterSpec("x", 2.0, upper=1.0)
        with self.assertRaisesRegex(ValueError, "must not be empty"):
            ObservationGroup([], np.empty((0, 0)), "toy")

    def test_complex_arrays_are_rejected_without_warnings(self):
        with warnings.catch_warnings(record=True) as caught:
            with self.assertRaisesRegex(ValueError, "complex values"):
                ObservationGroup(np.array([1.0 + 1.0j]), [[1.0]], "toy")
            with self.assertRaisesRegex(ValueError, "complex values"):
                ObservationGroup(np.array([1.0 + 1.0j], dtype=object), [[1.0]], "toy")
            with self.assertRaisesRegex(ValueError, "complex values"):
                ObservationGroup([1.0], np.array([[1.0 + 1.0j]]), "toy")
            with self.assertRaisesRegex(ValueError, "complex values"):
                ParameterPrior(("x",), np.array([1.0 + 1.0j]), [[1.0]])
            with self.assertRaisesRegex(ValueError, "complex values"):
                FixedOPEConstraint(("x",), np.array([[1.0 + 1.0j]]), [0.0], [[1.0]])
        self.assertEqual(caught, [])

        problem = FitProblem(
            (ParameterSpec("x", 0.0),),
            (ObservationGroup([0.0], [[1.0]], "complex"),),
        )
        registry = PredictorRegistry()
        registry.register("complex", lambda p, config: np.array([p["x"] + 1.0j]))
        with warnings.catch_warnings(record=True) as caught:
            with self.assertRaisesRegex(ValueError, "complex values"):
                fit_correlated(problem, registry)
        self.assertEqual(caught, [])

    def test_multiple_blocks_require_independence_declaration(self):
        parameters = (ParameterSpec("x", 0.0),)
        observations = (
            ObservationGroup([0.0], [[1.0]], "toy"),
            ObservationGroup([0.0], [[1.0]], "toy"),
        )
        with self.assertRaisesRegex(ValueError, "block_independence_assumption is required"):
            FitProblem(parameters, observations)
        with self.assertRaisesRegex(ValueError, "nonempty"):
            FitProblem(parameters, observations, block_independence_assumption="   ")
        with self.assertRaisesRegex(ValueError, "block_independence_assumption is required"):
            FitProblem(
                parameters,
                observations[:1],
                priors=(ParameterPrior(("x",), [0.0], [[1.0]]),),
            )

    def test_registry_and_prediction_validation(self):
        problem = FitProblem((ParameterSpec("x", 0.0),), (ObservationGroup([0.0], [[1.0]], "missing"),))
        with self.assertRaisesRegex(KeyError, "not registered"):
            fit_correlated(problem, PredictorRegistry())
        bad = PredictorRegistry()
        bad.register("missing", lambda p, config: np.array([np.inf]))
        with self.assertRaisesRegex(ValueError, "finite"):
            fit_correlated(problem, bad)
        duplicate = PredictorRegistry()
        duplicate.register("x", lambda p, config: [0.0])
        with self.assertRaisesRegex(ValueError, "already registered"):
            duplicate.register("x", lambda p, config: [0.0])

    def test_standardized_multiplicative_nuisance_response_and_manual_chi2(self):
        effect = NuisanceEffect("normalization", "eta_norm", "multiplicative", fraction=0.1)
        base = np.array([10.0, 20.0])
        np.testing.assert_array_equal(apply_nuisance_effects(base, (effect,), {"eta_norm": 0.0}), base)
        at_one = apply_nuisance_effects(base, (effect,), {"eta_norm": 1.0})
        np.testing.assert_allclose(at_one, 1.1 * base, rtol=2e-15)
        self.assertTrue(np.all(apply_nuisance_effects(base, (effect,), {"eta_norm": -100.0}) > 0.0))

        registry = PredictorRegistry()
        registry.register("toy.constant", lambda parameters, config: np.array([10.0]))
        problem = FitProblem(
            (ParameterSpec("eta_norm", 0.0),),
            (ObservationGroup([12.0], [[4.0]], "toy.constant", nuisance_effects=(effect,)),),
            priors=(ParameterPrior(("eta_norm",), [0.0], [[1.0]]),),
            block_independence_assumption="data and standardized nuisance prior are independent",
        )
        result = fit_correlated(problem, registry)
        eta = result.as_mapping()["eta_norm"]
        prediction = 10.0 * math.exp(math.log1p(0.1) * eta)
        manual = (prediction - 12.0) ** 2 / 4.0 + eta**2
        self.assertAlmostEqual(result.diagnostics.chi2.total, manual, places=12)

    def test_nuisance_linear_limit_and_additive_response(self):
        eta = 0.7
        errors = []
        for fraction in (1e-2, 5e-3):
            exact = apply_nuisance_effects(
                np.array([3.0]),
                (NuisanceEffect("n", "eta", "multiplicative", fraction=fraction),),
                {"eta": eta},
            )[0]
            linear = 3.0 * (1.0 + fraction * eta)
            errors.append(abs(exact - linear))
        self.assertAlmostEqual(errors[0] / errors[1], 4.0, delta=0.05)
        additive = NuisanceEffect("offset", "eta_offset", "additive", additive_vector=[1.0, -2.0], unit="nb")
        np.testing.assert_allclose(
            apply_nuisance_effects(np.array([4.0, 5.0]), (additive,), {"eta_offset": 1.0}),
            [5.0, 3.0],
        )

    def test_nuisance_metadata_rejections(self):
        additive = NuisanceEffect("offset", "eta", "additive", additive_vector=[1.0], unit="nb")
        with self.assertRaisesRegex(ValueError, "length"):
            ObservationGroup([1.0, 2.0], np.eye(2), "toy", unit="nb", nuisance_effects=(additive,))
        with self.assertRaisesRegex(ValueError, "unit"):
            ObservationGroup([1.0], [[1.0]], "toy", unit="GeV", nuisance_effects=(additive,))
        with self.assertRaisesRegex(ValueError, "both diagonal and nuisance"):
            ObservationGroup([1.0], [[1.0]], "toy", unit="nb", diagonal_source_ids=("offset",), nuisance_effects=(additive,))
        with self.assertRaisesRegex(ValueError, "nonnegative"):
            NuisanceEffect("bad", "eta", "multiplicative", fraction=-1.0)

        observation = ObservationGroup([1.0], [[1.0]], "toy", unit="nb", nuisance_effects=(additive,))
        with self.assertRaisesRegex(ValueError, "occur in parameters"):
            FitProblem((ParameterSpec("x", 0.0),), (observation,),
                       priors=(ParameterPrior(("x",), [0.0], [[1.0]]),),
                       block_independence_assumption="blocks independent")
        with self.assertRaisesRegex(ValueError, "exactly one prior"):
            FitProblem((ParameterSpec("eta", 0.0),), (observation,),
                       block_independence_assumption="blocks independent")
        with self.assertRaisesRegex(ValueError, "standardized"):
            FitProblem((ParameterSpec("eta", 0.0),), (observation,),
                       priors=(ParameterPrior(("eta",), [0.0], [[2.0]]),),
                       block_independence_assumption="blocks independent")

        other = NuisanceEffect("offset", "eta_other", "additive", additive_vector=[1.0], unit="nb")
        with self.assertRaisesRegex(ValueError, "consistent parameter"):
            FitProblem(
                (ParameterSpec("eta", 0.0), ParameterSpec("eta_other", 0.0)),
                (observation, ObservationGroup([1.0], [[1.0]], "toy", unit="nb", nuisance_effects=(other,))),
                priors=(ParameterPrior(("eta", "eta_other"), [0.0, 0.0], np.eye(2)),),
                block_independence_assumption="all blocks independent",
            )
        multiplicative = NuisanceEffect("global", "eta", "multiplicative", fraction=0.1)
        with self.assertRaisesRegex(ValueError, "both diagonal and nuisance"):
            FitProblem(
                (ParameterSpec("eta", 0.0),),
                (ObservationGroup([1.0], [[1.0]], "toy", diagonal_source_ids=("global",)),
                 ObservationGroup([1.0], [[1.0]], "toy", nuisance_effects=(multiplicative,))),
                priors=(ParameterPrior(("eta",), [0.0], [[1.0]]),),
                block_independence_assumption="all blocks independent",
            )


if __name__ == "__main__":
    unittest.main()
