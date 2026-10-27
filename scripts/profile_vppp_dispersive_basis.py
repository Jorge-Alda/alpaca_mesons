"""Conditional local VPPP identifiability with supplied finite-width spectra.

Run: .venv/bin/python scripts/profile_vppp_dispersive_basis.py INPUT.npz
Tests: .venv/bin/python -m unittest discover -s scripts -p profile_vppp_dispersive_basis.py

No physical inputs, widths, covariance, priors or amplitude data are provided.
Natural units GeV. Input spectral_s is in GeV^2; spectral_weights integrate
in that variable and include the caller's declared tail prescription. The
spectral grid must span all retained strength. Im resonance basis has units
GeV^-2 (unit residue); Im conformal basis is dimensionless per coefficient
of units GeV^-3; subtraction_map has units GeV^-3. residue_map has units GeV^-1. These maps differentiate with respect
to the literal combinations (g1-g2-g3), (2g4+g5).

NPZ fields: spectral_s (n,), spectral_weights (n,), resonance_spectral
(n,r), conformal_spectral (n,c), evaluation_q_squared (k,) complex off the
real cut, subtraction_map (2,), residue_map (r,2), observation_response
(m,2k), covariance (m,m), residual (m,). The real observation response acts
on [Re Delta a(q^2), Im Delta a(q^2)]. For measured cross sections it must
be the derivative of a declared complete bin-integrated KT prediction at
a supplied fit point, not a replacement of cross sections by amplitudes.
The residual is data minus that same reference prediction. All constraints
apply to variations about a reference already satisfying the spectral sum
rule. Results are quadratic local profiles, not nonlinear likelihood fits.
"""

import argparse
import json
import unittest

import numpy as np
from scipy import linalg
import sympy as sp

if __package__:
    from . import derive_vppp_subtraction_matching as matching
else:
    import derive_vppp_subtraction_matching as matching


def derive_zero_width_identity():
    """Derive a delta-spectrum dispersion integral and its existing pole limit."""
    invariant, mass_squared = sp.symbols('q_squared mass_squared', positive=True)
    spectral_variable = sp.Symbol('s_prime', real=True)
    residue = sp.Symbol('spectral_residue', real=True)
    spectral_density = sp.pi * residue * sp.DiracDelta(spectral_variable-mass_squared)
    integrated = sp.integrate(spectral_density / (spectral_variable*(spectral_variable-invariant)),
                             (spectral_variable, 0, sp.oo))*invariant/sp.pi
    _, _, matched_residue, matched_pole, _, _ = matching.derive_single_pole_matching()
    projected = integrated.subs({invariant: matching.projection.total_invariant,
        mass_squared: matching.projection.vector_mass**2, residue: matched_residue})
    return integrated, sp.cancel(projected-matched_pole)


def _rank(matrix, relative_tolerance):
    singular_values = linalg.svdvals(matrix)
    threshold = relative_tolerance * (singular_values[0] if singular_values.size else 0.0)
    return int(np.count_nonzero(singular_values > threshold)), singular_values


def profile_linear_model(coefficient_response, conformal_response, coefficient_constraint,
                         conformal_constraint, covariance, residual, relative_tolerance=1e-10):
    """Enforce sum rule exactly, whiten correlated observations, profile nuisances.

    Constraints read coefficient_constraint @ delta_g + conformal_constraint
    @ delta_c = 0. Infeasible coefficient directions are reported explicitly;
    they are constraints of the selected basis and are never called data rank.
    Columns retain the literal dimensionless coupling coordinates; singular
    values therefore depend on this convention. No theory prior is inserted.
    """
    coefficient_response = np.asarray(coefficient_response, float)
    conformal_response = np.asarray(conformal_response, float)
    coefficient_constraint = np.atleast_2d(coefficient_constraint).astype(float)
    conformal_constraint = np.atleast_2d(conformal_constraint).astype(float)
    covariance, residual = np.asarray(covariance, float), np.asarray(residual, float)
    if coefficient_response.ndim != 2 or coefficient_response.shape[1] != 2:
        raise ValueError('coefficient_response must have two literal-combination columns')
    observations = coefficient_response.shape[0]
    if conformal_response.ndim != 2 or conformal_response.shape[0] != observations:
        raise ValueError('conformal_response observation dimension does not match')
    if coefficient_constraint.shape[1] != 2 or conformal_constraint.shape != (
            coefficient_constraint.shape[0], conformal_response.shape[1]):
        raise ValueError('sum-rule constraint dimensions do not match')
    if covariance.shape != (observations, observations) or residual.shape != (observations,):
        raise ValueError('covariance or residual dimensions do not match observations')
    arrays = [coefficient_response, conformal_response, coefficient_constraint,
              conformal_constraint, covariance, residual]
    if not all(np.isfinite(array).all() for array in arrays):
        raise ValueError('all model and covariance entries must be finite')
    if not np.allclose(covariance, covariance.T, rtol=1e-12, atol=1e-14):
        raise ValueError('covariance must be symmetric')
    cholesky = linalg.cholesky(covariance, lower=True)
    conformal_solution = -linalg.pinv(conformal_constraint, rtol=relative_tolerance) @ coefficient_constraint
    unsatisfied = coefficient_constraint + conformal_constraint @ conformal_solution
    coefficient_allowed = linalg.null_space(unsatisfied, rcond=relative_tolerance) if (
        np.linalg.norm(unsatisfied) > relative_tolerance * max(1., np.linalg.norm(coefficient_constraint))) else np.eye(2)
    conformal_free = linalg.null_space(conformal_constraint, rcond=relative_tolerance)
    response = coefficient_response + conformal_response @ conformal_solution
    whitened_response = linalg.solve_triangular(cholesky, response, lower=True)
    whitened_conformal = linalg.solve_triangular(cholesky, conformal_response @ conformal_free, lower=True)
    whitened_residual = linalg.solve_triangular(cholesky, residual, lower=True)
    # Orthogonal subtraction uses an SVD pseudoinverse; no covariance inverse is formed.
    conformal_inverse = linalg.pinv(whitened_conformal, rtol=relative_tolerance)
    profiled_response = whitened_response-whitened_conformal @ (conformal_inverse @ whitened_response)
    profiled_residual = whitened_residual-whitened_conformal @ (conformal_inverse @ whitened_residual)
    roundoff_threshold = relative_tolerance * max(linalg.norm(whitened_response, 2), 1e-300)
    profiled_response[np.abs(profiled_response) < roundoff_threshold] = 0.0
    allowed_response = profiled_response @ coefficient_allowed
    rank, singular_values = _rank(allowed_response, relative_tolerance)
    coefficients = coefficient_allowed @ (linalg.pinv(allowed_response, rtol=relative_tolerance) @ profiled_residual)
    minimum_residual = profiled_residual-profiled_response @ coefficients
    information = allowed_response.T @ allowed_response
    # Individual curvature exists only if that literal coordinate cannot vary in a data null direction.
    data_null = coefficient_allowed @ linalg.null_space(allowed_response, rcond=relative_tolerance)
    covariance_allowed = linalg.pinv(information, rtol=relative_tolerance)
    coefficient_covariance = coefficient_allowed @ covariance_allowed @ coefficient_allowed.T
    curvatures = []
    for coordinate in range(2):
        if np.linalg.norm(coefficient_allowed[coordinate]) < relative_tolerance:
            curvatures.append(None)  # fixed by exact basis constraint
        elif np.linalg.norm(data_null[coordinate]) > relative_tolerance:
            curvatures.append(0.0)
        else:
            curvatures.append(float(1/coefficient_covariance[coordinate, coordinate]))
    return {'parameter_order': ['g1-g2-g3', '2g4+g5'],
        'interpretation': 'local quadratic profile at caller-supplied reference prediction; no theory priors',
        'allowed_coefficient_basis': coefficient_allowed.tolist(),
        'conformal_particular_map': conformal_solution.tolist(),
        'profiled_rank': rank, 'singular_values': singular_values.tolist(),
        'profiled_information_allowed_coordinates': information.tolist(),
        'minimum_norm_coefficient_shift': coefficients.tolist(),
        'minimum_chi_squared': float(minimum_residual @ minimum_residual),
        'individual_profile_curvatures': curvatures,
        'profile_curvature_definition': 'Delta chi_squared = curvature * (coordinate - minimum)^2; null means fixed by exact basis constraint',
        'unidentified_coefficient_directions': data_null.tolist()}


def calculate_from_spectral_grid(inputs):
    """Build one-subtraction responses and moment constraint from supplied spectra."""
    names = ['spectral_s', 'spectral_weights', 'resonance_spectral', 'conformal_spectral',
             'evaluation_q_squared', 'subtraction_map', 'residue_map', 'observation_response',
             'covariance', 'residual']
    missing = [name for name in names if name not in inputs]
    if missing:
        raise ValueError('Missing explicit physical inputs: '+', '.join(missing))
    spectral_s = np.asarray(inputs['spectral_s'], float)
    weights = np.asarray(inputs['spectral_weights'], float)
    resonance = np.asarray(inputs['resonance_spectral'], float)
    conformal = np.asarray(inputs['conformal_spectral'], float)
    virtualities = np.asarray(inputs['evaluation_q_squared'], complex)
    subtraction = np.asarray(inputs['subtraction_map'], float)
    residues = np.asarray(inputs['residue_map'], float)
    response = np.asarray(inputs['observation_response'], float)
    if spectral_s.ndim != 1 or spectral_s.size < 2 or np.any(spectral_s <= 0) or np.any(np.diff(spectral_s) <= 0):
        raise ValueError('spectral_s must be strictly increasing, positive, with at least two points')
    if weights.shape != spectral_s.shape or np.any(weights < 0):
        raise ValueError('spectral_weights must be nonnegative integration weights on spectral_s')
    if resonance.ndim != 2 or conformal.ndim != 2 or resonance.shape[0] != spectral_s.size or conformal.shape[0] != spectral_s.size:
        raise ValueError('spectral basis arrays must have spectral_s as their first dimension')
    if subtraction.shape != (2,) or residues.shape != (resonance.shape[1], 2):
        raise ValueError('subtraction_map and residue_map require two literal-combination columns')
    if virtualities.ndim != 1 or np.any(virtualities.imag == 0):
        raise ValueError('evaluation points must be explicitly complex off the real cut; no principal-value prescription is inferred')
    if response.ndim != 2 or response.shape[1] != 2*virtualities.size:
        raise ValueError('observation_response must act on concatenated real and imaginary amplitudes')
    if not all(np.isfinite(array).all() for array in [spectral_s, weights, resonance, conformal, virtualities, subtraction, residues, response]):
        raise ValueError('spectral grid inputs must be finite')
    dispersion = virtualities[:, None]*weights[None, :]/(np.pi*spectral_s[None, :]*(spectral_s[None, :]-virtualities[:, None]))
    resonance_basis = dispersion @ resonance
    conformal_basis = dispersion @ conformal
    coefficient_basis = subtraction[None, :] + resonance_basis @ residues
    moment_weights = weights/(np.pi*spectral_s)
    coefficient_constraint = subtraction - (moment_weights @ resonance) @ residues
    conformal_constraint = -(moment_weights @ conformal)
    coefficient_observations = response @ np.vstack([coefficient_basis.real, coefficient_basis.imag])
    conformal_observations = response @ np.vstack([conformal_basis.real, conformal_basis.imag])
    return profile_linear_model(coefficient_observations, conformal_observations,
        coefficient_constraint, conformal_constraint, inputs['covariance'], inputs['residual'])


class ProfileTests(unittest.TestCase):
    """Algebraic fixtures below have no interpretation as physical spectra or data."""

    def test_zero_width_matches_existing_derivation(self):
        _, difference = derive_zero_width_identity()
        self.assertEqual(difference, 0)

    def test_exact_constraint_and_conformal_profiling(self):
        # Sum rule: delta_c0=delta_g0; free c1 duplicates g1's observation column.
        result = profile_linear_model(np.eye(3, 2), np.array([[0., 0.], [0., 1.], [1., 0.]]),
            [[1., 0.]], [[-1., 0.]], np.eye(3), np.zeros(3))
        self.assertEqual(result['profiled_rank'], 1)
        np.testing.assert_allclose(result['individual_profile_curvatures'], [2., 0.], atol=1e-12)
        np.testing.assert_allclose(result['conformal_particular_map'], [[1., 0.], [0., 0.]])

    def test_constraint_without_completion_is_not_data_rank(self):
        result = profile_linear_model(np.eye(2), np.zeros((2, 0)), [[1., 0.]],
            np.zeros((1, 0)), np.eye(2), np.zeros(2))
        self.assertEqual(result['profiled_rank'], 1)
        self.assertIsNone(result['individual_profile_curvatures'][0])

    def test_correlated_covariance_and_minimum(self):
        covariance = np.array([[2., .5], [.5, 1.]])
        result = profile_linear_model(np.eye(2), np.zeros((2, 0)), [[0., 0.]],
            np.zeros((1, 0)), covariance, np.array([2., 3.]))
        np.testing.assert_allclose(result['minimum_norm_coefficient_shift'], [2., 3.])
        np.testing.assert_allclose(result['individual_profile_curvatures'], [.5, 1.])

    def test_spectral_quadrature_delta_limit(self):
        # Discrete delta at s=2 with unit integration weight: Im B=pi delta(s-2).
        # The exact sum rule ties c=2*g0-g1; the combined response is checked
        # against its independently assembled pole expression at complex q^2.
        virtuality = np.array([.5+.2j, 1.+.3j])
        inputs = dict(spectral_s=np.array([1., 2., 3.]), spectral_weights=np.ones(3),
            resonance_spectral=np.array([[0.], [np.pi], [0.]]),
            conformal_spectral=np.array([[0.], [np.pi], [0.]]),
            evaluation_q_squared=virtuality, subtraction_map=np.array([1., 0.]),
            residue_map=np.array([[0., 1.]]), observation_response=np.eye(4),
            covariance=np.eye(4), residual=np.zeros(4))
        result = calculate_from_spectral_grid(inputs)
        expected = 1.+virtuality/(2.-virtuality)
        expected_information = np.sum(np.abs(expected)**2)
        self.assertEqual(result['profiled_rank'], 1)
        np.testing.assert_allclose(result['profiled_information_allowed_coordinates'],
            [[expected_information, 0.], [0., 0.]], atol=1e-12)

    def test_fully_duplicated_nuisance_has_zero_data_rank(self):
        result = profile_linear_model(np.eye(2), np.eye(2), [[0., 0.]],
            [[0., 0.]], np.eye(2), np.zeros(2))
        self.assertEqual(result['profiled_rank'], 0)
        self.assertEqual(result['individual_profile_curvatures'], [0., 0.])

    def test_missing_inputs_fail(self):
        with self.assertRaisesRegex(ValueError, 'Missing explicit physical inputs'):
            calculate_from_spectral_grid({})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', nargs='?', help='NPZ containing explicit spectral basis, matching maps and observation covariance')
    parser.add_argument('--derive-zero-width', action='store_true', help='print the symbolic delta-spectrum limit and matching difference')
    arguments = parser.parse_args()
    if arguments.derive_zero_width:
        integrated, difference = derive_zero_width_identity()
        print('Once-subtracted delta-spectrum integral:\n'+sp.latex(integrated))
        print('Difference from existing matching pole:\n'+sp.latex(difference))
        return
    if arguments.input is None:
        parser.error('explicit physical input NPZ is required; no physical defaults are provided')
    with np.load(arguments.input, allow_pickle=False) as inputs:
        print(json.dumps(calculate_from_spectral_grid(inputs), indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
