#!/usr/bin/env python3
"""LO pseudoscalar and ALP masses for finite mixing in the invariant frame.

The coefficients are derivatives at a/f_a = 0, not Taylor coefficients.
All masses, f_a, f_pi, B_0 and the scalar and pseudoscalar source matrices
are in GeV; m_0_squared is in GeV^2. The sources are 3 by 3 matrices in
(u,d,s) order. The vacuum is assumed to have S_0 diagonal, P_0 = 0 and
Theta_0 = 0, with any ALP tadpole cancelled. P_1 must conserve electric
charge. S_0, P_1 and S_2 must be Hermitian. Only S_0, S_2, P_1 and Theta_1
enter the quadratic Lagrangian; other list elements are not used.

The six neutral poles are eigenvalues of the canonical mass matrix. Labels
are assigned in one-to-one fashion by maximizing the sum of squared overlaps
with the unmixed meson and ALP fields. K0_1 is the neutral-kaon combination
coupled directly to the ALP, and K0_2 is the orthogonal combination. For zero
coupling, they are the two real kaon fields. Labels at degeneracies and near
maximal mixing are conventional; compare the full set of poles there.

Example (illustrative inputs, not a fit):
    PYTHONPATH=scripts .venv/bin/python -c 'import numpy as np; from calculate_alp_meson_masses import calculate_masses; z=np.zeros((3,3)); print(calculate_masses(1., .6, [np.diag([.003,.005,.1]),z,z], [z,z], [0.,0.], f_pi=.092, B_0=2.5, m_0_squared=.75))'
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import linear_sum_assignment


def _matrix(value, name):
    result = np.asarray(value, dtype=complex)
    if result.shape != (3, 3) or not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must be a finite 3 by 3 matrix")
    if not np.allclose(result, result.conj().T, rtol=1e-12, atol=1e-14):
        raise ValueError(f"{name} must be Hermitian")
    return result


def calculate_masses(f_a, m_a0, hat_S, hat_P, hat_Theta, *, f_pi, B_0, m_0_squared):
    """Return a dict of LO masses in GeV, including finite ALP-meson mixing.

    The positional source arguments are lists indexed by derivative order.
    Required entries are hat_S[0], hat_S[2], hat_P[0], hat_P[1],
    hat_Theta[0] and hat_Theta[1]. Higher derivatives do not affect the LO
    mass matrix. The keyword inputs are the independent LO chiral constants.
    Neutral kaon keys denote real mass eigenstates, while charged-particle and
    antiparticle masses have separate keys. A negative pole raises ValueError.
    """
    if len(hat_S) < 3 or len(hat_P) < 2 or len(hat_Theta) < 2:
        raise ValueError("source lists must contain S[0:3], P[0:2] and Theta[0:2]")
    if not all(np.isfinite(x) for x in (f_a, m_a0, f_pi, B_0, m_0_squared)):
        raise ValueError("mass inputs must be finite real numbers")
    if f_a <= 0 or f_pi <= 0 or m_a0 < 0:
        raise ValueError("f_a and f_pi must be positive and m_a0 nonnegative")
    scalar_zero = _matrix(hat_S[0], "hat_S[0]")
    scalar_second = _matrix(hat_S[2], "hat_S[2]")
    pseudoscalar_zero = _matrix(hat_P[0], "hat_P[0]")
    pseudoscalar_first = _matrix(hat_P[1], "hat_P[1]")
    if not np.allclose(scalar_zero, np.diag(np.diag(scalar_zero)), rtol=0, atol=1e-14):
        raise ValueError("hat_S[0] must be diagonal in the quark mass basis")
    if not np.allclose(pseudoscalar_zero, 0, rtol=0, atol=1e-14) or abs(hat_Theta[0]) > 1e-14:
        raise ValueError("the aligned vacuum requires hat_P[0] = 0 and hat_Theta[0] = 0")
    if not np.isreal(hat_Theta[1]) or not np.isfinite(hat_Theta[1]):
        raise ValueError("hat_Theta[1] must be finite and real")
    if not np.allclose(pseudoscalar_first[0, 1:], 0, rtol=0, atol=1e-14):
        raise ValueError("hat_P[1] must conserve electric charge")

    root3 = np.sqrt(3.)
    generators = np.array([
        np.diag([1, -1, 0]),
        np.diag([1, 1, -2]) / root3,
        np.eye(3) * np.sqrt(2. / 3.),
        [[0, 0, 0], [0, 0, 1], [0, 1, 0]],
        [[0, 0, 0], [0, 0, -1j], [0, 1j, 0]],
    ], dtype=complex)
    chi_zero = 2 * B_0 * scalar_zero
    ratio = f_pi / f_a
    mass_squared = np.zeros((6, 6))
    for i, first in enumerate(generators):
        for j, second in enumerate(generators):
            mass_squared[i, j] = np.real(
                np.trace((first @ second + second @ first) @ chi_zero) / 4
                + m_0_squared * np.trace(first) * np.trace(second) / 6)
        mass_squared[i, 5] = mass_squared[5, i] = np.real(
            -ratio * B_0 * np.trace(first @ pseudoscalar_first)
            + ratio * m_0_squared * hat_Theta[1] * np.trace(first) / 3)
    mass_squared[5, 5] = (m_a0**2 + ratio**2 *
                          (2 * m_0_squared * hat_Theta[1]**2 / 3
                           - B_0 * np.trace(scalar_second).real))

    # Construct a reference basis from unmixed pi0, eta, eta-prime and kaons.
    _, diagonal_vectors = np.linalg.eigh(mass_squared[:3, :3])
    reference = np.zeros((6, 6))
    reference[:3, :3] = diagonal_vectors
    kaon_couplings = mass_squared[3:5, 5]
    kaon_norm = np.linalg.norm(kaon_couplings)
    if kaon_norm:
        reference[3:5, 3] = kaon_couplings / kaon_norm
        reference[3:5, 4] = [-reference[4, 3], reference[3, 3]]
    else:
        reference[3:5, 3:5] = np.eye(2)
    reference[5, 5] = 1
    poles, eigenvectors = np.linalg.eigh(mass_squared)
    # Permit roundoff at a zero eigenvalue; a genuinely negative pole is unstable.
    tolerance = 1e-12 * max(1., np.max(np.abs(poles)))
    charged_squared = B_0 * np.real(np.diag(scalar_zero)[0] + np.diag(scalar_zero)[1:])
    if np.min(poles) < -tolerance or np.min(charged_squared) < -tolerance:
        raise ValueError("the LO mass matrix has a negative squared mass")
    poles = np.maximum(poles, 0.)
    charged_squared = np.maximum(charged_squared, 0.)
    rows, columns = linear_sum_assignment(-np.abs(reference.T @ eigenvectors)**2)
    labels = ("pi0", "eta", "eta_prime", "K0_1", "K0_2", "a")
    result = {labels[row]: float(np.sqrt(poles[column])) for row, column in zip(rows, columns)}
    result.update({"pi+": float(np.sqrt(charged_squared[0])),
                   "pi-": float(np.sqrt(charged_squared[0])),
                   "K+": float(np.sqrt(charged_squared[1])),
                   "K-": float(np.sqrt(charged_squared[1]))})
    return result
