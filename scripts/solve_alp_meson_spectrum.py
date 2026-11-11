#!/usr/bin/env python3
"""Poles, pole vectors and residues of a general quadratic theory D(s) = s K - M^2.

Run:   .venv/bin/python scripts/solve_alp_meson_spectrum.py
Tests: .venv/bin/python -m unittest discover -s scripts -p test_solve_alp_meson_spectrum.py

Conventions (Alda:2026nrd, memory/alp/shared/finite-alp-meson-mixing.md and
memory/alp/shared/alp-inverse-propagator-congruence.md): the fields are
Phi_I = (P_i, a), the quadratic Lagrangian is
L_2 = 1/2 dPhi^T K dPhi - 1/2 Phi^T M^2 Phi and the inverse propagator is
D(s) = s K - M^2. K and M^2 are real symmetric constant matrices; K must be
positive definite. Nothing else is assumed about them, so the leading-order
blocks, or blocks including NLO local terms, can be supplied unchanged. A
momentum-dependent D(s) that is not linear in s is outside this module.

LO input. The default LO input is the neutral 6x6 system on
(P_3, P_8, P_0, P_6, P_7, a) built by
    from derive_alp_lo_neutral_blocks import closed_form_blocks, illustrative_lo_inputs, NEUTRAL_FIELD_LABELS
    kinetic, mass_squared = closed_form_blocks(**illustrative_lo_inputs("neutral"))
    spectrum = solve_spectrum_numerical(np.array(kinetic.evalf(30), dtype=float),
                                        np.array(mass_squared.evalf(30), dtype=float),
                                        field_labels=NEUTRAL_FIELD_LABELS)
For sources commuting with the quark charge matrix the charged mesons
(P_1, P_2, P_4, P_5) decouple exactly from this system
(scripts/derive_alp_charged_meson_decoupling.py). The ALP must be the last
field when the default sign convention below is used.

The poles m_n^2 are the zeros of det D(s), i.e. the solutions of the
generalized eigenproblem M^2 v = m^2 K v. The pole vector v^(n) is a null
vector of D(m_n^2), normalized by v^(n)T K v^(n) = 1; vectors belonging to
the same pole are chosen K-orthonormal. With this normalization the residue
of D(s)^{-1} at a simple pole is R_n = v^(n) v^(n)T, and
    sum_n v^(n) v^(n)T = K^{-1},
    D(s)^{-1} = sum_n v^(n) v^(n)T / (s - m_n^2).
The composition of pole n in the basis (P_i, a) is the list of components
v^(n)_I. In the normalization (w^(n), 1) of Alda:2026nrd one has
v^(n) = (w^(n), 1) / sqrt(v^T K v) up to a sign.

Sign convention (a choice, not physics): each v^(n) is multiplied by -1 if
needed so that its component along `sign_reference_index` (default: the last
field, the ALP) is positive; if that component vanishes, the component with
the largest modulus (first one in case of ties) is made positive. For a
degenerate pole the individual vectors are an arbitrary K-orthonormal basis
of the null space; only the sum of their residues is unique.

Two implementations are provided:
  * `solve_spectrum_numerical` (NumPy/SciPy floats): Cholesky reduction
    K = L L^T, symmetric eigenproblem of L^{-1} M^2 L^{-T} with eigenvectors y,
    v = L^{-T} y.
  * `solve_spectrum_exact` (SymPy, exact arithmetic): the zeros of the
    characteristic polynomial det(s K - M^2) = det K det(s - L^{-1} M^2 L^{-T})
    are found with sympy.roots, the null space of D(m_n^2) is computed exactly
    and orthonormalized with the inner product given by K (Gram-Schmidt). This
    is the same eigenproblem as the Cholesky reduction; the Cholesky factor is
    not formed explicitly because it introduces square roots for rational K.
    Positivity of K is checked exactly with the leading principal minors.
    The exact solver requires that sympy.roots finds all poles in closed
    form (e.g. rational poles); otherwise it raises ValueError.

Units: K is dimensionless; M^2, s and m_n^2 are squared masses in a common
unit (GeV^2 in the project); pole vectors have the dimension of K^{-1/2},
i.e. they are dimensionless.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.linalg
import sympy as sp


class KineticMatrixNotPositiveDefiniteError(ValueError):
    """Raised when the kinetic matrix K is not positive definite."""


@dataclass
class PoleSpectrum:
    """Poles m_n^2, pole vectors v^(n) and residues R_n = v^(n) v^(n)T.

    `pole_vectors[:, n]` (NumPy array or SymPy Matrix) is
    v^(n); `residues[n]` is R_n. Poles are sorted in increasing order.
    """

    field_labels: tuple[str, ...]
    kinetic: object
    mass_squared: object
    pole_squared_masses: list
    pole_vectors: object
    residues: list
    exact: bool

    @property
    def size(self) -> int:
        return len(self.field_labels)

    def pole_vector(self, index):
        """Return v^(index) as a column (SymPy Matrix or 1-d NumPy array)."""
        return self.pole_vectors[:, index]

    def composition(self, index) -> dict:
        """Components v^(index)_I of pole `index` in the basis (P_i, a)."""
        vector = self.pole_vector(index)
        return {label: vector[position] for position, label in enumerate(self.field_labels)}

    def completeness_sum(self):
        """sum_n v^(n) v^(n)T, which must equal K^{-1}."""
        if self.exact:
            return sum(self.residues, sp.zeros(self.size))
        return np.sum(self.residues, axis=0)

    def propagator_from_poles(self, squared_momentum):
        """sum_n R_n / (s - m_n^2) at s = `squared_momentum` (not a pole)."""
        if self.exact:
            return sum((residue / (squared_momentum - pole)
                        for pole, residue in zip(self.pole_squared_masses, self.residues)),
                       sp.zeros(self.size))
        return np.sum([residue / (squared_momentum - pole)
                       for pole, residue in zip(self.pole_squared_masses, self.residues)], axis=0)


def inverse_propagator(kinetic, mass_squared, squared_momentum):
    """D(s) = s K - M^2 for SymPy or NumPy matrices."""
    return squared_momentum * kinetic - mass_squared


def _default_labels(size, field_labels):
    if field_labels is None:
        return tuple([f"P_{index}" for index in range(size - 1)] + ["a"])
    if len(field_labels) != size:
        raise ValueError(f"{len(field_labels)} field labels given for a {size}x{size} matrix")
    return tuple(field_labels)


def _check_square_pair(kinetic_shape, mass_shape):
    if kinetic_shape != mass_shape or kinetic_shape[0] != kinetic_shape[1]:
        raise ValueError(f"K and M^2 must be square matrices of equal size, got {kinetic_shape} and {mass_shape}")


# ----------------------------------------------------------------- numerical


def _fix_sign_numerical(vector, sign_reference_index, tolerance):
    reference = vector[sign_reference_index]
    if abs(reference) <= tolerance * np.max(np.abs(vector)):
        reference = vector[int(np.argmax(np.abs(vector)))]
    return -vector if reference < 0 else vector


def solve_spectrum_numerical(kinetic, mass_squared, field_labels=None,
                             sign_reference_index=-1, symmetry_tolerance=1e-12):
    """Solve M^2 v = m^2 K v in floating point by Cholesky reduction.

    Parameters: real symmetric `kinetic` (positive definite) and
    `mass_squared` arrays of equal shape; optional `field_labels`.
    Returns a `PoleSpectrum` with NumPy arrays.
    Raises KineticMatrixNotPositiveDefiniteError if K is not positive definite
    and ValueError if the inputs are not square, equal-sized and symmetric.
    """
    kinetic = np.asarray(kinetic, dtype=float)
    mass_squared = np.asarray(mass_squared, dtype=float)
    _check_square_pair(kinetic.shape, mass_squared.shape)
    for name, matrix in (("K", kinetic), ("M^2", mass_squared)):
        scale = max(1.0, float(np.max(np.abs(matrix))))
        if np.max(np.abs(matrix - matrix.T)) > symmetry_tolerance * scale:
            raise ValueError(f"{name} is not symmetric")
    try:
        cholesky_factor = scipy.linalg.cholesky(kinetic, lower=True)
    except np.linalg.LinAlgError as error:
        smallest = float(np.min(np.linalg.eigvalsh(kinetic)))
        raise KineticMatrixNotPositiveDefiniteError(
            f"K is not positive definite (smallest eigenvalue {smallest:.6g}); "
            "the quadratic theory has a ghost or a non-propagating direction") from error
    half_reduced = scipy.linalg.solve_triangular(cholesky_factor, mass_squared, lower=True)
    reduced = scipy.linalg.solve_triangular(cholesky_factor, half_reduced.T, lower=True)
    reduced = (reduced + reduced.T) / 2
    poles, orthonormal_vectors = np.linalg.eigh(reduced)
    pole_vectors = scipy.linalg.solve_triangular(cholesky_factor.T, orthonormal_vectors, lower=False)
    for column in range(pole_vectors.shape[1]):
        pole_vectors[:, column] = _fix_sign_numerical(pole_vectors[:, column], sign_reference_index, 1e-12)
    residues = [np.outer(pole_vectors[:, column], pole_vectors[:, column]) for column in range(pole_vectors.shape[1])]
    return PoleSpectrum(
        field_labels=_default_labels(kinetic.shape[0], field_labels),
        kinetic=kinetic, mass_squared=mass_squared,
        pole_squared_masses=list(poles), pole_vectors=pole_vectors,
        residues=residues, exact=False)


# --------------------------------------------------------------------- exact


def exact_positive_definite(kinetic) -> bool:
    """Sylvester criterion: all leading principal minors of K are positive."""
    return all(sp.simplify(kinetic[:order, :order].det()) > 0 for order in range(1, kinetic.rows + 1))


def _fix_sign_exact(vector, sign_reference_index):
    reference = vector[sign_reference_index]
    if reference == 0:
        magnitudes = [abs(sp.N(component)) for component in vector]
        reference = vector[magnitudes.index(max(magnitudes))]
    return -vector if sp.N(reference) < 0 else vector


def _kinetic_orthonormal_basis(vectors, kinetic):
    basis = []
    for vector in vectors:
        for previous in basis:
            vector = vector - (previous.T * kinetic * vector)[0] * previous
        norm_squared = sp.simplify((vector.T * kinetic * vector)[0])
        basis.append((vector / sp.sqrt(norm_squared)).applyfunc(sp.radsimp))
    return basis


def solve_spectrum_exact(kinetic, mass_squared, field_labels=None, sign_reference_index=-1):
    """Solve M^2 v = m^2 K v in exact SymPy arithmetic.

    Parameters: real symmetric SymPy-compatible `kinetic` (positive definite)
    and `mass_squared` of equal shape. Returns a `PoleSpectrum` with SymPy
    matrices. Raises KineticMatrixNotPositiveDefiniteError if K is not
    positive definite and ValueError if not all poles are found in closed form.
    """
    kinetic = sp.Matrix(kinetic)
    mass_squared = sp.Matrix(mass_squared)
    _check_square_pair(kinetic.shape, mass_squared.shape)
    for name, matrix in (("K", kinetic), ("M^2", mass_squared)):
        if (matrix - matrix.T).applyfunc(sp.simplify) != sp.zeros(*matrix.shape):
            raise ValueError(f"{name} is not symmetric")
    if not exact_positive_definite(kinetic):
        minors = [sp.simplify(kinetic[:order, :order].det()) for order in range(1, kinetic.rows + 1)]
        raise KineticMatrixNotPositiveDefiniteError(
            f"K is not positive definite (leading principal minors {minors}); "
            "the quadratic theory has a ghost or a non-propagating direction")
    squared_momentum = sp.Symbol("s")
    characteristic = sp.Poly(sp.expand(
        inverse_propagator(kinetic, mass_squared, squared_momentum).det(method="berkowitz")), squared_momentum)
    pole_multiplicities = sp.roots(characteristic)
    if sum(pole_multiplicities.values()) != kinetic.rows:
        raise ValueError("not all zeros of det(s K - M^2) were found in closed form; "
                         "use solve_spectrum_numerical")
    poles, vectors = [], []
    for pole in sorted(pole_multiplicities, key=lambda value: sp.N(value)):
        null_space = inverse_propagator(kinetic, mass_squared, pole).nullspace(simplify=True)
        if len(null_space) != pole_multiplicities[pole]:
            raise ValueError(f"null space of D({pole}) has dimension {len(null_space)}, "
                             f"expected multiplicity {pole_multiplicities[pole]}")
        for vector in _kinetic_orthonormal_basis(null_space, kinetic):
            poles.append(pole)
            vectors.append(_fix_sign_exact(vector, sign_reference_index))
    pole_vectors = sp.Matrix.hstack(*vectors)
    residues = [(vector * vector.T).applyfunc(sp.radsimp) for vector in vectors]
    return PoleSpectrum(
        field_labels=_default_labels(kinetic.rows, field_labels),
        kinetic=kinetic, mass_squared=mass_squared,
        pole_squared_masses=poles, pole_vectors=pole_vectors,
        residues=residues, exact=True)


# ------------------------------------------------------------ demonstration


def main():
    """Solve a 4x4 exact example with rational poles, then the LO neutral 6x6 system."""
    from derive_alp_lo_neutral_blocks import NEUTRAL_FIELD_LABELS
    from verify_alp_lo_spectrum import lo_blocks_numerical, random_rational_quadratic_form

    data = random_rational_quadratic_form(size=4, seed=20261111)
    exact = solve_spectrum_exact(data["kinetic"], data["mass_squared"])
    numerical = solve_spectrum_numerical(np.array(data["kinetic"], dtype=float),
                                         np.array(data["mass_squared"], dtype=float))
    print("Exact rational example on (P_0, P_1, P_2, a); illustrative numbers.")
    print(f"poles (exact)     : {exact.pole_squared_masses}")
    print(f"poles (numerical) : {[float(pole) for pole in numerical.pole_squared_masses]}")
    completeness = (exact.completeness_sum() - exact.kinetic.inv()).applyfunc(sp.simplify)
    print(f"sum_n v v^T = K^-1 exactly: {completeness == sp.zeros(4)}")
    for index in range(4):
        print(f"pole {index}: m^2 = {exact.pole_squared_masses[index]}, composition = {exact.composition(index)}")
    print()
    lo = solve_spectrum_numerical(*lo_blocks_numerical("neutral"), field_labels=NEUTRAL_FIELD_LABELS)
    print(f"LO neutral system {NEUTRAL_FIELD_LABELS}, generic neutral sources; illustrative inputs (GeV^2).")
    for index in range(lo.size):
        composition = {label: round(float(value), 6) for label, value in lo.composition(index).items()}
        print(f"pole {index}: m^2 = {lo.pole_squared_masses[index]:.9g}, composition = {composition}")


if __name__ == "__main__":
    main()
