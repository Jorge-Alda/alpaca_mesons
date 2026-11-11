#!/usr/bin/env python3
"""Stability and near-degeneracy diagnostics of D(s) = s K - M^2 (tasks 2.2 and 2.3).

Run:   .venv/bin/python scripts/diagnose_alp_meson_spectrum.py
Tests: .venv/bin/python -m unittest discover -s scripts -p test_diagnose_alp_meson_spectrum.py

Conventions as in scripts/solve_alp_meson_spectrum.py (Alda:2026nrd;
memory/alp/shared/finite-alp-meson-mixing.md): Phi = (P_i, a) with the ALP
last, L_2 = 1/2 dPhi^T K dPhi - 1/2 Phi^T M^2 Phi, poles m_n^2 with
M^2 v = m^2 K v, v^(n)T K v^(n) = 1, residues R_n = v^(n) v^(n)T. The LO
input is the neutral 6x6 system (P_3, P_8, P_0, P_6, P_7, a) of
derive_alp_lo_neutral_blocks.py.

Stability (task 2.2). Conditions for a physical spectrum:
  (i)   K positive definite (no ghost or non-propagating field);
  (ii)  M^2 positive semidefinite (stationary point is a minimum: vacuum stability);
  (iii) all poles real and nonnegative.
For K > 0 the poles are real and, by Sylvester's law of inertia applied to
L^{-1} M^2 L^{-T} (K = L L^T), the number of negative poles equals the number
of negative eigenvalues of M^2; so (ii) <=> (iii) once (i) holds. All three
are evaluated independently as a consistency check. `solve_physical_spectrum`
raises UnstableVacuumError, naming every violated condition, instead of
returning a spectrum. At LO, K_PP = 1, K_aa = 1 + K_ia K_ia, and
  det K = K_aa - K_ia K_ia = 1   (Schur complement of the unit block),
  det(K - t) = (1 - t)^(N-1) (t^2 - (2 + kappa) t + 1),  kappa = K_ia K_ia,
so the two nontrivial eigenvalues have product 1 and positive sum: K > 0 for
every real K_ia (checked symbolically for N = 5 and for symbolic sources).
Stability therefore reduces at LO to M^2 >= 0. With M^2_PP > 0 this is the
Schur condition M^2_aa - M^2_aP (M^2_PP)^{-1} M^2_Pa >= 0, i.e.
  m_a0^2 >= M^2_aP (M^2_PP)^{-1} M^2_Pa - 2 f_pi^2 m_0^2 c_G^2/(N_f f_a^2)
(lower bound on the bare ALP mass; computed exactly for the illustrative inputs).

Degeneracy (task 2.3). Sorted poles are chained into groups when neighbouring
poles differ by at most `relative_tolerance` times max(|m_n^2|) (a stated,
adjustable tolerance; default 1e-2). For each group the output contains:
  * the spread and the minimal gap to the rest of the spectrum;
  * a K-orthonormal basis of the group (for an exactly degenerate group, with
    gap below `exact_relative_tolerance`, the basis is rotated so that at most
    one vector has a nonzero ALP component; the individual vectors are not
    unique, only the sum of residues R_G = sum_{n in G} v v^T and the
    projector P_G = R_G K are);
  * the ALP weight (R_G)_aa = sum_{n in G} (v^(n)_a)^2 (basis independent) and
    the number of states carrying an ALP component above
    `alp_component_tolerance`;
  * whether the group contains the ALP-like pole (largest |v_a| overall);
  * for two resolved ALP-carrying states, the mixing angle
    theta = arctan(|v_a^lower| / |v_a^upper|) in [0, pi/2] (theta = 0: lower
    state without ALP component; pi/4: maximal mixing; pi/2: lower state
    carries all the ALP component of the pair);
  * if K_PP is the unit matrix (LO), the two-level estimate of
    finite-alp-meson-mixing.md with canonical blocks (P = P' + v a, v = -K_Pa):
      M^2_{P'a,can} = u.(M^2_Pa + M^2_PP v),
      M^2_{aa,can} = M^2_aa + 2 v.M^2_Pa + v.M^2_PP v,
      m_pm^2 = (m_P^2 + M^2_aa,can)/2 +- sqrt((M^2_aa,can - m_P^2)^2 + 4 (M^2_{P'a,can})^2)/2,
    for the eigenvector u of M^2_PP with eigenvalue m_P^2 closest to the group
    (a degenerate meson eigenspace is rotated so that only one combination
    couples to the ALP). At M^2_aa,can = m_P^2 the splitting is 2|M^2_{P'a,can}|,
    and theta = atan2(2|M^2_{P'a,can}|, M^2_aa,can - m_P^2)/2 in [0, pi/2]
    (same convention as above; exact for an isolated two-level system).
A group is flagged as an ALP-meson near-degeneracy (where the small-mixing
expansion of the ALP pole fails) only if it contains the ALP-like pole and at
least two states with a nonzero ALP component. Groups whose members other
than one have zero ALP component (e.g. the K0 doublet at m_K0^2 for
flavour-diagonal sources, or the exact m_K0^2 pole near the ALP-like pole)
are reported as degeneracies without that flag.

The single-pole formula w^(a)(s) = (M^2_PP - s K_PP)^{-1}(s K_Pa - M^2_Pa)
fails where M^2_PP - m_n^2 K_PP is singular (exact m_K0^2 pole, ALP-like
pole exactly degenerate with a decoupled meson); `single_pole_formula_singular`
reports this for every pole.

Units: K dimensionless; M^2, s, m_n^2, gaps and splittings in GeV^2 for the
LO blocks; angles in radians. Inputs are illustrative, not fits.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass, field

import numpy as np
import scipy.linalg
import sympy as sp

from derive_alp_charged_meson_decoupling import symbolic_hermitian
from derive_alp_lo_neutral_blocks import NEUTRAL_FIELD_LABELS, closed_form_blocks, illustrative_lo_inputs
from derive_alp_neutral_block_reduction import diagonal_quark_spurion
from solve_alp_meson_spectrum import (
    PoleSpectrum,
    solve_spectrum_exact,
    solve_spectrum_numerical,
)

KINETIC_CONDITION = "K is not positive definite"
MASS_CONDITION = "M^2 has a negative eigenvalue (unstable vacuum)"
POLE_CONDITION = "not all poles are real and nonnegative"


# ------------------------------------------------------------- stability


class UnstableVacuumError(ValueError):
    """Raised instead of returning a spectrum when a stability condition fails.

    `violated_conditions` lists the failed conditions by name and
    `diagnostics` holds the full StabilityDiagnostics.
    """

    def __init__(self, diagnostics):
        self.diagnostics = diagnostics
        self.violated_conditions = list(diagnostics.violated_conditions)
        super().__init__("no physical spectrum: " + "; ".join(diagnostics.describe()))


@dataclass
class StabilityDiagnostics:
    """Positivity of K, of M^2 and of the poles (floats, or SymPy numbers if exact)."""

    kinetic_smallest_eigenvalue: object
    kinetic_determinant: object
    kinetic_positive_definite: bool
    mass_squared_smallest_eigenvalue: object
    mass_squared_positive_semidefinite: bool
    poles: list
    poles_real_and_nonnegative: bool
    exact: bool = False
    violated_conditions: list = field(default_factory=list)

    @property
    def stable(self) -> bool:
        return not self.violated_conditions

    def describe(self):
        """One line per violated condition, with the number that violates it."""
        lines = []
        for condition in self.violated_conditions:
            if condition == KINETIC_CONDITION:
                lines.append(f"{condition} (smallest eigenvalue {self.kinetic_smallest_eigenvalue}, "
                             f"det K = {self.kinetic_determinant})")
            elif condition == MASS_CONDITION:
                lines.append(f"{condition} (smallest eigenvalue {self.mass_squared_smallest_eigenvalue})")
            else:
                lines.append(f"{condition} (poles {self.poles})")
        return lines


def _collect_violations(diagnostics):
    diagnostics.violated_conditions = [
        condition for condition, satisfied in (
            (KINETIC_CONDITION, diagnostics.kinetic_positive_definite),
            (MASS_CONDITION, diagnostics.mass_squared_positive_semidefinite),
            (POLE_CONDITION, diagnostics.poles_real_and_nonnegative))
        if not satisfied]
    return diagnostics


def stability_diagnostics(kinetic, mass_squared, relative_tolerance=1e-12):
    """Floating-point stability diagnostics of (K, M^2); never raises for unstable input.

    A tolerance relative to max|M^2| (or max|K|) is allowed for zero
    eigenvalues and zero poles (e.g. a massless state in the chiral limit).
    """
    kinetic = np.asarray(kinetic, dtype=float)
    mass_squared = np.asarray(mass_squared, dtype=float)
    kinetic_eigenvalues = np.linalg.eigvalsh(kinetic)
    mass_eigenvalues = np.linalg.eigvalsh(mass_squared)
    kinetic_scale = max(1.0, float(np.max(np.abs(kinetic))))
    mass_scale = max(1.0, float(np.max(np.abs(mass_squared))))
    kinetic_positive = bool(kinetic_eigenvalues[0] > relative_tolerance * kinetic_scale)
    poles = scipy.linalg.eigvals(mass_squared, kinetic)
    finite = np.all(np.isfinite(poles))
    real = finite and bool(np.max(np.abs(poles.imag)) <= 1e3 * relative_tolerance * mass_scale)
    nonnegative = real and bool(np.min(poles.real) >= -1e3 * relative_tolerance * mass_scale)
    return _collect_violations(StabilityDiagnostics(
        kinetic_smallest_eigenvalue=float(kinetic_eigenvalues[0]),
        kinetic_determinant=float(np.linalg.det(kinetic)),
        kinetic_positive_definite=kinetic_positive,
        mass_squared_smallest_eigenvalue=float(mass_eigenvalues[0]),
        mass_squared_positive_semidefinite=bool(mass_eigenvalues[0] >= -relative_tolerance * mass_scale),
        poles=sorted(poles.real) if real else list(poles),
        poles_real_and_nonnegative=bool(kinetic_positive and nonnegative),
    ))


def _principal_minors(matrix):
    size = matrix.rows
    for order in range(1, size + 1):
        for rows in itertools.combinations(range(size), order):
            yield matrix.extract(list(rows), list(rows)).det()


def exact_stability_diagnostics(kinetic, mass_squared):
    """Exact (SymPy) diagnostics: Sylvester criteria and the real roots of det(s K - M^2).

    K > 0: all leading principal minors positive. M^2 >= 0: all principal
    minors nonnegative. Poles: real roots of the characteristic polynomial
    (with multiplicity) must number dim K and be >= 0.
    """
    kinetic, mass_squared = sp.Matrix(kinetic), sp.Matrix(mass_squared)
    size = kinetic.rows
    leading_minors = [sp.simplify(kinetic[:order, :order].det()) for order in range(1, size + 1)]
    kinetic_positive = all(minor > 0 for minor in leading_minors)
    mass_positive = all(sp.simplify(minor) >= 0 for minor in _principal_minors(mass_squared))
    squared_momentum = sp.Symbol("s")
    polynomial = sp.Poly(sp.expand((squared_momentum * kinetic - mass_squared).det(method="berkowitz")), squared_momentum)
    real_roots = polynomial.real_roots() if polynomial.degree() > 0 else []
    poles_ok = polynomial.degree() == size and len(real_roots) == size and all(root >= 0 for root in real_roots)
    kinetic_eigenvalues = sp.Poly(kinetic.charpoly()).real_roots()
    mass_eigenvalues = sp.Poly(mass_squared.charpoly()).real_roots()
    return _collect_violations(StabilityDiagnostics(
        kinetic_smallest_eigenvalue=min(kinetic_eigenvalues) if kinetic_eigenvalues else None,
        kinetic_determinant=leading_minors[-1],
        kinetic_positive_definite=bool(kinetic_positive),
        mass_squared_smallest_eigenvalue=min(mass_eigenvalues),
        mass_squared_positive_semidefinite=bool(mass_positive),
        poles=list(real_roots),
        poles_real_and_nonnegative=bool(kinetic_positive and poles_ok),
        exact=True,
    ))


def solve_physical_spectrum(kinetic, mass_squared, field_labels=None, relative_tolerance=1e-12):
    """Solve in floating point only if (K, M^2) passes the stability conditions.

    Returns (PoleSpectrum, StabilityDiagnostics). Raises UnstableVacuumError
    naming the violated conditions otherwise; the spectrum is then not returned.
    """
    diagnostics = stability_diagnostics(kinetic, mass_squared, relative_tolerance)
    if not diagnostics.stable:
        raise UnstableVacuumError(diagnostics)
    return solve_spectrum_numerical(kinetic, mass_squared, field_labels=field_labels), diagnostics


def solve_physical_spectrum_exact(kinetic, mass_squared, field_labels=None):
    """Exact counterpart of solve_physical_spectrum (SymPy)."""
    diagnostics = exact_stability_diagnostics(kinetic, mass_squared)
    if not diagnostics.stable:
        raise UnstableVacuumError(diagnostics)
    return solve_spectrum_exact(kinetic, mass_squared, field_labels=field_labels), diagnostics


# ----------------------------------------------- LO kinetic positivity proof


def lo_kinetic_matrix_checks(meson_count=5):
    """Symbolic: det K = 1, leading minors 1 and the characteristic polynomial for the LO K."""
    mixing = sp.Matrix(sp.symbols(f"K_1:{meson_count + 1}a", real=True))
    size = meson_count + 1
    kinetic = sp.eye(size)
    kinetic[:meson_count, meson_count] = mixing
    kinetic[meson_count, :meson_count] = mixing.T
    kappa = (mixing.T * mixing)[0]
    kinetic[meson_count, meson_count] = 1 + kappa
    variable = sp.Symbol("t")
    characteristic = sp.expand((kinetic - variable * sp.eye(size)).det(method="berkowitz"))
    expected = sp.expand((1 - variable)**(size - 2) * (variable**2 - (2 + kappa) * variable + 1))
    smaller, larger = sp.solve(variable**2 - (2 + kappa) * variable + 1, variable)
    symbolic_inputs = {
        "axial_source": symbolic_hermitian("cA"), "pseudoscalar_source": symbolic_hermitian("cP"),
        "quark_spurion": diagonal_quark_spurion(), "gluon_coupling": sp.Symbol("c_G", real=True),
        "pion_decay_constant": sp.Symbol("f_pi", positive=True), "alp_decay_constant": sp.Symbol("f_a", positive=True),
        "anomaly_mass_squared": sp.Symbol("m_0_sq", real=True), "bare_alp_mass_squared": sp.Symbol("m_a0_sq", real=True),
    }
    lo_kinetic, _ = closed_form_blocks(**symbolic_inputs)
    return {"checks": {
        "det K = 1 for symbolic K_ia": sp.expand(kinetic.det(method="berkowitz")) == 1,
        "all leading principal minors of K equal 1": all(
            sp.expand(kinetic[:order, :order].det()) == 1 for order in range(1, size + 1)),
        "det(K - t) = (1 - t)^(N-1) (t^2 - (2 + kappa) t + 1)": sp.expand(characteristic - expected) == 0,
        "nontrivial eigenvalues: product 1, sum 2 + kappa > 0, so both positive": (
            sp.simplify(smaller * larger - 1) == 0 and sp.simplify(smaller + larger - 2 - kappa) == 0),
        "det K = 1 for the LO blocks with symbolic charge-conserving sources":
            sp.simplify(lo_kinetic.det(method="berkowitz")) == 1,
    }, "smallest_eigenvalue": smaller}


def lo_bare_mass_stability_bound(inputs):
    """Exact lower bound on m_a0^2 for M^2 >= 0 at LO (requires M^2_PP > 0); returns a dict."""
    kinetic, mass_squared = closed_form_blocks(**inputs)
    meson = mass_squared[:-1, :-1]
    mixing = mass_squared[:-1, -1]
    induced = sp.simplify(mass_squared[-1, -1] - inputs["bare_alp_mass_squared"])
    bound = sp.nsimplify(sp.simplify((mixing.T * meson.inv() * mixing)[0] - induced))
    meson_positive = all(sp.simplify(meson[:order, :order].det()) > 0 for order in range(1, meson.rows + 1))
    return {"bound": bound, "meson_block_positive_definite": meson_positive, "kinetic": kinetic}


# ------------------------------------------------------------- degeneracy


@dataclass
class NearDegenerateGroup:
    """A group of (near-)degenerate poles and its basis-independent content."""

    pole_indices: list
    poles: list
    spread: float
    gap_to_rest: float
    exactly_degenerate: bool
    basis: object
    residue_sum: object
    projector: object
    alp_weight: float
    alp_components: list
    alp_carrying_state_count: int
    contains_alp_like_pole: bool
    mixing_angle: object = None
    two_level: object = None

    @property
    def alp_meson_near_degeneracy(self) -> bool:
        """ALP-like pole near-degenerate with at least one other ALP-carrying state."""
        return self.contains_alp_like_pole and self.alp_carrying_state_count >= 2


@dataclass
class DegeneracyReport:
    """All groups of near-degenerate poles and the minimal gap of the spectrum."""

    relative_tolerance: float
    minimal_gap: float
    minimal_gap_indices: tuple
    groups: list

    @property
    def alp_meson_near_degeneracies(self):
        return [group for group in self.groups if group.alp_meson_near_degeneracy]

    @property
    def other_degeneracies(self):
        """Degeneracies not involving ALP mixing of the ALP-like pole (e.g. the K0 doublet)."""
        return [group for group in self.groups if not group.alp_meson_near_degeneracy]


def _chain_groups(values, tolerance):
    groups = [[0]]
    for index in range(1, len(values)):
        if abs(values[index] - values[groups[-1][-1]]) <= tolerance:
            groups[-1].append(index)
        else:
            groups.append([index])
    return groups


def alp_adapted_basis(basis, alp_index=-1):
    """Rotate the columns of `basis` (orthogonally) so that only the first has an ALP component.

    The K-orthonormality of the columns is preserved. Used for exactly
    degenerate groups, whose basis is otherwise arbitrary.
    """
    row = basis[alp_index, :]
    norm = float(np.linalg.norm(row))
    if basis.shape[1] == 1 or norm == 0.0:
        return basis
    # Householder reflection H with H e_1 = row/|row|; H is symmetric orthogonal.
    target = row / norm
    difference = np.zeros_like(target)
    difference[0] = 1.0
    difference = difference - target
    if np.linalg.norm(difference) < 1e-15:
        rotation = np.eye(basis.shape[1])
    else:
        difference = difference / np.linalg.norm(difference)
        rotation = np.eye(basis.shape[1]) - 2 * np.outer(difference, difference)
    adapted = basis @ rotation
    if adapted[alp_index, 0] < 0:
        adapted[:, 0] = -adapted[:, 0]
    return adapted


def two_level_estimates(kinetic, mass_squared, alp_index=-1, relative_tolerance=1e-9):
    """Two-level formulas of finite-alp-meson-mixing.md for every meson eigenstate (needs K_PP = 1).

    Returns a list of dicts (one per eigenvalue of M^2_PP; a degenerate
    eigenspace is rotated so that only one combination couples to the ALP and
    the others appear with zero coupling), or None if K_PP is not the unit matrix.
    """
    size = kinetic.shape[0]
    alp = alp_index % size
    meson = [index for index in range(size) if index != alp]
    kinetic_meson = kinetic[np.ix_(meson, meson)]
    if np.max(np.abs(kinetic_meson - np.eye(len(meson)))) > 1e-12:
        return None
    meson_mass = mass_squared[np.ix_(meson, meson)]
    shift = -kinetic[meson, alp]
    canonical_mixing = mass_squared[meson, alp] + meson_mass @ shift
    canonical_alp = mass_squared[alp, alp] + 2 * shift @ mass_squared[meson, alp] + shift @ meson_mass @ shift
    values, vectors = np.linalg.eigh(meson_mass)
    scale = max(1.0, float(np.max(np.abs(values))))
    estimates = []
    for group in _chain_groups(values, relative_tolerance * scale):
        mass = float(np.mean(values[group]))
        couplings = vectors[:, group].T @ canonical_mixing
        coupling = float(np.linalg.norm(couplings))
        for position in range(len(group)):
            value = coupling if position == 0 else 0.0
            difference = canonical_alp - mass
            root = np.sqrt(difference**2 + 4 * value**2)
            estimates.append({
                "meson_mass_squared": mass,
                "canonical_mixing": value,
                "canonical_alp_mass_squared": float(canonical_alp),
                "poles": ((mass + canonical_alp - root) / 2, (mass + canonical_alp + root) / 2),
                "splitting": float(root),
                "splitting_at_degeneracy": 2 * abs(value),
                "mixing_angle": 0.5 * float(np.arctan2(2 * abs(value), difference)),
            })
    return estimates


def find_near_degeneracies(spectrum: PoleSpectrum, relative_tolerance=1e-2, exact_relative_tolerance=1e-9,
                           alp_component_tolerance=1e-8, alp_index=-1):
    """Groups of near-degenerate poles of a numerical PoleSpectrum (see module docstring)."""
    poles = np.array(spectrum.pole_squared_masses, dtype=float)
    vectors = np.array(spectrum.pole_vectors, dtype=float)
    kinetic = np.array(spectrum.kinetic, dtype=float)
    mass_squared = np.array(spectrum.mass_squared, dtype=float)
    scale = max(1e-300, float(np.max(np.abs(poles))))
    gaps = np.diff(poles)
    minimal = int(np.argmin(gaps)) if len(gaps) else 0
    alp_like = int(np.argmax(np.abs(vectors[alp_index, :])))
    estimates = two_level_estimates(kinetic, mass_squared, alp_index)
    groups = []
    for members in _chain_groups(poles, relative_tolerance * scale):
        if len(members) < 2:
            continue
        spread = float(poles[members[-1]] - poles[members[0]])
        outside = [abs(poles[index] - poles[member]) for member in members
                   for index in range(len(poles)) if index not in members]
        exactly = spread <= exact_relative_tolerance * scale
        basis = vectors[:, members]
        if exactly:
            basis = alp_adapted_basis(basis, alp_index)
        residue_sum = basis @ basis.T
        alp_components = [float(value) for value in basis[alp_index, :]]
        carrying = [position for position, value in enumerate(alp_components) if abs(value) > alp_component_tolerance]
        angle = None
        if not exactly and len(carrying) == 2:
            lower, upper = carrying
            angle = float(np.arctan2(abs(alp_components[lower]), abs(alp_components[upper])))
        two_level = None
        if estimates is not None and len(carrying) >= 2:
            centre = float(np.mean(poles[members]))
            coupled = [entry for entry in estimates if entry["canonical_mixing"] != 0.0] or estimates
            two_level = min(coupled, key=lambda entry: abs(entry["meson_mass_squared"] - centre))
        groups.append(NearDegenerateGroup(
            pole_indices=list(members), poles=[float(poles[index]) for index in members],
            spread=spread, gap_to_rest=float(min(outside)) if outside else float("inf"),
            exactly_degenerate=bool(exactly), basis=basis, residue_sum=residue_sum,
            projector=residue_sum @ kinetic, alp_weight=float(residue_sum[alp_index, alp_index]),
            alp_components=alp_components, alp_carrying_state_count=len(carrying),
            contains_alp_like_pole=alp_like in members, mixing_angle=angle, two_level=two_level))
    return DegeneracyReport(
        relative_tolerance=relative_tolerance,
        minimal_gap=float(gaps[minimal]) if len(gaps) else float("inf"),
        minimal_gap_indices=(minimal, minimal + 1) if len(gaps) else (),
        groups=groups)


def exact_degenerate_groups(spectrum: PoleSpectrum, alp_index=-1):
    """Exactly equal poles of an exact PoleSpectrum: basis, residue sum, projector, ALP weight (SymPy)."""
    poles = spectrum.pole_squared_masses
    groups = []
    for members in _chain_groups(poles, 0):
        if len(members) < 2:
            continue
        basis = sp.Matrix.hstack(*[spectrum.pole_vector(index) for index in members])
        residue_sum = (basis * basis.T).applyfunc(sp.simplify)
        groups.append({
            "pole": poles[members[0]], "pole_indices": members, "basis": basis,
            "residue_sum": residue_sum,
            "projector": (residue_sum * spectrum.kinetic).applyfunc(sp.simplify),
            "alp_weight": sp.simplify(residue_sum[alp_index, alp_index]),
        })
    return groups


def single_pole_formula_singular(spectrum: PoleSpectrum, relative_tolerance=1e-9, alp_index=-1):
    """For each pole: smallest singular value of M^2_PP - m_n^2 K_PP and whether it vanishes (relative)."""
    kinetic = np.array(spectrum.kinetic, dtype=float)
    mass_squared = np.array(spectrum.mass_squared, dtype=float)
    size = kinetic.shape[0]
    alp = alp_index % size
    meson = [index for index in range(size) if index != alp]
    scale = max(1.0, float(np.max(np.abs(mass_squared))))
    result = []
    for pole in spectrum.pole_squared_masses:
        block = mass_squared[np.ix_(meson, meson)] - float(pole) * kinetic[np.ix_(meson, meson)]
        smallest = float(np.linalg.svd(block, compute_uv=False)[-1])
        result.append({"pole": float(pole), "smallest_singular_value": smallest,
                       "singular": smallest <= relative_tolerance * scale})
    return result


# ------------------------------------------------------- illustrative inputs


def isolated_two_level_inputs(ratio=sp.Rational(1, 10), axial=sp.Rational(1, 2), pseudoscalar=sp.Rational(-1, 4)):
    """LO inputs where the ALP mixes with P_3 only: exact two-level system inside the 6x6.

    chi_0 = diag(chi, chi, chi_s) (isospin limit, so P_3 is a mass eigenstate
    with m_pi^2 = chi), c_q^A = axial lambda_3, c_q^P = pseudoscalar lambda_3,
    c_G = 0. Then K_ia and M^2_ia vanish except for i = 3, and
    M^2_{P'a,can} = 2 epsilon chi (axial - pseudoscalar),
    M^2_aa,can = m_a0^2 + 4 epsilon^2 chi axial (axial - 2 pseudoscalar).
    Returns the inputs with m_a0^2 chosen so that M^2_aa,can = chi exactly.
    """
    chi, chi_strange = sp.Rational(1, 50), sp.Rational(47, 100)
    pion_decay_constant = sp.Rational(23, 250)
    canonical_shift = 4 * ratio**2 * chi * axial * (axial - 2 * pseudoscalar)
    lambda_3 = sp.diag(1, -1, 0)
    return {
        "axial_source": axial * lambda_3,
        "pseudoscalar_source": pseudoscalar * lambda_3,
        "quark_spurion": sp.diag(chi, chi, chi_strange),
        "gluon_coupling": sp.Integer(0),
        "pion_decay_constant": pion_decay_constant,
        "alp_decay_constant": pion_decay_constant / ratio,
        "anomaly_mass_squared": sp.Rational(4, 5),
        "bare_alp_mass_squared": chi - canonical_shift,
    }


def isolated_two_level_expectation(inputs):
    """Exact m_pi^2, M^2_{P'a,can} and the two poles at degeneracy for isolated_two_level_inputs."""
    chi = inputs["quark_spurion"][0, 0]
    ratio = inputs["pion_decay_constant"] / inputs["alp_decay_constant"]
    axial = inputs["axial_source"][0, 0]
    pseudoscalar = inputs["pseudoscalar_source"][0, 0]
    mixing = 2 * ratio * chi * (axial - pseudoscalar)
    return {"meson_mass_squared": chi, "canonical_mixing": mixing,
            "poles": (chi - abs(mixing), chi + abs(mixing)), "splitting": 2 * abs(mixing)}


def lo_blocks_float(inputs):
    kinetic, mass_squared = closed_form_blocks(**inputs)
    return np.array(kinetic.evalf(30), dtype=float), np.array(mass_squared.evalf(30), dtype=float)


# --------------------------------------------------------------------- run


def derive():
    """Run the stability and degeneracy checks; returns named groups of results."""
    groups = {"2.2 LO kinetic matrix": lo_kinetic_matrix_checks()["checks"]}
    for source_structure in ("flavor_diagonal", "neutral"):
        inputs = illustrative_lo_inputs(source_structure)
        kinetic, mass_squared = lo_blocks_float(inputs)
        _, diagnostics = solve_physical_spectrum(kinetic, mass_squared, field_labels=NEUTRAL_FIELD_LABELS)
        bound = lo_bare_mass_stability_bound(inputs)
        checks = {
            "illustrative point is stable": diagnostics.stable,
            "det K = 1 (floating point)": abs(diagnostics.kinetic_determinant - 1) < 1e-12,
            "M^2_PP positive definite (exact)": bound["meson_block_positive_definite"],
            "m_a0^2 bound (exact)": bound["bound"],
        }
        for label, offset in (("just above", sp.Rational(1, 10**6)), ("just below", -sp.Rational(1, 10**6))):
            shifted = dict(inputs, bare_alp_mass_squared=bound["bound"] + offset)
            try:
                solve_physical_spectrum(*lo_blocks_float(shifted))
                stable = True
            except UnstableVacuumError as error:
                stable = False
                checks[f"{label} the bound: violated conditions"] = error.violated_conditions
            checks[f"{label} the bound: stable = {offset > 0}"] = stable == (offset > 0)
        groups[f"2.2 stability, {source_structure} sources"] = checks
    inputs = isolated_two_level_inputs()
    expectation = isolated_two_level_expectation(inputs)
    kinetic, mass_squared = closed_form_blocks(**inputs)
    exact_spectrum, _ = solve_physical_spectrum_exact(kinetic, mass_squared, field_labels=NEUTRAL_FIELD_LABELS)
    report = find_near_degeneracies(solve_spectrum_numerical(*lo_blocks_float(inputs), field_labels=NEUTRAL_FIELD_LABELS))
    alarm = report.alp_meson_near_degeneracies
    groups["2.3 isolated two-level system at degeneracy"] = {
        "exact poles contain m_pi^2 -+ |M^2_P'a,can|": all(pole in exact_spectrum.pole_squared_masses
                                                           for pole in expectation["poles"]),
        "exact splitting = 2|M^2_P'a,can|": expectation["splitting"],
        "one ALP-meson near-degeneracy found": len(alarm) == 1,
        "numerical splitting equals 2|M^2_P'a,can|": bool(alarm) and abs(alarm[0].spread - float(expectation["splitting"])) < 1e-12,
        "mixing angle pi/4 at degeneracy": bool(alarm) and abs(alarm[0].mixing_angle - np.pi / 4) < 1e-9,
        "K0 doublet reported without ALP flag": any(
            group.exactly_degenerate and group.alp_weight < 1e-20 and not group.alp_meson_near_degeneracy
            for group in report.groups),
    }
    return groups


def main():
    for name, checks in derive().items():
        print(name)
        for label, value in checks.items():
            if isinstance(value, bool):
                print(f"  [{'pass' if value else 'FAIL'}] {label}")
            else:
                print(f"  {label} = {value}")
    spectrum = solve_spectrum_numerical(*lo_blocks_float(illustrative_lo_inputs("flavor_diagonal")),
                                        field_labels=NEUTRAL_FIELD_LABELS)
    report = find_near_degeneracies(spectrum)
    print("Degeneracies, flavour-diagonal illustrative sources (relative tolerance 1e-2):")
    print(f"  minimal gap {report.minimal_gap:.6g} GeV^2 between poles {report.minimal_gap_indices}")
    for group in report.groups:
        print(f"  poles {group.poles}: exact {group.exactly_degenerate}, ALP weight {group.alp_weight:.3g}, "
              f"ALP-carrying states {group.alp_carrying_state_count}, ALP-meson flag {group.alp_meson_near_degeneracy}")


if __name__ == "__main__":
    main()
