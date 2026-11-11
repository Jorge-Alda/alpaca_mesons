"""Field-redefinition invariant photon coefficient dressed by the axial connection.

Run (from the repository root):
    .venv/bin/python scripts/derive_alp_dressed_photon_invariant.py

Conventions and units
---------------------
Natural units; x = a/f_a is dimensionless and real. Flavor matrices are
3x3 in the basis (u, d, s) and noncommuting. The quark charge matrix is
Q = diag(2/3, -1/3, -1/3) in units of the positron charge e, and N_c = 3.
The quark mass function M(x) multiplies bar(q_R) M q_L. C_R(x), C_L(x) are
the dimensionless Hermitian coefficients of the derivative currents
(d_mu a / f_a) bar(q_eta) gamma^mu C_eta q_eta, C_A = (C_R - C_L)/2.
Theta_gamma(x) is the coefficient of the photon operator
    L_gg = (alpha_em / 4 pi) Theta_gamma(x) F_{mu nu} tilde(F)^{mu nu},
so Theta_{gamma,1} = c_gamma (Bauer:2020jbp normalization, Alda:2026nrd rotation
signs). Taylor coefficients X_n are derivatives d^n X/dx^n at x = 0.

The photon couples minimally through bar(q) gamma^mu e Q A_mu q with fixed Q.
Electric-charge conservation of the ALP couplings requires
[Q, M(x)] = [Q, C_R(x)] = [Q, C_L(x)] = 0, i.e. all matrices are block
diagonal in the charge sectors b in {u} and {d, s}. Field redefinitions
q_eta' = g_eta(x) q_eta, g_eta(0) = 1, are restricted to [Q, g_eta] = 0
(the commutant U(1)_u x U(2)_ds), so the photon coupling is unchanged.
Within the ds block the generators do not commute.

Laws used:
    C_eta' = g_eta C_eta g_eta^dagger + i g_eta d_x g_eta^dagger,
    M' = g_R M g_L^dagger,
    Theta_gamma'(x) = Theta_gamma(x) - 2 N_c sum_b Q_b^2 (1/2i) log(det_b g_R / det_b g_L),
with det_b the determinant inside charge block b and the logarithm continuous
from x = 0. For g_eta = exp(i G_eta) this equals
Theta_gamma - N_c tr(Q^2 (G_R - G_L)).

Derivation of the finite anomaly law (assumption extending the linear law of
the Alda:2026nrd photon-normalization extract). The extract gives, for
g_{R,L} = exp[i x (delta +- kappa)], c_gamma -> c_gamma - 2 N_c tr(Q^2 kappa),
i.e. the infinitesimal law
    d Theta_gamma = -N_c tr[Q^2 (d alpha_R - d alpha_L)],
for a rotation g -> (1 + i d alpha) g with [Q, d alpha] = 0. Writing a finite
path of rotations g(x) and d alpha_eta = -i (d g_eta) g_eta^dagger, the
identity tr_b(g d g^dagger) = d log det_b g^dagger in each block gives
tr(Q^2 d alpha) = sum_b Q_b^2 (-i) d log det_b g, and integrating along x
yields the finite law above. This uses that Q (and hence the anomaly
coefficient) does not change along the path, which holds only for [Q, g] = 0.

Claim verified here: the dressed photon coefficient
    hat(Theta_gamma)(x) = Theta_gamma(x) + 2 N_c int_0^x tr(Q^2 C_A(x')) dx'
is invariant under all such commutant redefinitions, with Taylor
coefficients hat(Theta)_{gamma,n} = Theta_{gamma,n} + 2 N_c tr(Q^2 C_{A,n-1}) and
hat(Theta)_{gamma,1} = c_gamma + 2 N_c tr(Q^2 c^A) = c_gammaA.

Inputs
------
* memory/extracts/Alda:2026nrd/photon-anomaly-normalization.md: photon
  operator normalization, linear anomaly law, and c_gammaA.
* memory/alp/shared/alp-general-chiral-redefinition-invariants.md:
  connection law, transport W_eta (dW/dx = i C W, W(0) = 1), hat(M).
* memory/alp/shared/alp-dressed-topological-invariant.md: hat(Theta) and the
  helper functions imported from derive_alp_dressed_topological_invariant.py.
* memory/alp/photon/alp-photon-invariant-from-mixing.md: the application
  C_agg = c_gamma + 2 N_c tr[Q^2 (c^A + M_a^{-1}(C))].

Assumptions and approximations
------------------------------
* All series are truncated at O(x^4) (SERIES_ORDER): series checks, not
  all-orders proofs.
* The finite anomaly law is the integrated form of the linear law (see above).
* Only the photon F tilde(F) term with constant Q is tracked; the weak
  interactions and any Z-photon terms are ignored.

Each check prints PASS or FAIL; the process exits with status 1 on failure.
"""

import random
import sys
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from derive_alp_dressed_mass_first_order import random_rational_hermitian  # noqa: E402
from derive_alp_dressed_topological_invariant import (  # noqa: E402
    ALP_FIELD_RATIO,
    SERIES_ORDER,
    connection_transform,
    determinant_series,
    dressed_theta,
    integrate_from_zero,
    logarithm_series,
    matrix_exponential_series,
    matrix_transport_series,
    positive_mass_matrix,
    random_rational_complex,
    report,
    symbolic_hermitian_matrix,
    taylor_coefficients,
    taylor_polynomial,
    truncate,
    truncate_matrix,
)

COLOR_NUMBER = 3
FLAVOR_NUMBER = 3
QUARK_CHARGES = (sp.Rational(2, 3), sp.Rational(-1, 3), sp.Rational(-1, 3))
CHARGE_MATRIX = sp.diag(*QUARK_CHARGES)
CHARGE_BLOCKS = ((0,), (1, 2))  # u and (d, s) sectors
DEFAULT_SEED = 20261111


# ---------------------------------------------------------------------------
# Charge-sector structure
# ---------------------------------------------------------------------------

def project_to_charge_commutant(matrix, charges=QUARK_CHARGES):
    """Remove entries connecting quarks of different charge, so [Q, matrix] = 0."""
    return sp.Matrix(matrix.shape[0], matrix.shape[1],
                     lambda row, column: matrix[row, column]
                     if charges[row] == charges[column] else 0)


def commutator(first, second):
    return first * second - second * first


def is_zero_matrix(matrix):
    return matrix.applyfunc(sp.expand) == sp.zeros(*matrix.shape)


def block(matrix, indices):
    return matrix.extract(list(indices), list(indices))


def random_charge_conserving_theory(seed=DEFAULT_SEED, size=FLAVOR_NUMBER):
    """C_R(x), C_L(x), M(x) commuting with Q, plus Theta_gamma(x) and Theta(x) coefficients."""
    generator = random.Random(seed)

    def hermitian():
        return project_to_charge_commutant(random_rational_hermitian(generator, size))

    right_coefficients = [hermitian() for _ in range(SERIES_ORDER + 1)]
    left_coefficients = [hermitian() for _ in range(SERIES_ORDER + 1)]
    mass_coefficients = ([project_to_charge_commutant(positive_mass_matrix(generator, size))]
                         + [project_to_charge_commutant(random_rational_complex(generator, size))
                            for _ in range(SERIES_ORDER)])

    def scalars():
        return [sp.Rational(generator.randint(-9, 9), generator.randint(1, 5))
                for _ in range(SERIES_ORDER + 1)]

    return {"right": right_coefficients, "left": left_coefficients,
            "mass": mass_coefficients, "theta_gamma": scalars(), "theta": scalars()}


def random_redefinition(seed, mode="nonlinear", size=FLAVOR_NUMBER):
    """Hermitian G_R(x), G_L(x) = sum_{n>=1} G_n x^n / n!.

    mode:
      "nonlinear"       generators in the commutant of Q, with extra
                        block-identity pieces at x^2..x^4 (different for R, L);
      "traceless"       commutant generators with tr G_eta = 0 (tr(Q^2 G) != 0);
      "block_traceless" commutant generators traceless inside each charge block,
                        so Theta_gamma itself is unchanged;
      "charge_mixing"   unrestricted generators, with u-d entries ([Q, g] != 0).
    """
    generator = random.Random(seed)
    generators = {}
    for chirality in ("right", "left"):
        coefficients = [sp.zeros(size, size)]
        for degree in range(1, SERIES_ORDER + 1):
            coefficient = random_rational_hermitian(generator, size)
            if mode != "charge_mixing":
                coefficient = project_to_charge_commutant(coefficient)
            if mode == "nonlinear" and degree >= 2:
                for indices in CHARGE_BLOCKS:
                    weight = sp.Rational(generator.randint(-9, 9), generator.randint(1, 5))
                    for index in indices:
                        coefficient[index, index] += weight
            elif mode == "traceless":
                coefficient = coefficient - coefficient.trace() / size * sp.eye(size)
            elif mode == "block_traceless":
                for indices in CHARGE_BLOCKS:
                    block_trace = sum(coefficient[index, index] for index in indices)
                    for index in indices:
                        coefficient[index, index] -= block_trace / len(indices)
            coefficients.append(coefficient)
        generators[chirality] = coefficients
    return generators["right"], generators["left"]


# ---------------------------------------------------------------------------
# Photon coefficient laws
# ---------------------------------------------------------------------------

def axial_charge_trace_integral(right_series, left_series, charge_matrix=CHARGE_MATRIX):
    """int_0^x tr(Q^2 C_A)."""
    return sp.expand(integrate_from_zero(
        (charge_matrix**2 * (right_series - left_series)).trace() / 2))


def dressed_theta_gamma(theta_gamma_series, right_series, left_series,
                        charge_matrix=CHARGE_MATRIX):
    """hat(Theta_gamma) = Theta_gamma + 2 N_c int_0^x tr(Q^2 C_A)."""
    return truncate(theta_gamma_series + 2 * COLOR_NUMBER
                    * axial_charge_trace_integral(right_series, left_series, charge_matrix))


def block_anomaly_shift(right_group, left_group):
    """2 N_c sum_b Q_b^2 (1/2i) log(det_b g_R / det_b g_L), continuous branch."""
    shift = 0
    for indices in CHARGE_BLOCKS:
        charge = QUARK_CHARGES[indices[0]]
        log_ratio = (logarithm_series(determinant_series(block(right_group, indices)))
                     - logarithm_series(determinant_series(block(left_group, indices))))
        shift += 2 * COLOR_NUMBER * charge**2 * log_ratio / (2 * sp.I)
    return sp.expand(shift)


def linear_trace_anomaly_shift(right_generator, left_generator):
    """N_c tr(Q^2 (G_R - G_L)), the trace form of the anomaly shift."""
    return truncate(COLOR_NUMBER * (CHARGE_MATRIX**2 * (right_generator - left_generator)).trace())


# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------

def check_charge_preservation(seed=DEFAULT_SEED):
    """Check 1: [Q, W_eta] = 0 and W_eta^dagger Q W_eta = Q, also [Q, hat M] = 0, O(x^4)."""
    theory = random_charge_conserving_theory(seed)
    transports = {chirality: matrix_transport_series(theory[chirality])
                  for chirality in ("right", "left")}
    results = {}
    for chirality, transport in transports.items():
        results[f"{chirality}_commutes_with_charge"] = is_zero_matrix(
            truncate_matrix(commutator(CHARGE_MATRIX, transport)))
        results[f"{chirality}_preserves_charge"] = is_zero_matrix(
            truncate_matrix(transport.H * CHARGE_MATRIX * transport) - CHARGE_MATRIX)
    dressed_mass = truncate_matrix(transports["right"].H * taylor_polynomial(theory["mass"])
                                   * transports["left"])
    results["dressed_mass_commutes_with_charge"] = is_zero_matrix(
        truncate_matrix(commutator(CHARGE_MATRIX, dressed_mass)))
    return results


def check_redefinition_invariance(seed=DEFAULT_SEED, redefinition_seed=7, mode="nonlinear"):
    """Check 2: hat(Theta_gamma) invariant under commutant g_eta = exp(i G_eta), O(x^4)."""
    theory = random_charge_conserving_theory(seed)
    right_series = taylor_polynomial(theory["right"])
    left_series = taylor_polynomial(theory["left"])
    theta_gamma_series = taylor_polynomial(theory["theta_gamma"])
    theta_series = taylor_polynomial(theory["theta"])
    right_coefficients, left_coefficients = random_redefinition(redefinition_seed, mode)
    right_generator = taylor_polynomial(right_coefficients)
    left_generator = taylor_polynomial(left_coefficients)
    right_group = matrix_exponential_series(right_generator)
    left_group = matrix_exponential_series(left_generator)

    right_series_new = connection_transform(right_group, right_series)
    left_series_new = connection_transform(left_group, left_series)
    anomaly_shift = block_anomaly_shift(right_group, left_group)
    theta_gamma_series_new = truncate(theta_gamma_series - anomaly_shift)
    theta_series_new = truncate(theta_series - sp.expand(
        (logarithm_series(determinant_series(right_group))
         - logarithm_series(determinant_series(left_group))) / (2 * sp.I)))

    # i tr_b(g d_x g^dagger) = d_x tr_b G in each charge block.
    block_trace_residuals = []
    for group, generator_series in ((right_group, right_generator),
                                    (left_group, left_generator)):
        inhomogeneous = truncate_matrix(sp.I * group * sp.diff(group.H, ALP_FIELD_RATIO),
                                        SERIES_ORDER - 1)
        for indices in CHARGE_BLOCKS:
            block_trace_residuals.append(truncate(
                block(inhomogeneous, indices).trace()
                - sp.diff(block(generator_series, indices).trace(), ALP_FIELD_RATIO),
                SERIES_ORDER - 1))

    right_transport_new = matrix_transport_series(taylor_coefficients(right_series_new))
    left_transport_new = matrix_transport_series(taylor_coefficients(left_series_new))
    covariance = [truncate_matrix(right_transport_new
                                  - right_group * matrix_transport_series(theory["right"])),
                  truncate_matrix(left_transport_new
                                  - left_group * matrix_transport_series(theory["left"]))]
    ds_indices = CHARGE_BLOCKS[1]
    return {
        "generators_commute_with_charge": all(
            is_zero_matrix(commutator(CHARGE_MATRIX, coefficient))
            for coefficient in right_coefficients + left_coefficients),
        "ds_block_noncommuting": not is_zero_matrix(commutator(
            block(right_coefficients[1], ds_indices), block(right_coefficients[2], ds_indices))),
        "anomaly_law_residual": sp.expand(
            anomaly_shift - linear_trace_anomaly_shift(right_generator, left_generator)),
        "block_trace_residuals": block_trace_residuals,
        "theta_gamma_changes": sp.expand(theta_gamma_series_new - theta_gamma_series) != 0,
        "transport_covariance_is_zero": all(is_zero_matrix(matrix) for matrix in covariance),
        "dressed_theta_gamma_residual": sp.expand(
            dressed_theta_gamma(theta_gamma_series_new, right_series_new, left_series_new)
            - dressed_theta_gamma(theta_gamma_series, right_series, left_series)),
        "dressed_theta_residual": sp.expand(
            dressed_theta(theta_series_new, right_series_new, left_series_new)
            - dressed_theta(theta_series, right_series, left_series)),
    }


def check_taylor_form(order=6):
    """Check 3: Taylor coefficients, hat(Theta)_{gamma,1} = c_gammaA, its invariance,
    explicit form."""
    theta_gamma_coefficients = sp.symbols(f"Theta_gamma_0:{order + 1}")
    charge_trace_coefficients = sp.symbols(f"trQ2C_A0:{order + 1}")
    dressed = (taylor_polynomial(list(theta_gamma_coefficients)) + 2 * COLOR_NUMBER
               * integrate_from_zero(taylor_polynomial(list(charge_trace_coefficients))))
    taylor_residuals = [sp.expand(
        sp.diff(dressed, ALP_FIELD_RATIO, degree).subs(ALP_FIELD_RATIO, 0)
        - theta_gamma_coefficients[degree]
        - (2 * COLOR_NUMBER * charge_trace_coefficients[degree - 1] if degree else 0))
        for degree in range(order + 1)]

    # First order from the series: Theta_{gamma,1} = c_gamma, C_A(x) = c^A + O(x).
    photon_coupling = sp.Symbol("c_gamma", real=True)
    axial_coupling = symbolic_hermitian_matrix("cA")
    axial_generator = symbolic_hermitian_matrix("kappa")
    theta_gamma_series = photon_coupling * ALP_FIELD_RATIO
    first_order = sp.diff(dressed_theta_gamma(theta_gamma_series, axial_coupling, -axial_coupling),
                          ALP_FIELD_RATIO).subs(ALP_FIELD_RATIO, 0)
    photon_axial_invariant = (photon_coupling
                              + 2 * COLOR_NUMBER * (CHARGE_MATRIX**2 * axial_coupling).trace())
    shifted = ((photon_coupling
                - 2 * COLOR_NUMBER * (CHARGE_MATRIX**2 * axial_generator).trace())
               + 2 * COLOR_NUMBER * (CHARGE_MATRIX**2 * (axial_coupling + axial_generator)).trace())
    explicit = (photon_coupling + sp.Rational(8, 3) * axial_coupling[0, 0]
                + sp.Rational(2, 3) * (axial_coupling[1, 1] + axial_coupling[2, 2]))
    diagonal_symbols = {axial_coupling[index, index] for index in range(FLAVOR_NUMBER)}
    off_diagonal_symbols = axial_coupling.free_symbols - diagonal_symbols
    return {"taylor_residuals": taylor_residuals,
            "first_order_residual": sp.expand(first_order - photon_axial_invariant),
            "shift_residual": sp.expand(shifted - photon_axial_invariant),
            "explicit_form_residual": sp.expand(first_order - explicit),
            "only_diagonal_entries": not (sp.expand(first_order).free_symbols
                                          & off_diagonal_symbols)}


def check_charge_mixing_failure(seed=DEFAULT_SEED, redefinition_seed=7):
    """Check 4: a redefinition with u-d mixing ([Q, g] != 0) breaks the construction.

    The anomaly law in block form is undefined for such g; the naive
    candidate uses the trace form Theta_gamma' = Theta_gamma - N_c tr(Q^2 (G_R - G_L)).
    The residual of hat(Theta_gamma) is split into the part from the homogeneous
    term g C g^dagger and the part from the inhomogeneous term i g d g^dagger.
    """
    theory = random_charge_conserving_theory(seed)
    right_series = taylor_polynomial(theory["right"])
    left_series = taylor_polynomial(theory["left"])
    theta_gamma_series = taylor_polynomial(theory["theta_gamma"])
    theta_series = taylor_polynomial(theory["theta"])
    right_coefficients, left_coefficients = random_redefinition(redefinition_seed,
                                                                "charge_mixing")
    right_generator = taylor_polynomial(right_coefficients)
    left_generator = taylor_polynomial(left_coefficients)
    right_group = matrix_exponential_series(right_generator)
    left_group = matrix_exponential_series(left_generator)

    # q = g^dagger q', so bar(q) Q q = bar(q') g Q g^dagger q': new chiral charges.
    right_charge_new = truncate_matrix(right_group * CHARGE_MATRIX * right_group.H)
    left_charge_new = truncate_matrix(left_group * CHARGE_MATRIX * left_group.H)

    right_series_new = connection_transform(right_group, right_series)
    left_series_new = connection_transform(left_group, left_series)
    theta_gamma_series_new = truncate(theta_gamma_series
                                - linear_trace_anomaly_shift(right_generator, left_generator))
    theta_series_new = truncate(theta_series - sp.expand(
        (logarithm_series(determinant_series(right_group))
         - logarithm_series(determinant_series(left_group))) / (2 * sp.I)))

    homogeneous_residual = 2 * COLOR_NUMBER * axial_charge_trace_integral(
        truncate_matrix(right_group * right_series * right_group.H),
        truncate_matrix(left_group * left_series * left_group.H))
    homogeneous_residual = truncate(homogeneous_residual - 2 * COLOR_NUMBER
                                    * axial_charge_trace_integral(right_series, left_series))
    inhomogeneous = [truncate_matrix(sp.I * group * sp.diff(group.H, ALP_FIELD_RATIO))
                     for group in (right_group, left_group)]
    inhomogeneous_residual = truncate(
        2 * COLOR_NUMBER * axial_charge_trace_integral(*inhomogeneous)
        - linear_trace_anomaly_shift(right_generator, left_generator))
    dressed_theta_gamma_residual = sp.expand(
        dressed_theta_gamma(theta_gamma_series_new, right_series_new, left_series_new)
        - dressed_theta_gamma(theta_gamma_series, right_series, left_series))
    return {
        "right_charge_changed": not is_zero_matrix(right_charge_new - CHARGE_MATRIX),
        "left_charge_changed": not is_zero_matrix(left_charge_new - CHARGE_MATRIX),
        "charges_chiral": not is_zero_matrix(right_charge_new - left_charge_new),
        "charge_x_dependent": sp.diff(right_charge_new, ALP_FIELD_RATIO) != sp.zeros(
            FLAVOR_NUMBER, FLAVOR_NUMBER),
        "couplings_violate_charge": not is_zero_matrix(
            truncate_matrix(commutator(CHARGE_MATRIX, right_series_new))),
        "dressed_theta_gamma_residual": dressed_theta_gamma_residual,
        "lowest_residual_power": min(
            degree for (degree,)
            in sp.Poly(dressed_theta_gamma_residual, ALP_FIELD_RATIO).monoms()),
        "homogeneous_residual": homogeneous_residual,
        "inhomogeneous_residual": inhomogeneous_residual,
        "dressed_theta_residual": sp.expand(
            dressed_theta(theta_series_new, right_series_new, left_series_new)
            - dressed_theta(theta_series, right_series, left_series)),
    }


def check_independence_from_gluons(seed=DEFAULT_SEED):
    """Check 5: hat(Theta_gamma) and hat(Theta) weight the diagonal of C_A differently."""
    quark_charges = sp.symbols("Q_u Q_d Q_s", real=True)
    weight_map = sp.Matrix([[1, 1, 1],
                            [2 * COLOR_NUMBER * charge**2 for charge in QUARK_CHARGES]])
    symbolic_weight_map = sp.Matrix([[1, 1, 1],
                                     [2 * COLOR_NUMBER * charge**2 for charge in quark_charges]])
    minors = [sp.factor(symbolic_weight_map.extract([0, 1], [first, second]).det())
              for first, second in ((0, 1), (0, 2), (1, 2))]
    ds_weight_map = weight_map.extract([0, 1], [1, 2])

    # Redefinition shifts of (c_G, c_gamma) for diagonal kappa: (-tr kappa, -2N_c tr(Q^2 kappa)).
    shift_map = -weight_map

    # With C_A supported in the ds block only,
    # hat(Theta_gamma) - Theta_gamma = 2 N_c Q_d^2 (hat(Theta) - Theta).
    theory = random_charge_conserving_theory(seed)
    right_series = taylor_polynomial(theory["right"])
    left_series = taylor_polynomial(theory["left"])
    theta_gamma_series = taylor_polynomial(theory["theta_gamma"])
    theta_series = taylor_polynomial(theory["theta"])
    ds_projector = sp.diag(0, 1, 1)
    ds_right = ds_projector * right_series * ds_projector
    ds_left = ds_projector * left_series * ds_projector
    ds_proportionality = sp.expand(
        (dressed_theta_gamma(theta_gamma_series, ds_right, ds_left) - theta_gamma_series)
        - 2 * COLOR_NUMBER * QUARK_CHARGES[1]**2
        * (dressed_theta(theta_series, ds_right, ds_left) - theta_series))
    general_difference = sp.expand(
        (dressed_theta_gamma(theta_gamma_series, right_series, left_series) - theta_gamma_series)
        - 2 * COLOR_NUMBER * QUARK_CHARGES[1]**2
        * (dressed_theta(theta_series, right_series, left_series) - theta_series))
    return {"weight_rank": weight_map.rank(),
            "symbolic_minors": minors,
            "ds_weight_rank": ds_weight_map.rank(),
            "shift_rank": shift_map.rank(),
            "ds_proportionality_residual": ds_proportionality,
            "full_not_proportional": general_difference != 0}


def check_application_identity():
    """Check 6: c_gamma + 2N_c tr[Q^2 (c^A + X)] = hat(Theta)_{gamma,1} + 2N_c tr(Q^2 X)."""
    photon_coupling = sp.Symbol("c_gamma", real=True)
    axial_coupling = symbolic_hermitian_matrix("cA")
    axial_generator = symbolic_hermitian_matrix("kappa")
    invariant_matrix = symbolic_hermitian_matrix("X")
    dressed_one = photon_coupling + 2 * COLOR_NUMBER * (CHARGE_MATRIX**2 * axial_coupling).trace()

    def physical(photon, axial):
        return photon + 2 * COLOR_NUMBER * (CHARGE_MATRIX**2 * (axial + invariant_matrix)).trace()

    shifted = physical(photon_coupling
                       - 2 * COLOR_NUMBER * (CHARGE_MATRIX**2 * axial_generator).trace(),
                       axial_coupling + axial_generator)
    return {"identity_residual": sp.expand(
                physical(photon_coupling, axial_coupling) - dressed_one
                - 2 * COLOR_NUMBER * (CHARGE_MATRIX**2 * invariant_matrix).trace()),
            "shift_residual": sp.expand(shifted - physical(photon_coupling, axial_coupling))}


def main():
    results = []
    charge = check_charge_preservation()
    for chirality in ("right", "left"):
        results.append(report(f"1 {chirality}: [Q, W] = 0 to O(x^4)",
                              charge[f"{chirality}_commutes_with_charge"]))
        results.append(report(f"1 {chirality}: W^dagger Q W = Q to O(x^4)",
                              charge[f"{chirality}_preserves_charge"]))
    results.append(report("1 [Q, hat M] = 0 to O(x^4)", charge["dressed_mass_commutes_with_charge"]))

    for mode, theta_gamma_expected in (("nonlinear", True), ("traceless", True),
                                 ("block_traceless", False)):
        values = check_redefinition_invariance(mode=mode)
        results.append(report(f"2 {mode}: [Q, G_n] = 0", values["generators_commute_with_charge"]))
        results.append(report(f"2 {mode}: ds-block generators noncommuting",
                              values["ds_block_noncommuting"]))
        results.append(report(f"2 {mode}: i tr_b(g d g^dagger) = d tr_b G in each block",
                              all(residual == 0 for residual in values["block_trace_residuals"])))
        results.append(report(f"2 {mode}: 2N_c sum_b Q_b^2 (1/2i) log(det_b g_R/det_b g_L)"
                              " = N_c tr(Q^2(G_R-G_L))", values["anomaly_law_residual"] == 0))
        results.append(report(f"2 {mode}: Theta_gamma changes = {values['theta_gamma_changes']}"
                              f" (expected {theta_gamma_expected})",
                              values["theta_gamma_changes"] == theta_gamma_expected))
        results.append(report(f"2 {mode}: W' = g W (R and L) to O(x^4)",
                              values["transport_covariance_is_zero"]))
        results.append(report(f"2 {mode}: hat Theta_gamma' = hat Theta_gamma to O(x^4)",
                              values["dressed_theta_gamma_residual"] == 0))
        results.append(report(f"2 {mode}: hat Theta' = hat Theta to O(x^4)",
                              values["dressed_theta_residual"] == 0))

    taylor = check_taylor_form()
    results.append(report("3 hat Theta_{gamma,n} = Theta_{gamma,n} + 2N_c tr(Q^2 C_{A,n-1})"
                          " (n = 0..6)",
                          all(residual == 0 for residual in taylor["taylor_residuals"])))
    results.append(report("3 hat Theta_{gamma,1} = c_gamma + 2N_c tr(Q^2 c^A) = c_gammaA",
                          taylor["first_order_residual"] == 0))
    results.append(report("3 hat Theta_{gamma,1} invariant under"
                          " c_gamma -> c_gamma - 2N_c tr(Q^2 kappa), c^A -> c^A + kappa",
                          taylor["shift_residual"] == 0))
    results.append(report("3 hat Theta_{gamma,1} = c_gamma + (8/3) c^A_uu + (2/3)(c^A_dd + c^A_ss)",
                          taylor["explicit_form_residual"] == 0))
    results.append(report("3 only diagonal entries of c^A enter", taylor["only_diagonal_entries"]))

    failure = check_charge_mixing_failure()
    results.append(report("4 u-d mixing: hat Theta_gamma residual starts at x^3"
                          " (hat Theta_{gamma,1} = c_gammaA and hat Theta_{gamma,2} unaffected)",
                          failure["lowest_residual_power"] == 3))
    results.append(report("4 u-d mixing: g_R Q g_R^dagger != Q", failure["right_charge_changed"]))
    results.append(report("4 u-d mixing: g_L Q g_L^dagger != Q", failure["left_charge_changed"]))
    results.append(report("4 u-d mixing: new charges differ for R and L (chiral photon coupling)",
                          failure["charges_chiral"]))
    results.append(report("4 u-d mixing: new charges depend on x", failure["charge_x_dependent"]))
    results.append(report("4 u-d mixing: [Q, C_R'] != 0", failure["couplings_violate_charge"]))
    results.append(report("4 u-d mixing: naive hat Theta_gamma not invariant",
                          failure["dressed_theta_gamma_residual"] != 0))
    results.append(report("4 u-d mixing: homogeneous part tr(Q^2 g C g^dagger) - tr(Q^2 C) != 0",
                          failure["homogeneous_residual"] != 0))
    results.append(report("4 u-d mixing: inhomogeneous part != N_c tr(Q^2(G_R-G_L))",
                          failure["inhomogeneous_residual"] != 0))
    results.append(report("4 u-d mixing: hat Theta still invariant",
                          failure["dressed_theta_residual"] == 0))

    independence = check_independence_from_gluons()
    results.append(report("5 rank of weight map {(1,1,1), 2N_c Q^2} = 2",
                          independence["weight_rank"] == 2))
    print(f"    5 2x2 minors for general charges: {independence['symbolic_minors']}")
    results.append(report("5 rank of (c_G, c_gamma) shift map = 2", independence["shift_rank"] == 2))
    results.append(report("5 ds block alone: weight map rank 1 (proportional)",
                          independence["ds_weight_rank"] == 1))
    results.append(report("5 ds support: hat Theta_gamma - Theta_gamma"
                          " = (2/3)(hat Theta - Theta) to O(x^4)",
                          independence["ds_proportionality_residual"] == 0))
    results.append(report("5 full support: not proportional", independence["full_not_proportional"]))

    application = check_application_identity()
    results.append(report("6 c_gamma + 2N_c tr[Q^2(c^A+X)] = hat Theta_{gamma,1} + 2N_c tr(Q^2 X)",
                          application["identity_residual"] == 0))
    results.append(report("6 invariant under the delta, kappa shifts with X fixed",
                          application["shift_residual"] == 0))

    if not all(results):
        sys.exit(1)


if __name__ == "__main__":
    main()
