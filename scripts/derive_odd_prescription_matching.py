"""Derive neutral Aloni/Bai pole projections and independent-source matching.

Run with ``.venv/bin/python scripts/derive_odd_prescription_matching.py``.
All eleven odd coefficients are dimensionless; F, F_V, M_0 and M_V are in
GeV. The chiral octet OPE is tree level, leading in 1/N_C, with one vector
nonet. M_V normalizes the VJP operators and M_0 is the chiral octet pole.
This script does not identify VMD with the QCD short-distance prescription.
The phenomenological suppression factor of Aloni:2018vki is set to one.
Neutral projections assume commuting flavor matrices, common vector mass,
and independent pseudoscalar momentum and quark-mass-spurion inputs.
Use --derive-only to perform algebra before writing the recorded JSON.
"""

from pathlib import Path
import hashlib
import json
import argparse

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
COEFFICIENT_NAMES = [f"c{index}" for index in range(1, 8)] + [
    f"d{index}" for index in range(1, 5)
]


def derive_neutral_pole_matching():
    """Solve neutral pole projections with independent mass-spurion inputs.

    The physical inputs G_PVgamma and G_PVV multiply the canonical tensor
    epsilon(epsilon_1,epsilon_2,q_1,q_2); they have units GeV**(-1).
    This operation does not equate off-shell Proca and tensor vertices.
    """
    coefficients = sp.symbols(" ".join(COEFFICIENT_NAMES))
    c1, c2, c3, c4, c5, c6, c7, d1, d2, d3, d4 = coefficients
    mass_squared, first_virtuality, second_virtuality, spurion_ratio = sp.symbols(
        "m_P_squared s t F_H_over_F_Q"
    )
    F, M_V, first_mass, second_mass, electric_charge = sp.symbols(
        "F M_V M_i M_j e", nonzero=True
    )
    radiative_coupling, vector_coupling = sp.symbols("G_PVgamma G_PVV")
    vjp = ((c1+c2-c5)*mass_squared
           +(-c1+c2+c5-2*c6)*first_virtuality
           +(c1-c2+c5)*second_virtuality+4*c3*spurion_ratio)
    vvp = ((d1-d3)*mass_squared+d3*(first_virtuality+second_virtuality)
           +4*d2*spurion_ratio)
    # F_Q = tr({T_i,T_j}T_P), and both tensor Wick assignments are included.
    direct_radiative = 2*sp.sqrt(2)*electric_charge*vjp.subs(
        {first_virtuality: first_mass**2, second_virtuality: 0}
    )/(F*first_mass*M_V)
    direct_vectors = -2*sp.sqrt(2)*vvp.subs(
        {first_virtuality: first_mass**2, second_virtuality: second_mass**2}
    )/(F*first_mass*second_mass)
    equations = []
    for expression in (direct_radiative-radiative_coupling,
                       direct_vectors-vector_coupling):
        equations.extend(sp.Poly(expression, mass_squared, spurion_ratio).coeffs())
    matrix, right_hand_side = sp.linear_eq_to_matrix(equations, coefficients)
    solution = next(iter(sp.linsolve((matrix, right_hand_side), coefficients)))
    substitutions = dict(zip(coefficients, solution))
    residuals = [sp.simplify(value.subs(substitutions, simultaneous=True))
                 for value in equations]
    # An arbitrary nonzero constant Proca PVV vertex cannot be represented by
    # the tensor numerator H as an off-shell polynomial without local terms.
    proca_numerator = sp.symbols("H_Proca", nonzero=True)
    offshell_equations = sp.Poly(vvp-proca_numerator, mass_squared,
                                first_virtuality, second_virtuality,
                                spurion_ratio).coeffs()
    offshell_solution = sp.solve(offshell_equations, [d1, d2, d3], dict=True)
    return {
        "coefficients": coefficients, "vjp": vjp, "vvp": vvp,
        "direct_radiative": direct_radiative, "direct_vectors": direct_vectors,
        "equations": equations, "matrix": matrix,
        "solution": solution, "residuals": residuals,
        "offshell_equations": offshell_equations,
        "offshell_solution": offshell_solution,
    }


def derive_bai_omega_gamma_components():
    """Combine the source and anomaly terms printed in Bai:2025fvl.

    The paper's c_eff multiplies a epsilon(partial omega,partial gamma).
    No vector-photon conversion contribution is added to these local terms.
    """
    e, g, N_C, f_a, c_gg = sp.symbols("e g N_C f_a c_gg", nonzero=True)
    Q_d, Q_u, kappa_d, kappa_u = sp.symbols("Q_d Q_u kappa_d kappa_u")
    c_d_L, c_d_R, c_u_L, c_u_R = sp.symbols("c_d_L c_d_R c_u_L c_u_R")
    wzw = e*g*N_C*(Q_d*(c_d_L-c_d_R)+Q_u*(c_u_L-c_u_R)
                      -2*c_gg*(Q_d*kappa_d+Q_u*kappa_u))/(16*sp.pi**2*f_a)
    anomaly = -e*g*N_C*c_gg*(Q_d*kappa_d+Q_u*kappa_u)/(8*sp.pi**2*f_a)
    return {"wzw": wzw, "anomaly": anomaly,
            "sum_before_meson_poles": sp.factor(wzw+anomaly),
            "unrotated_local_term": sp.factor((wzw+anomaly).subs(
                {kappa_d: 0, kappa_u: 0}))}


def derive_complete_vmd_tff():
    """Match the complete neutral tensor TFF to anomaly-normalized double VMD.

    Common vector mass and residue are explicit approximations. F_H/F_Q
    remains an independent spurion variable. This matches an amplitude;
    it does not establish equality for independent axial external sources.
    """
    matching = derive_neutral_pole_matching()
    coefficients = matching["coefficients"]
    mass_squared, s, t, spurion_ratio = sp.symbols("m_P_squared s t F_H_over_F_Q")
    F, F_V, M_V, N_C = sp.symbols("F F_V M_V N_C", nonzero=True)
    anomaly = -sp.sqrt(2)*N_C/(8*sp.pi**2*F)
    first_denominator, second_denominator = M_V**2-s, M_V**2-t
    first_c = matching["vjp"]
    second_c = first_c.xreplace({s: t, t: s})
    numerator = (anomaly*first_denominator*second_denominator
                 +4*F_V/(F*M_V)*(first_c*second_denominator
                                         +second_c*first_denominator)
                 -4*sp.sqrt(2)*F_V**2/F*matching["vvp"])
    requested_numerator = anomaly*M_V**4
    equations = sp.Poly(sp.expand(numerator-requested_numerator),
                        mass_squared, s, t, spurion_ratio).coeffs()
    matrix, rhs = sp.linear_eq_to_matrix(equations, coefficients)
    solutions = sp.linsolve((matrix, rhs), coefficients)
    solution = next(iter(solutions))
    residual = sp.factor((numerator-requested_numerator).subs(
        dict(zip(coefficients, solution)), simultaneous=True))
    return {"numerator": numerator, "requested_numerator": requested_numerator,
            "equations": equations, "rank": matrix.rank(),
            "solution": solution, "residual": residual}


def derive_neutral_source_actions():
    """Expand the commuting neutral cubic forms in the two source actions.

    Differential-form signs use SymPy permutation signatures. The first
    WZW group is i(dv dU† v-dv dU v), the second is
    i[(dv v+v dv) alpha+(dv v+v dv) beta]. The only pure
    background contribution linear in b and quadratic in v comes from
    -(dv A_L+A_L dv)A_R+(dv A_R+A_R dv)A_L.
    """
    from sympy.combinatorics import Permutation

    C, F, g, f_a, N_C = sp.symbols("C F g f_a N_C", nonzero=True)
    goldstone_derivative = sp.sqrt(2)/F
    # One swaps two one-forms when ordering v dPhi as dPhi v.
    swap_one_forms = Permutation([1, 0]).signature()
    first_wzw_group = sp.I*(-sp.I*goldstone_derivative
                              -sp.I*goldstone_derivative)
    second_wzw_group = sp.I*4*sp.I*goldstone_derivative*swap_one_forms
    coefficient_dv_dphi_v = sp.simplify(C*(first_wzw_group+second_wzw_group))
    # d(Phi dv v)=dPhi dv v+Phi dv dv, neglect a total derivative.
    coefficient_phi_dv_dv = -coefficient_dv_dphi_v
    background_vector_normalization = g/sp.sqrt(2)
    bai_pvv = sp.simplify(coefficient_phi_dv_dv*background_vector_normalization**2)
    # A_L=v+b, A_R=v-b; v^2=b^2=0 in the commuting neutral cubic sector.
    coefficient_al_ar_vb = swap_one_forms-1
    coefficient_dv_v_b = sp.expand(C*(-2*coefficient_al_ar_vb
                                       +2*(-coefficient_al_ar_vb)))
    # b=(k_L-k_R) da/(2 f_a); moving da through v changes sign,
    # and da dv v=-a dv dv after integration by parts.
    bai_avv = sp.simplify(coefficient_dv_v_b/(2*f_a)
                         *swap_one_forms*(-1)*background_vector_normalization**2)
    g_vvp = 3*g**2/(8*sp.pi**2*F)
    # Aloni: P=Phi/sqrt2, V_Aloni=V_ideal/sqrt2, and
    # V_munu tilde V^munu d4x=2 dV wedge dV.
    aloni_pvv = sp.simplify(g_vvp/sp.Integer(4)
                            /sp.sqrt(2)/2*2)
    anomaly_substitution = {C: N_C/(48*sp.pi**2)}
    return {
        "bai_dv_dphi_v": coefficient_dv_dphi_v,
        "bai_phi_dV_dV": sp.simplify(bai_pvv.subs(anomaly_substitution)),
        "bai_a_dV_dV_axial_matrix": sp.simplify(bai_avv.subs(anomaly_substitution)),
        "aloni_phi_dV_dV": aloni_pvv,
        "bai_to_aloni_pvv_ratio": sp.simplify((bai_pvv/aloni_pvv).subs(anomaly_substitution)),
    }


def derive_source_component_normalization():
    """Contract explicit components, including the two Fourier derivatives.

    One nonzero polarization component suffices because the neutral parity
    odd two-vector tensor is unique. The trace F_Q contains both flavor
    assignments, so a second identical Lorentz contraction is not inserted.
    """
    from itertools import product

    first_momentum = sp.Matrix([1, 0, 0, 0])
    first_polarization = sp.Matrix([0, 1, 0, 0])
    second_momentum = sp.Matrix([0, 0, 1, 0])
    second_polarization = sp.Matrix([0, 0, 0, 1])
    first_strength = (first_momentum*first_polarization.T
                      -first_polarization*first_momentum.T)
    second_strength = (second_momentum*second_polarization.T
                       -second_polarization*second_momentum.T)
    # epsilon_lower=-epsilon_upper for metric (+---).
    canonical_tensor = sum(-sp.LeviCivita(i, j, k, l)
        *first_polarization[i]*second_polarization[j]
        *first_momentum[k]*second_momentum[l]
        for i, j, k, l in product(range(4), repeat=4))
    strength_dual_strength = sum(-sp.LeviCivita(i, j, k, l)
        *first_strength[i, j]*second_strength[k, l]/2
        for i, j, k, l in product(range(4), repeat=4))
    # i**2 accounts for the derivatives in the interaction vertex.
    component_amplitude_ratio = sp.simplify(sp.I**2
                                            *strength_dual_strength/canonical_tensor)
    return {"canonical_tensor": canonical_tensor,
            "strength_dual_strength": strength_dual_strength,
            "component_action_to_G_per_F_Q": component_amplitude_ratio,
            "form_action_to_G_per_F_Q": component_amplitude_ratio/2}


def derive_full_alp_vector_coefficients():
    """Assemble contacts, neutral pseudoscalar poles and rotation anomaly.

    Expressions are per flavor-matrix coefficient; traces with the two
    physical vector generators are taken after this assembly. The scalar
    theta_Pa*T_P denotes one summand, and the expression is linear in it.
    """
    source = derive_neutral_source_actions()
    N_C, g, f_a, F, c_gg = sp.symbols("N_C g f_a F c_gg", nonzero=True)
    meson_overlap = sp.Symbol("sum_P_theta_Pa_T_P")
    axial_coupling = sp.Symbol("k_L_minus_k_R")
    rotation_phase = sp.Symbol("theta_L_minus_theta_R")
    source_contact = source["bai_a_dV_dV_axial_matrix"]
    anomaly = N_C*g**2*c_gg*rotation_phase/(48*sp.pi**2*f_a)
    unrotated = (source["bai_phi_dV_dV"]*meson_overlap
                 +source_contact*axial_coupling)
    rotated = unrotated+anomaly
    # A rotation shifts k_L-k_R by c_gg(theta_L-theta_R) and
    # the pseudoscalar pole overlap by F*c_gg(theta_L-theta_R)/(sqrt2*f_a).
    rotated_substituted = rotated.subs({
        axial_coupling: axial_coupling+c_gg*rotation_phase,
        meson_overlap: meson_overlap+F*c_gg*rotation_phase/(sp.sqrt(2)*f_a),
    }, simultaneous=True)
    return {"aloni_full_coefficient_matrix": source["aloni_phi_dV_dV"]*meson_overlap,
            "bai_unrotated_full_coefficient_matrix": sp.factor(unrotated),
            "bai_rotated_full_coefficient_matrix": sp.factor(rotated),
            "bai_rotation_residual": sp.simplify(rotated_substituted-unrotated),
            "bai_anomalous_rotation_term": anomaly}


def derive_source_pole_families():
    """Specialize the pole family to the literal source normalizations.

    Aloni's photon residue is obtained from its displayed mass mixing and
    canonical Proca mass term, giving V_ideal -> e Q A/(sqrt2 g).
    Bai's fundamental photon occurs in the total connection with e Q A,
    while its vector connection is g V_ideal/sqrt2.
    """
    neutral = derive_neutral_pole_matching()
    sources = derive_neutral_source_actions()
    e, g, F, F_V = sp.symbols("e g F F_V", nonzero=True)
    G_PVgamma, G_PVV = sp.symbols("G_PVgamma G_PVV")
    M_i, M_j, M_V = sp.symbols("M_i M_j M_V", nonzero=True)
    source_amplitudes = {
        "aloni": {G_PVV: sources["aloni_phi_dV_dV"],
                  G_PVgamma: sources["aloni_phi_dV_dV"]*e/(sp.sqrt(2)*g)},
        "bai": {G_PVV: sources["bai_phi_dV_dV"],
                G_PVgamma: sources["bai_phi_dV_dV"]*sp.sqrt(2)*e/g},
    }
    result = {}
    for name, amplitudes in source_amplitudes.items():
        result[name+"_G_PVV_per_F_Q"] = amplitudes[G_PVV]
        result[name+"_G_PVgamma_per_F_Q"] = sp.simplify(amplitudes[G_PVgamma])
        pole_family = sp.Tuple(*[
            sp.simplify(value.subs(amplitudes, simultaneous=True).subs(
                {M_i: M_V, M_j: M_V}, simultaneous=True))
            for value in neutral["solution"]])
        result[name+"_direct_vertex_family"] = pole_family
        exchange = -4*e*F_V*pole_family[9]/(F*M_V)
        physical_amplitudes = dict(amplitudes)
        physical_amplitudes[G_PVgamma] = amplitudes[G_PVgamma]-exchange
        result[name+"_radiative_exchange_G_per_F_Q"] = exchange
        result[name+"_common_mass_pole_family"] = sp.Tuple(*[
            sp.simplify(value.subs(physical_amplitudes, simultaneous=True).subs(
                {M_i: M_V, M_j: M_V}, simultaneous=True))
            for value in neutral["solution"]])
    # Independent ALP source: U_a=-(k_L-k_R)/f_a. F_Q(U_a)=2D_ij.
    N_C = sp.Symbol("N_C", nonzero=True)
    d3 = neutral["coefficients"][9]
    direct_source_d3 = sp.solve(sp.Eq(8*d3, -N_C*g**2/(12*sp.pi**2)), d3)[0]
    result["bai_direct_axial_source_d3"] = direct_source_d3
    result["bai_pole_minus_source_d3"] = sp.simplify(
        result["bai_common_mass_pole_family"][9]-direct_source_d3)
    return result


def derive_octet_family():
    """Solve five sourced canonical OPE equations without fixing null directions."""
    coefficients = sp.symbols(" ".join(COEFFICIENT_NAMES))
    c1, c2, c3, c4, c5, c6, c7, d1, d2, d3, d4 = coefficients
    F, F_V, M_0, M_V, N_C = sp.symbols("F F_V M_0 M_V N_C", nonzero=True)
    equations = [
        4*c3+c1,
        c1-c2+c5,
        c5-c6+N_C*M_V/(64*sp.sqrt(2)*sp.pi**2*F_V),
        d1+8*d2-N_C*M_0**2/(64*sp.pi**2*F_V**2)+F**2/(4*F_V**2),
        d3-N_C*M_0**2/(64*sp.pi**2*F_V**2)+F**2/(8*F_V**2),
    ]
    matrix, right_hand_side = sp.linear_eq_to_matrix(equations, coefficients)
    solution = next(iter(sp.linsolve((matrix, right_hand_side), coefficients)))
    residuals = [sp.simplify(equation.subs(dict(zip(coefficients, solution)),
                                          simultaneous=True))
                 for equation in equations]
    return coefficients, equations, matrix, right_hand_side, solution, residuals


def derive_photon_vector_normalization():
    """Match the canonical tensor source to Aloni's Proca mass bilinear.

    The tensor auxiliary-field elimination gives -(P_strength+mixing*F)^2
    divided by 4 M^2. After P=M R and W=R+epsilon A, the strength is
    canonical and the mass term generates the photon-vector bilinear.
    Flavor overlaps q=Tr(TQ), h=Tr(T{chi_plus,Q}) are independent.
    Aloni's matrix is W/sqrt(2). No even quadratic term is supplied by
    Bai's anomalous WZW background functional.
    """
    M, e, g, F_V, z_V, q, h = sp.symbols(
        "M e g F_V z_V q h", nonzero=True)
    tensor, auxiliary_strength, photon_strength = sp.symbols("V_tensor P_strength F_strength")
    source_residue = F_V*q+z_V*h
    mixing = sp.sqrt(2)*e*source_residue
    first_order = (M**2*tensor**2/4-auxiliary_strength*tensor/2
                   -mixing*tensor*photon_strength/2)
    eliminated_tensor = sp.solve(sp.diff(first_order, tensor), tensor)[0]
    reduced = sp.factor(first_order.subs(tensor, eliminated_tensor))
    epsilon = mixing/M
    vector, photon = sp.symbols("W A")
    mass_term = sp.expand(M**2*(vector-epsilon*photon)**2/2)
    rcht_bilinear = mass_term.coeff(vector).coeff(photon)
    aloni_bilinear = -M**2*e*q/(sp.sqrt(2)*g)
    solution = sp.solve(sp.Eq(rcht_bilinear.subs(z_V, 0), aloni_bilinear), F_V)[0]
    return {"eliminated_tensor": eliminated_tensor,
            "reduced_first_order": reduced,
            "rcht_mass_bilinear": rcht_bilinear,
            "aloni_mass_bilinear": aloni_bilinear,
            "aloni_F_V": solution}


def build_artifact():
    """Evaluate the derived affine families at documented numerical inputs.

    Free coefficients are set to zero to select reproducible representatives;
    this choice has no statistical meaning and does not complete source matching.
    Algebraic expressions remain in the derivation output and memory record.
    """
    photon_matching = derive_photon_vector_normalization()
    source_families = derive_source_pole_families()
    source_actions = derive_neutral_source_actions()
    neutral = derive_neutral_pole_matching()
    coefficients = neutral["coefficients"]
    free_names = ["c4", "c5", "c6", "c7", "d4"]
    records = []
    for name, bibtex_key, pion_decay_constant, chiral_decay_constant in (
            ("aloni", "Aloni:2018vki", sp.Rational(92, 1000),
             sp.Rational(92, 1000)),
            ("bai", "Bai:2025fvl", sp.Rational(130, 1000),
             sp.Rational(130, 1000)/sp.sqrt(2))):
        vector_mass = (float(sp.N(sp.sqrt(2)*sp.sqrt(12*sp.pi)*pion_decay_constant, 17))
                       if name == "aloni" else 0.770)
        vector_decay_constant = (float(sp.N(photon_matching["aloni_F_V"].subs({
            sp.Symbol("M", nonzero=True): vector_mass,
            sp.Symbol("g", nonzero=True): sp.sqrt(12*sp.pi)}), 17))
                                 if name == "aloni" else 0.0)
        substitutions = {
            sp.Symbol("g", nonzero=True): sp.sqrt(12*sp.pi),
            sp.Symbol("N_C", nonzero=True): 3,
            sp.Symbol("F", nonzero=True): chiral_decay_constant,
            **({sp.Symbol("F_V", nonzero=True): vector_decay_constant}
               if vector_decay_constant is not None else {}),
            sp.Symbol("M_V", nonzero=True): vector_mass,
            **{coefficient: 0 for coefficient in coefficients
               if str(coefficient) in free_names},
        }
        family = source_families[name+"_common_mass_pole_family"]
        numerical = {}
        for coefficient, value in zip(coefficients, family):
            evaluated = value.subs(substitutions)
            numerical[str(coefficient)] = (None if evaluated.free_symbols
                                           else float(sp.N(evaluated, 17)))
        pole_d3 = numerical["d3"]
        source_d3 = (0.0 if name == "aloni" else float(sp.N(
            source_families["bai_direct_axial_source_d3"].subs(substitutions), 17)))
        record = {
            "label": bibtex_key.lower().replace(":", "_"),
            "scenario": {"rho_input": "common_vector_mass",
                         "kstar_input": "common_vector_mass",
                         "phi_input": "common_vector_mass",
                         "branch": "ideal_flavour_sources",
                         "source_hypothesis": name,
                         "phenomenological_form_factor": 1.0,
                         "F_V_input_source": ("Aloni:2018vki photon-vector mass mixing"
                                              if name == "aloni" else "zero for omitted quadratic mixing operator")},
            "model": "original_11",
            "matching_status": "neutral_pole_amplitude_representative_only",
            "bibtex_keys": [bibtex_key],
            "chi2": None, "success": None, "observations": [],
            "coefficient_units": "dimensionless",
            "canonical_representative_coefficients": numerical,
            "reference_backend_representative_coefficients": {
                coefficient: (-value if value is not None else None)
                for coefficient, value in numerical.items()},
            "numerical_matching_inputs": {
                "g": float(sp.N(sp.sqrt(12*sp.pi), 17)), "N_C": 3,
                "source_f_pi_GeV": float(pion_decay_constant),
                "RChT_F_GeV": float(sp.N(chiral_decay_constant, 17)),
                "common_M_V_GeV": vector_mass, "common_F_V_GeV": vector_decay_constant,
                "neutral_pion_log_U_denominator_GeV": 0.092,
                **({"Fujiwara_a": 2} if name == "aloni" else {})},
            "input_provenance": {
                "g_and_f_pi": (bibtex_key+" supplementary chiral Lagrangian"
                    if name == "aloni" else bibtex_key+" Sec. 2, Eqs. (2.8), (2.12)"),
                "F_conversion": "Distinct rounded source inputs: Aloni f_pi=0.092 GeV; Bai f_pi=0.130 GeV. Their definitions differ by sqrt(2), and both normalize the neutral-pion coefficient in log U with a denominator approximately 0.092 GeV. The shared displayed neutral-pion denominator is rounded to 0.092 GeV; Bai exact 0.130/sqrt(2) rounds to that value and remains unrounded in the coefficient matching. Aloni F uses its independent paper input.",
                "F_V": ("Derived from Aloni:2018vki displayed photon-vector mass bilinear with its vector matrix normalization." if name == "aloni" else "Bai:2025fvl displayed prescription omits the dynamical even-parity photon-vector quadratic mixing operator. Its coefficient F_V is set to zero under the user instruction for omitted operators; this is not a measured physical normalization."),
                "M_V": ("Fujiwara:1984mp relation with a=2 and the independently chosen Aloni f_pi=0.092 GeV; user-specified normalization." if name == "aloni" else "M_V=0.770 GeV repository operator normalization; no paper-defined Bai mass is established here.")},
            "constraint_rank": int(neutral["matrix"].rank()),
            "free_dimension": 11-int(neutral["matrix"].rank()),
            "free_coefficients": free_names,
            "coefficient_null_dimension": len(free_names),
            "coefficient_null_basis_columns": [[
                float(sp.N(sp.diff(value, sp.Symbol(free_name)).subs(substitutions), 17))
                for free_name in free_names] for value in family],
            "representative_policy": "Set the five otherwise-unfixed listed coefficients to zero. These are coordinate representatives of the neutral pole family, with no coefficient errors or prior assigned.",
            "neutral_pole_scope": "Mass-independent neutral pseudoscalar-vector-vector and pseudoscalar-vector-photon amplitudes with independent pseudoscalar momentum and mass-spurion inputs, including the tensor vector-photon exchange contribution.",
            "independent_source_obstruction": {
                "pole_d3": pole_d3, "direct_axial_source_d3": source_d3,
                "difference": pole_d3-source_d3,
                "conclusion": "The neutral pole representative does not reproduce the independent axial-source functional; it is not a complete eleven-operator matching."},
            "verification": {"symbolic_substitution_residuals": [
                float(residual) for residual in neutral["residuals"]]},
        }
        record["source_action_coefficient_GeV_inverse"] = float(sp.N(
            source_actions[name+"_phi_dV_dV"].subs(substitutions), 17))
        # Matching metadata is nested in scenario to preserve the fit schema.
        expected_keys = ["label", "scenario", "model", "control", "parameters",
                         "observables", "fit_statistics", "coefficient_diagnostics"]
        record.update({
            "control": "neutral_pole_matching",
            "archive_source": "scripts/derive_odd_prescription_matching.py",
            "archive_pair_index": len(records),
            "odd_chi2": None, "ancillary_chi2": None,
            "global_minimum_proven": False,
            "fixed_nuisance_measured_rank": int(neutral["matrix"].rank()),
            "fixed_nuisance_coefficient_null_dimension": len(free_names),
            "fixed_nuisance_null_basis_columns": record["coefficient_null_basis_columns"],
            "projection_relative_defect": 0.0,
            "null_annihilation_norm": 0.0,
            "even_parameters": [vector_mass, vector_mass, 0.0, 0.0,
                                vector_decay_constant, 0.0],
            "even_parameters_by_name": {
                "m_rho_GeV": vector_mass, "m_Kstar_GeV": vector_mass,
                "Delta_V_GeV2": 0.0, "kappa_V": 0.0,
                "F_V_GeV": vector_decay_constant, "z_V_GeV_inverse": 0.0},
            "shared_total_widths_gev": {},
            "local_prediction_uncertainty_scope": "No fit or prediction uncertainty is assigned. The rank and null diagnostics describe exact neutral pole matching constraints rather than experimental measurements.",
        })
        record['parameters'] = dict(record['even_parameters_by_name'])
        record['parameters'].update(record['canonical_representative_coefficients'])
        record['parameters']['omega_phi_mixing_angle_rad'] = 0.0
        record['fit_statistics'] = {'chi2': None, 'odd_chi2': None,
                                    'ancillary_chi2': None, 'success': None}
        record['observables'] = []
        record['coefficient_diagnostics'] = dict(
            measured_rank=record['fixed_nuisance_measured_rank'],
            null_dimension=record['fixed_nuisance_coefficient_null_dimension'],
            null_basis_columns=record['fixed_nuisance_null_basis_columns'],
            projection_relative_defect=0.0, null_annihilation_norm=0.0)
        record["scenario"]["even_sector_policy"] = "Common vector mass; mass-spurion and double-trace source coefficients are explicitly zero because these prescriptions add no such terms. F_V is derived from the Aloni quadratic mixing term; Bai F_V is zero for the omitted quadratic mixing operator, not a measured physical value."
        record["scenario"]["matching_metadata"] = {
            key: record.pop(key) for key in list(record) if key not in expected_keys}
        records.append({key: record[key] for key in expected_keys})
    inputs = [
        "scripts/derive_odd_prescription_matching.py",
        "memory/rcht/odd-parity/odd-prescription-matching.md",
        "memory/extracts/Aloni:2018vki/vector-meson-dominance.md",
        "memory/extracts/Bai:2025fvl/no-axial-vector-and-rcht-matching.md",
    ]
    return {
        "scope": "Numerical neutral pole representatives of the Aloni and Bai prescriptions in the provisional fit coefficient dictionaries; no fit was performed.",
        "selection_policy": "Distinct rounded source decay constants; Aloni mass from the user-specified Fujiwara relation; Bai mass from repository normalization; Aloni F_V from its quadratic mixing term and Bai F_V zero for its omitted quadratic mixing operator; free coefficients set to zero.",
        "conventions": "epsilon^{0123}=+1; metric (+---); canonical c_i=-c_i^R and d_i=-d_i^R; dimensionless odd coefficients.",
        "uncertainty_policy": "Numerical representatives carry no fit quality or uncertainties. Five unconstrained directions remain; independent-source matching is obstructed.",
        "results": records,
        "sha256": {relative: hashlib.sha256((ROOT/relative).read_bytes()).hexdigest()
                   for relative in inputs},
    }


def derive_aloni_vector_two_pion_matching():
    """Match the displayed rho two-pion vertex on its vector pole.

    The paper matrices are Phi/sqrt(2) and V_ideal/sqrt(2). The
    positive G_V convention is selected after matching squared couplings.
    This pole match does not equate off-shell tensor and Proca vertices.
    """
    pion_plus_generator = sp.zeros(3)
    pion_plus_generator[0, 1] = 1
    pion_minus_generator = pion_plus_generator.T
    rho_generator = sp.diag(1, -1, 0)/sp.sqrt(2)
    flavor_trace = sp.trace(rho_generator*(
        pion_plus_generator*pion_minus_generator
        -pion_minus_generator*pion_plus_generator))
    vector_coupling, chiral_decay_constant, vector_mass = sp.symbols(
        "g F M_V", positive=True)
    tensor_coupling = sp.Symbol("G_V", positive=True)
    source_coupling = sp.simplify(vector_coupling*flavor_trace/(2*sp.sqrt(2)))
    tensor_pole_coupling = tensor_coupling*vector_mass/chiral_decay_constant**2
    equation = sp.Eq(tensor_pole_coupling**2, source_coupling**2)
    solution = sp.solve(equation, tensor_coupling)[0]
    return {"rho_pion_flavor_trace": flavor_trace,
            "source_rho_pion_coupling": source_coupling,
            "squared_coupling_equation": equation,
            "positive_G_V": solution}


def build_simplified_artifacts():
    """Return one flat parameter document per literature prescription.

    Widths omitted by the prescription are explicit zeroes, as documented
    in metadata; they are not estimates of measured total widths.
    """
    artifacts = {}
    for record in build_artifact()["results"]:
        metadata = {
            "bibtex_references": record["scenario"]["matching_metadata"]["bibtex_keys"],
        }
        parameters = dict(record["parameters"])
        if record["label"] == "aloni_2018vki":
            matching_inputs = record["scenario"]["matching_metadata"]["numerical_matching_inputs"]
            parameters["G_V_GeV"] = float(derive_aloni_vector_two_pion_matching()[
                "positive_G_V"].subs({
                    sp.Symbol("g", positive=True): sp.sqrt(12*sp.pi),
                    sp.Symbol("F", positive=True): sp.Rational(92, 1000),
                    sp.Symbol("M_V", positive=True): matching_inputs["common_M_V_GeV"]}))
        else:
            parameters["G_V_GeV"] = 0.0
        for vector_name in ("Kstar0", "Kstarplus", "etap", "omega", "phi", "rho0", "rhoplus"):
            parameters[f"Gamma_total_{vector_name}_GeV"] = 0.0
        artifacts[record["label"]] = {"metadata": metadata, "parameters": parameters}
    return artifacts


def main():
    """Write the deterministic JSON and print the actual SymPy solution steps."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--derive-only", action="store_true", help="Print algebra without writing JSON.")
    arguments = parser.parse_args()
    for label, calculation in (("Vector two-pion normalization", derive_aloni_vector_two_pion_matching()),
                               ("Photon vector normalization", derive_photon_vector_normalization()),
                               ("Neutral pole", derive_neutral_pole_matching()),
                               ("Complete double VMD TFF", derive_complete_vmd_tff()),
                               ("Neutral source actions", derive_neutral_source_actions()),
                               ("Component normalization", derive_source_component_normalization()),
                               ("Full ALP coefficients", derive_full_alp_vector_coefficients()),
                               ("Source pole families", derive_source_pole_families())):
        for name, value in calculation.items():
            if isinstance(value, (sp.Basic, sp.MatrixBase)):
                print(f"{label} {name}: {sp.latex(value)}")
            elif isinstance(value, (list, tuple)):
                print(f"{label} {name}: {sp.latex(sp.Tuple(*value))}")
    _, equations, matrix, right_hand_side, solution, _ = derive_octet_family()
    for index, equation in enumerate(equations, 1):
        print(f"Input OPE equation {index}: {sp.latex(sp.Eq(equation, 0))}")
    print(f"Coefficient matrix: {sp.latex(matrix)}")
    print(f"Right-hand side: {sp.latex(right_hand_side)}")
    print(f"Solved family: {sp.latex(sp.Tuple(*solution))}")
    if arguments.derive_only:
        return
    for label, artifact in build_simplified_artifacts().items():
        filename = f"odd_prescription_matching_{label}.json"
        for directory in ("data/vector_sector", "alpaca_chiral"):
            path = ROOT/directory/filename
            path.write_text(json.dumps(artifact, indent=2)+"\n")
            print(path.relative_to(ROOT))


if __name__ == "__main__":
    main()
