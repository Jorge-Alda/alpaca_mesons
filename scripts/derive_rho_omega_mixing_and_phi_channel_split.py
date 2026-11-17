#!/usr/bin/env python3
"""rho0-omega mixing and the phi -> K+K- versus K0 K0bar split at first order in m_u - m_d.

Invocation:
    .venv/bin/python scripts/derive_rho_omega_mixing_and_phi_channel_split.py
Tests:
    .venv/bin/python -m unittest scripts.test_derive_rho_omega_mixing_and_phi_channel_split

Derivation only (task 5.1 of extend-joint-vector-model-diagnostics). No fit,
refit, scan or optimizer is run; numbers are evaluated at the fixed accepted
endpoint `basin_probes/even_shared_widths` of the joint vector likelihood.

Conventions (Alda:2026nrd, Ecker:1988te, and the memory records
rcht/foundations/vector-nonideal-mass-action.md, vector-mass-extension.md,
rcht-photon-vector-diagonalization.md):
  * natural units, GeV; e_m^V and kappa_V dimensionless, Delta_V and the
    spurion entries x_q in GeV^2, F_V in GeV, z_V and y_V in GeV^-1;
  * vacuum spurion chi_+ = diag(x_u, x_d, x_s), x_q = 4 B_0 m_q, with the
    Alda:2026nrd arXiv-v1 central values; the joint fit uses x_u = x_d = x_l;
  * neutral generators T_rho = diag(1,-1,0)/sqrt2, T_omegaI = diag(1,1,0)/sqrt2,
    T_phiI = diag(0,0,1); physical T_omega = cos(theta_V) T_omegaI - sin(theta_V) T_phiI,
    T_phi = sin(theta_V) T_omegaI + cos(theta_V) T_phiI;
  * mass Lagrangian (order e^0, linear in chi_+):
        M_V^2/4 <V V> - e_m^V <V V chi_+> + Delta_V/12 <V><V> + kappa_V/4 <V><V chi_+>;
    for real neutral components V = sum_i v_i T_i the mass matrix is
    2 d^2 L / dv_i dv_j; for complex charged fields it is 2 d^2 L / dV dV^dagger;
  * photon source matrix S = F_V Q + z_V {chi_+, Q} + y_V <chi_+ Q> 1, residues
    R_i = <T_i S>;
  * propagator matrix (M^2 - s 1)^(-1) in the transverse spin-one sector, with
    constant on-shell widths m^2 -> m^2 - i m Gamma (assumption).

Order: tree level, first order in x_u - x_d, order e^0 for the mass terms and
order e for the photon couplings, as in the joint fit. The order-e^2
photon-exchange mixing is evaluated only as a size estimate and is not part of
the result. No new coefficient is introduced.

Steps (each checked by assertions or by the unit tests):
  1. Neutral 3x3 and charged mass matrices with x_u != x_d; isospin-limit
     reduction to the records; the identity
        (M^2)_{rho omegaI} - kappa_V (x_u - x_d)/2 = m^2_{K*+} - m^2_{K*0}.
  2. First-order pole residue of a two-state propagator: the physical state V
     carries a component (M^2)_{rho V} / (mt_V^2 - mt_rho^2) of rho0.
  3. Isospin-breaking parts of the photon residues (z_V and y_V terms).
  4. phi -> K Kbar flavor traces with the rho0 component of phi: first-order
     charged/neutral split, sum unchanged.
  5. Numbers at the endpoint: admixtures, real-diagonalization failure for
     rho-omega, omega -> pi+ pi- check against PDG, leptonic and radiative
     effects, phi split, new rotated phi residuals, and the decomposition of
     the task-1.2 residual into sum and ratio directions.

Output summary (2026-11-17 run, endpoint basin_probes/even_shared_widths):
  * (M^2)_{rho omegaI} = (x_u - x_d)(kappa_V/2 - 2 e_m^V) = -4.865e-3 GeV^2,
    (M^2)_{rho phiI} = kappa_V (x_u - x_d)/(2 sqrt2) = -3.38e-4 GeV^2;
    physical-basis entries -4.810e-3 (rho-omega) and -7.99e-4 (rho-phi) GeV^2.
  * rho0 component of omega: -0.0047 + 0.0443 i (modulus 0.0445); a real
    zero-width diagonalization would give -0.42, outside first order.
  * omega -> pi+ pi- branching fraction 0.028 (fit G_V) or 0.036 (measured
    rho0 width), PDG 0.0153 +- 0.0012; K*+ - K*0 mass 2.46 MeV lighter
    (measured 3.88 MeV, electromagnetic part not included).
  * Gamma(omega -> e+e-) +15% at fixed coordinates (z_V source term +16%,
    mixing -0.9%), Gamma(rho -> e+e-) +1.8%; radiative rows <= 0.3 sigma.
  * phi -> K+K-/K0K0bar ratio +0.52% (preferred about -4.4%); phi block
    chi^2 17.88 -> 19.62 at the fixed endpoint; the order-e^2 photon-exchange
    estimate changes the sign (-0.43%).
  * task-1.2 residual: a common normalization removes 12.95 of 17.88, a
    charged/neutral split alone 9.64; 4.92 remains in the split direction
    after the normalization (p = 0.026 for one, 0.085 for two degrees of freedom).
"""
import json
import math
import sys
from pathlib import Path

import numpy as np
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import diagnose_joint_phi_kaon_residuals as phi_diagnostic  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
ENDPOINT = 'basin_probes/even_shared_widths'

x_u, x_d, x_s, x_l = sp.symbols('x_u x_d x_s x_l', real=True)
vector_mass_squared, mass_splitting_coefficient = sp.symbols('M_V2 e_m_V', real=True)
singlet_shift, double_trace_mass_coefficient = sp.symbols('Delta_V kappa_V', real=True)
F_V, z_V, y_V = sp.symbols('F_V z_V y_V', real=True)
mixing_angle = sp.symbols('theta_V', real=True)

CHARGE_MATRIX = sp.diag(sp.Rational(2, 3), -sp.Rational(1, 3), -sp.Rational(1, 3))
GENERATOR_RHO = sp.diag(1, -1, 0)/sp.sqrt(2)
GENERATOR_OMEGA_IDEAL = sp.diag(1, 1, 0)/sp.sqrt(2)
GENERATOR_PHI_IDEAL = sp.diag(0, 0, 1)


def matrix_unit(row, column):
    """Flavor matrix with one unit entry."""
    result = sp.zeros(3)
    result[row, column] = 1
    return result


def spurion(isospin_symmetric=False):
    """Vacuum chi_+; x_u = x_d = x_l when isospin_symmetric."""
    return sp.diag(x_l, x_l, x_s) if isospin_symmetric else sp.diag(x_u, x_d, x_s)


def mass_lagrangian(vector, chi):
    """Order-e^0 mass terms of the joint fit, one Lorentz contraction suppressed."""
    return sp.expand(vector_mass_squared/4*(vector*vector).trace()
                     - mass_splitting_coefficient*(vector*vector*chi).trace()
                     + singlet_shift/12*vector.trace()**2
                     + double_trace_mass_coefficient/4*vector.trace()*(vector*chi).trace())


def neutral_mass_matrix(chi):
    """Mass matrix in the (rho0, omegaI, phiI) basis: 2 d^2 L / dv_i dv_j."""
    components = sp.symbols('v_rho v_omega v_phi', real=True)
    generators = (GENERATOR_RHO, GENERATOR_OMEGA_IDEAL, GENERATOR_PHI_IDEAL)
    vector = sum((c*g for c, g in zip(components, generators)), sp.zeros(3))
    lagrangian = mass_lagrangian(vector, chi)
    return sp.Matrix(3, 3, lambda i, j: sp.simplify(
        2*sp.diff(lagrangian, components[i], components[j])))


def charged_mass_squared(row, column, chi):
    """m^2 = 2 d^2 L / dV dV^dagger for the complex field in entry (row, column)."""
    field, conjugate = sp.symbols('V Vbar')
    vector = field*matrix_unit(row, column) + conjugate*matrix_unit(column, row)
    return sp.simplify(2*sp.diff(mass_lagrangian(vector, chi), field, conjugate))


def source_matrix(chi):
    """Photon source matrix S_gamma of the photon-vector diagonalization record."""
    return F_V*CHARGE_MATRIX + z_V*(chi*CHARGE_MATRIX + CHARGE_MATRIX*chi) \
        + y_V*(chi*CHARGE_MATRIX).trace()*sp.eye(3)


def photon_residues(chi):
    """R_i = <T_i S_gamma> for rho0, omegaI, phiI."""
    source = source_matrix(chi)
    return {name: sp.simplify((generator*source).trace()) for name, generator in
            (('rho', GENERATOR_RHO), ('omega_I', GENERATOR_OMEGA_IDEAL),
             ('phi_I', GENERATOR_PHI_IDEAL))}


def first_order_pole_component():
    """Residue of a two-state propagator at the pole of state 'b' to first order in mixing.

    With inverse propagator [[a - s, delta], [delta, b - s]] (complex a, b allowed),
    the residue matrix at s = b factorizes as w w^T with w = (component, 1); the
    function returns the a-component at first order in delta.
    """
    a, b, delta, s, small = sp.symbols('a b delta s small')
    propagator = sp.Matrix([[a - s, small*delta], [small*delta, b - s]]).inv()
    first_order = sp.simplify(sp.diff(propagator[0, 1], small).subs(small, 0))
    zeroth_order = sp.simplify(propagator[1, 1].subs(small, 0))
    component = sp.simplify(sp.limit((b - s)*first_order, s, b)
                            / sp.limit((b - s)*zeroth_order, s, b))
    assert sp.simplify(component - delta/(b - a)) == 0
    return component, (a, b, delta)


def phi_kaon_traces(component_rho):
    """Traces tr(T [K, Kbar]) for phi -> K+K- and K0 K0bar with a rho0 component.

    T_phi + component_rho * T_rho; the fit squares (T_qq - T_ss) for q = u, d.
    """
    phi = sp.sin(mixing_angle)*GENERATOR_OMEGA_IDEAL + sp.cos(mixing_angle)*GENERATOR_PHI_IDEAL
    generator = phi + component_rho*GENERATOR_RHO
    kaon_plus, kaon_minus = matrix_unit(0, 2), matrix_unit(2, 0)
    kaon_zero, antikaon_zero = matrix_unit(1, 2), matrix_unit(2, 1)
    charged = (generator*(kaon_plus*kaon_minus - kaon_minus*kaon_plus)).trace()
    neutral = (generator*(kaon_zero*antikaon_zero - antikaon_zero*kaon_zero)).trace()
    return sp.simplify(charged), sp.simplify(neutral)


def derive():
    """Run the symbolic steps 1-4 and return all intermediate expressions."""
    chi, chi_symmetric = spurion(), spurion(isospin_symmetric=True)
    neutral = neutral_mass_matrix(chi)
    neutral_symmetric = neutral_mass_matrix(chi_symmetric)
    record_A = vector_mass_squared - 4*mass_splitting_coefficient*x_l + 2*singlet_shift/3 \
        + 2*double_trace_mass_coefficient*x_l
    record_B = sp.sqrt(2)*singlet_shift/3 + double_trace_mass_coefficient*(x_l + x_s)/sp.sqrt(2)
    record_D = vector_mass_squared - 4*mass_splitting_coefficient*x_s + singlet_shift/3 \
        + double_trace_mass_coefficient*x_s
    isospin_limit_differences = [
        sp.simplify(neutral_symmetric[1, 1] - record_A),
        sp.simplify(neutral_symmetric[1, 2] - record_B),
        sp.simplify(neutral_symmetric[2, 2] - record_D),
        sp.simplify(neutral_symmetric[0, 1]), sp.simplify(neutral_symmetric[0, 2]),
        sp.simplify(neutral_symmetric[0, 0] - (vector_mass_squared - 4*mass_splitting_coefficient*x_l))]
    assert all(value == 0 for value in isospin_limit_differences)
    charged = {'rho_plus': charged_mass_squared(0, 1, chi),
               'kstar_plus': charged_mass_squared(0, 2, chi),
               'kstar_zero': charged_mass_squared(1, 2, chi)}
    # Off-diagonal entries linear in x_u - x_d (exact for this Lagrangian).
    rho_omega_ideal, rho_phi_ideal = neutral[0, 1], neutral[0, 2]
    expected_rho_omega = (x_u - x_d)*(-2*mass_splitting_coefficient + double_trace_mass_coefficient/2)
    expected_rho_phi = double_trace_mass_coefficient*(x_u - x_d)/(2*sp.sqrt(2))
    assert sp.simplify(rho_omega_ideal - expected_rho_omega) == 0
    assert sp.simplify(rho_phi_ideal - expected_rho_phi) == 0
    kstar_identity = sp.simplify(rho_omega_ideal - double_trace_mass_coefficient*(x_u - x_d)/2
                                 - (charged['kstar_plus'] - charged['kstar_zero']))
    assert kstar_identity == 0
    # rho0 and rho+ diagonal masses stay equal at this order; isospin-average diagonals unchanged.
    rho_degeneracy = sp.simplify(neutral[0, 0] - charged['rho_plus'])
    diagonal_shifts = [sp.simplify((neutral[i, i] - neutral_symmetric[i, i]).subs(x_l, (x_u + x_d)/2))
                       for i in range(3)]
    assert rho_degeneracy == 0 and all(value == 0 for value in diagonal_shifts)
    # Rotation to physical omega, phi.
    rho_omega_physical = sp.cos(mixing_angle)*rho_omega_ideal - sp.sin(mixing_angle)*rho_phi_ideal
    rho_phi_physical = sp.sin(mixing_angle)*rho_omega_ideal + sp.cos(mixing_angle)*rho_phi_ideal
    component, _ = first_order_pole_component()
    residues = photon_residues(chi)
    residues_symmetric = photon_residues(chi_symmetric)
    residue_isospin_breaking = {name: sp.simplify(
        (residues[name] - residues_symmetric[name]).subs(x_l, (x_u + x_d)/2)) for name in residues}
    expected_residue_breaking = {'rho': z_V*(x_u - x_d)/(3*sp.sqrt(2)),
                                 'omega_I': (z_V + y_V)*(x_u - x_d)/sp.sqrt(2),
                                 'phi_I': y_V*(x_u - x_d)/2}
    for name, value in expected_residue_breaking.items():
        assert sp.simplify(residue_isospin_breaking[name] - value) == 0
    # phi -> K Kbar traces with a rho0 component c.
    component_rho = sp.symbols('c', real=True)
    charged_trace, neutral_trace = phi_kaon_traces(component_rho)
    zero_charged, zero_neutral = phi_kaon_traces(0)
    assert sp.simplify(zero_charged - zero_neutral) == 0
    ratio_first_order = sp.simplify(sp.diff(charged_trace**2/neutral_trace**2, component_rho).subs(component_rho, 0))
    isospin_limit_trace = sp.sin(mixing_angle)/sp.sqrt(2) - sp.cos(mixing_angle)
    assert sp.simplify(zero_charged - isospin_limit_trace) == 0
    assert sp.simplify(ratio_first_order - 2*sp.sqrt(2)/isospin_limit_trace) == 0
    ratio_first_order = 2*sp.sqrt(2)*component_rho/isospin_limit_trace
    # Calibration of e_m^V from the rho and charged K* masses when x_u != x_d.
    charged_calibration = sp.solve(sp.Eq(charged['rho_plus'] - charged['kstar_plus'],
                                         sp.Symbol('rho_minus_kstar_plus')), mass_splitting_coefficient)[0]
    assert sp.simplify(charged_calibration - sp.Symbol('rho_minus_kstar_plus')/(2*(x_s - x_d))) == 0
    sum_first_order = sp.simplify(sp.diff(charged_trace**2 + neutral_trace**2, component_rho).subs(component_rho, 0))
    assert sum_first_order == 0
    return dict(neutral_mass_matrix=neutral, isospin_limit_differences=isospin_limit_differences,
                charged_mass_squared=charged, rho_omega_ideal=rho_omega_ideal,
                rho_phi_ideal=rho_phi_ideal, kstar_identity=kstar_identity,
                rho_degeneracy=rho_degeneracy, diagonal_shifts=diagonal_shifts,
                rho_omega_physical=rho_omega_physical, rho_phi_physical=rho_phi_physical,
                pole_component=component, residues=residues,
                residue_isospin_breaking=residue_isospin_breaking,
                phi_traces=(charged_trace, neutral_trace), phi_ratio_first_order=ratio_first_order,
                phi_sum_first_order=sum_first_order, component_symbol=component_rho,
                charged_calibration=charged_calibration)


def endpoint_inputs():
    """Read the accepted endpoint and the measured inputs used by the joint fit."""
    import fit_babar_kkpi_global_vector_sector as joint_fit
    ledger = json.loads((ROOT/'data/vector_sector/joint_vector_fit_adequacy.json').read_text())
    endpoint = np.array(next(entry for entry in ledger['endpoints'] if entry['id'] == ENDPOINT)['endpoint'])
    preliminary, _ = joint_fit.load_inputs()
    even = np.array([preliminary[key] for key in joint_fit.EVEN_NAMES]) + joint_fit.EVEN_SCALES*endpoint[:6]
    alda = json.loads((ROOT/'data/vector_sector/alda_2026nrd_arxiv_v1.json').read_text())
    assumptions = json.loads((ROOT/'data/vector_sector/stage3_results.json').read_text())['assumptions']
    odd = json.loads((ROOT/'data/odd_parity/odd_low_energy_measurements.json').read_text())
    pdg = json.loads((ROOT/'data/vector_sector/pdg2024.json').read_text())
    sector, _ = joint_fit.calibration.J.action_sector(even, 'single_trace_negative', assumptions)
    totals = {name: odd['total_widths'][name]['value_gev'] + odd['total_widths'][name]['sigma_gev']*endpoint[7 + index]
              for index, name in enumerate(joint_fit.FREE_TOTALS)}
    branching = {row['channel']: row for row in odd['observations'] if row.get('status') != 'deferred'}
    spurion_values = 4*np.array(alda['B0mq_gev2'])
    assert abs((spurion_values[0] + spurion_values[1])/2 - assumptions['xl_gev2']) < 1e-12
    assert abs(spurion_values[2] - assumptions['xs_gev2']) < 1e-12
    return dict(rho_mass=even[0], kstar_mass=even[1], singlet_shift=even[2], kappa=even[3],
                F_V=even[4], z_V=even[5], G_V=endpoint[6], sector=sector, totals=totals,
                spurion=spurion_values, alpha=assumptions['alpha_em'],
                physical=alda['physical_masses_gev'], branching=branching, pdg=pdg,
                decay_constant=joint_fit.F)


def two_body_momentum(parent, first, second):
    """Decay momentum in GeV."""
    return math.sqrt((parent**2 - (first + second)**2)*(parent**2 - (first - second)**2))/(2*parent)


def radiative_momentum(parent, daughter):
    """Photon momentum of parent -> daughter gamma in GeV."""
    return (parent**2 - daughter**2)/(2*parent)


def projection_remainder(difference, weight, direction):
    """chi^2 left after freeing one linear direction: d^T W d - (v^T W d)^2 / v^T W v."""
    return float(difference @ weight @ difference
                 - (direction @ weight @ difference)**2/(direction @ weight @ direction))


def evaluate(result=None, inputs=None):
    """Step 5: numbers at the accepted endpoint; no fit is performed."""
    result = derive() if result is None else result
    inputs = endpoint_inputs() if inputs is None else inputs
    xu, xd, xs = inputs['spurion']
    rho_mass, kstar_mass = inputs['rho_mass'], inputs['kstar_mass']
    em = (rho_mass**2 - kstar_mass**2)/(2*(xs - (xu + xd)/2))
    rho_generator, omega_generator, phi_generator = (np.real(g) for g in inputs['sector'].generators)
    sin_theta = phi_generator[0, 0]*math.sqrt(2)
    theta = math.atan2(sin_theta, phi_generator[2, 2])
    substitution = {x_u: xu, x_d: xd, x_s: xs, x_l: (xu + xd)/2,
                    vector_mass_squared: rho_mass**2 + 4*em*(xu + xd)/2,
                    mass_splitting_coefficient: em, singlet_shift: inputs['singlet_shift'],
                    double_trace_mass_coefficient: inputs['kappa'], F_V: inputs['F_V'],
                    z_V: inputs['z_V'], y_V: 0., mixing_angle: theta}
    numeric = lambda expression: float(sp.N(sp.sympify(expression).subs(substitution)))
    # Cross-check: the isospin-symmetric neutral matrix reproduces the fit's masses and angle.
    neutral_symmetric = np.array(sp.N(neutral_mass_matrix(spurion(True)).subs(substitution)), dtype=float)[1:, 1:]
    eigenvalues = np.linalg.eigvalsh(neutral_symmetric)
    assert np.allclose(np.sqrt(eigenvalues), inputs['sector'].masses[1:], atol=1e-12)
    assert abs(0.5*math.atan2(2*neutral_symmetric[0, 1], neutral_symmetric[1, 1] - neutral_symmetric[0, 0]) - theta) < 1e-12
    masses = dict(zip(('rho', 'omega', 'phi'), inputs['sector'].masses))
    totals = inputs['totals']
    complex_mass = {name: masses[name]**2 - 1j*masses[name]*totals[key]
                    for name, key in (('rho', 'rho0'), ('omega', 'omega'), ('phi', 'phi'))}
    rho_omega = numeric(result['rho_omega_physical'])
    rho_phi = numeric(result['rho_phi_physical'])
    rho_omega_ideal = numeric(result['rho_omega_ideal'])
    rho_phi_ideal = numeric(result['rho_phi_ideal'])
    kstar_split = numeric(result['charged_mass_squared']['kstar_plus'] - result['charged_mass_squared']['kstar_zero'])
    # rho0 component of the physical omega and phi; omega component of rho is the negative.
    rho_in_omega = rho_omega/(complex_mass['omega'] - complex_mass['rho'])
    rho_in_phi = rho_phi/(complex_mass['phi'] - complex_mass['rho'])
    rho_in_omega_zero_width = rho_omega/(masses['omega']**2 - masses['rho']**2)
    rho_in_phi_zero_width = rho_phi/(masses['phi']**2 - masses['rho']**2)
    real_diagonalization_angle = 0.5*math.atan2(2*rho_omega, masses['omega']**2 - masses['rho']**2)
    # Order-e^2 photon-exchange entry s * 2 e^2 R_i R_j / (M_i M_j) (size estimate only).
    residues = dict(zip(('rho', 'omega', 'phi'), inputs['sector'].residues))
    electric_charge_squared = 4*math.pi*inputs['alpha']
    photon_rho_omega = masses['omega']**2*2*electric_charge_squared*residues['rho']*residues['omega']/(masses['rho']*masses['omega'])
    photon_rho_phi = masses['phi']**2*2*electric_charge_squared*residues['rho']*residues['phi']/(masses['rho']*masses['phi'])
    # omega -> pi+ pi- through the rho0 component (no direct coupling at this order).
    pion = inputs['physical']['pi_charged']
    momentum_rho = two_body_momentum(masses['rho'], pion, pion)
    momentum_omega = two_body_momentum(masses['omega'], pion, pion)
    tree_rho_width = inputs['G_V']**2*2*masses['rho']**2*momentum_rho**3/(
        12*math.pi*inputs['decay_constant']**4*masses['rho']**2)
    width_scaling = (masses['omega']/masses['rho'])**2*(momentum_omega/momentum_rho)**3
    omega_pipi_branching = {
        'with_fit_G_V': abs(rho_in_omega)**2*tree_rho_width*width_scaling/totals['omega'],
        'with_measured_rho0_width': abs(rho_in_omega)**2*totals['rho0']*width_scaling/totals['omega']}
    with_photon = (rho_omega + photon_rho_omega)/(complex_mass['omega'] - complex_mass['rho'])
    omega_pipi_branching['with_fit_G_V_and_photon_exchange_estimate'] = (
        abs(with_photon)**2*tree_rho_width*width_scaling/totals['omega'])
    # Photon residues: isospin-breaking parts (z_V) and rho0 admixture at the poles.
    breaking = {name: numeric(value) for name, value in result['residue_isospin_breaking'].items()}
    cos_theta = math.cos(theta)
    omega_breaking = cos_theta*breaking['omega_I'] - math.sin(theta)*breaking['phi_I']
    omega_pole_residue = residues['omega'] + omega_breaking + rho_in_omega*residues['rho']
    rho_pole_residue = residues['rho'] + breaking['rho'] - rho_in_omega*residues['omega']
    leptonic_relative_change = {
        'omega': abs(omega_pole_residue)**2/residues['omega']**2 - 1,
        'omega_source_term_only': (1 + omega_breaking/residues['omega'])**2 - 1,
        'omega_mixing_only': abs(1 + rho_in_omega*residues['rho']/residues['omega'])**2 - 1,
        'rho': abs(rho_pole_residue)**2/residues['rho']**2 - 1}
    leptonic_sigma = {name: inputs['pdg']['leptonic_widths'][name]['sigma_kev']/inputs['pdg']['leptonic_widths'][name]['value_kev']
                      for name in ('rho', 'omega')}
    # Radiative rows: bound 2 r |Re c| + r^2 |c|^2 with r = |A_other/A_own| from measured widths.
    def measured_width(label, total):
        return inputs['branching'][label]['value']*totals[total]
    radiative = {}
    for own, other, own_total, other_total, own_mass, other_mass, daughter in (
            ('rho0_pi0_gamma', 'omega_pi0_gamma', 'rho0', 'omega', masses['rho'], masses['omega'], 'pi0'),
            ('omega_pi0_gamma', 'rho0_pi0_gamma', 'omega', 'rho0', masses['omega'], masses['rho'], 'pi0'),
            ('rho0_eta_gamma', 'omega_eta_gamma', 'rho0', 'omega', masses['rho'], masses['omega'], 'eta'),
            ('omega_eta_gamma', 'rho0_eta_gamma', 'omega', 'rho0', masses['omega'], masses['rho'], 'eta')):
        daughter_mass = inputs['physical'][daughter]
        amplitude_ratio = math.sqrt(measured_width(other, other_total)/measured_width(own, own_total)
                                    * radiative_momentum(own_mass, daughter_mass)**3
                                    / radiative_momentum(other_mass, daughter_mass)**3)
        bound = 2*amplitude_ratio*abs(rho_in_omega.real) + amplitude_ratio**2*abs(rho_in_omega)**2
        relative_sigma = inputs['branching'][own]['sigma']/inputs['branching'][own]['value']
        radiative[own] = dict(amplitude_ratio=amplitude_ratio, relative_rate_bound=bound,
                              bound_in_row_sigma=bound/relative_sigma)
    # BaBar ground-state propagators at the lowest bin centre (estimate).
    lowest_energy = 1.37
    babar_relative_offdiagonal = abs(rho_omega/(complex_mass['rho'] - lowest_energy**2))
    # phi -> K Kbar split from the rho0 component of phi.
    charged_trace, neutral_trace = result['phi_traces']
    trace_value = lambda expression, value: complex(sp.N(expression.subs({mixing_angle: theta, result['component_symbol']: 0}))) \
        + value*complex(sp.N(sp.diff(expression, result['component_symbol']).subs(mixing_angle, theta)))
    charged_factor = abs(trace_value(charged_trace, rho_in_phi))**2/abs(trace_value(charged_trace, 0))**2
    neutral_factor = abs(trace_value(neutral_trace, rho_in_phi))**2/abs(trace_value(neutral_trace, 0))**2
    phi_ratio_change = charged_factor/neutral_factor - 1
    linear_ratio = lambda component: numeric(result['phi_ratio_first_order'].subs(result['component_symbol'], component))
    rho_in_phi_with_photon = (rho_phi + photon_rho_phi)/(complex_mass['phi'] - complex_mass['rho'])
    phi_ratio_change_linear = linear_ratio(rho_in_phi.real)
    phi_ratio_change_with_photon_estimate = linear_ratio(rho_in_phi_with_photon.real)
    calibration_change = (xs - (xu + xd)/2)/(xs - xd) - 1
    diagnostic = phi_diagnostic.diagnostic()
    covariance = np.array(diagnostic['covariance'])
    weight = np.linalg.inv(covariance)
    observed = np.array(diagnostic['observed_branching_fractions'])
    predicted = np.array(diagnostic['predicted_branching_fractions'])
    shifted = predicted*np.array([charged_factor, neutral_factor])
    cholesky = np.linalg.cholesky(covariance)
    shifted_cholesky = np.linalg.solve(cholesky, shifted - observed)
    _, shifted_symmetric, _ = phi_diagnostic.channel_residuals(shifted_cholesky, covariance)
    difference = predicted - observed
    sum_direction = predicted
    split_direction = np.array([predicted[0], -predicted[1]])
    decomposition = dict(
        block_chi2=float(difference @ weight @ difference),
        remaining_after_common_normalization=projection_remainder(difference, weight, sum_direction),
        remaining_after_charged_neutral_split=projection_remainder(difference, weight, split_direction))
    decomposition['removed_by_common_normalization'] = decomposition['block_chi2'] - decomposition['remaining_after_common_normalization']
    decomposition['removed_by_split'] = decomposition['block_chi2'] - decomposition['remaining_after_charged_neutral_split']
    from scipy.stats import chi2 as chi_square
    decomposition['remaining_after_common_normalization_p_one_dof'] = float(
        chi_square.sf(decomposition['remaining_after_common_normalization'], 1))
    decomposition['remaining_after_common_normalization_p_two_rows'] = float(
        chi_square.sf(decomposition['remaining_after_common_normalization'], 2))
    split_coefficient = -(split_direction @ weight @ difference)/(split_direction @ weight @ split_direction)
    # Comparison: VPP isospin-breaking option, phi coupling split (r2 - r3)(x_u - x_d) (operators record),
    # width-ratio change about 2 (r2 - r3)(x_u - x_d) at the task-3.2 best fit r2 = 1.59, r3 = 0.35 GeV^-2.
    vpp_ratio_change_task_3_2 = 2*(1.59 - 0.35)*(xu - xd)
    return dict(
        endpoint=ENDPOINT, mixing_angle_rad=theta, mass_splitting_coefficient=em,
        spurion_gev2=[xu, xd, xs],
        rho_omega_ideal_gev2=rho_omega_ideal, rho_phi_ideal_gev2=rho_phi_ideal,
        rho_omega_physical_gev2=rho_omega, rho_phi_physical_gev2=rho_phi,
        kstar_plus_minus_zero_mass_squared_gev2=kstar_split,
        kstar_plus_minus_zero_mass_mev=1000*kstar_split/(2*kstar_mass),
        measured_kstar_plus_minus_zero_mass_mev=inputs['pdg']['masses']['kstar_charged']['value_mev']
        - inputs['pdg']['masses']['kstar_neutral']['value_mev'],
        rho_in_omega=[rho_in_omega.real, rho_in_omega.imag], rho_in_omega_modulus=abs(rho_in_omega),
        rho_in_omega_zero_width=rho_in_omega_zero_width,
        real_diagonalization_angle_rad=real_diagonalization_angle,
        rho_in_phi=[rho_in_phi.real, rho_in_phi.imag], rho_in_phi_zero_width=rho_in_phi_zero_width,
        photon_exchange_rho_omega_gev2_estimate=photon_rho_omega,
        photon_exchange_rho_phi_gev2_estimate=photon_rho_phi,
        omega_pipi_branching=omega_pipi_branching, tree_rho_width_gev=tree_rho_width,
        residue_isospin_breaking_gev={'rho': breaking['rho'], 'omega': omega_breaking},
        isospin_limit_residues_gev=residues, leptonic_relative_change=leptonic_relative_change,
        leptonic_relative_sigma=leptonic_sigma, radiative_rows=radiative,
        babar_relative_offdiagonal_at_1p37=babar_relative_offdiagonal,
        phi_charged_factor=charged_factor, phi_neutral_factor=neutral_factor,
        phi_ratio_change=phi_ratio_change, phi_ratio_change_linear=phi_ratio_change_linear,
        phi_ratio_change_with_photon_exchange_estimate=phi_ratio_change_with_photon_estimate,
        mass_splitting_coefficient_relative_change_charged_calibration=calibration_change,
        phi_shifted_predictions=shifted.tolist(),
        phi_shifted_cholesky_residuals=shifted_cholesky.tolist(),
        phi_shifted_symmetric_residuals=shifted_symmetric.tolist(),
        phi_shifted_chi2=float(shifted_symmetric @ shifted_symmetric),
        phi_endpoint_symmetric_residuals=diagnostic['symmetric_whitened_channel_components'],
        phi_endpoint_chi2=diagnostic['rotated_chi2'],
        phi_residual_decomposition=decomposition,
        phi_best_split_branching_coefficient=float(split_coefficient),
        phi_best_split_ratio_change=float((1 + split_coefficient)/(1 - split_coefficient) - 1),
        vpp_option_ratio_change_at_task_3_2=vpp_ratio_change_task_3_2)


def main():
    """Print the symbolic results and the endpoint numbers."""
    result = derive()
    print('Neutral mass matrix (rho0, omegaI, phiI):')
    print(sp.latex(result['neutral_mass_matrix']))
    print('(M^2)_{rho omegaI} =', sp.latex(sp.factor(result['rho_omega_ideal'])))
    print('(M^2)_{rho phiI} =', sp.latex(sp.factor(result['rho_phi_ideal'])))
    for name, value in result['charged_mass_squared'].items():
        print(f'm^2_{name} =', sp.latex(value))
    print('pole component of rho0 in state b:', result['pole_component'])
    for name, value in result['residue_isospin_breaking'].items():
        print(f'isospin-breaking residue {name}:', sp.latex(sp.factor(value)))
    print('phi K+K-/K0K0bar ratio - 1 at first order in c:', sp.latex(result['phi_ratio_first_order']))
    numbers = evaluate(result)
    print(json.dumps(numbers, indent=2, default=float))


if __name__ == '__main__':
    main()
