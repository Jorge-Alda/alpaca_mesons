"""Fixed-scenario measured-width Jacobians of the bounded odd VJP/VVP action.

Run .venv/bin/python scripts/derive_odd_low_energy_jacobians.py.
GeV, tree level; no physical odd fit, OPE row, naturalness prior or branch
pooling. LO scalar masses enter local kernels; physical masses enter phase
space as a labelled kinematic replacement. Data and input constraints are
used once. All mass/source alternatives remain independent scenarios.
"""
from pathlib import Path
import argparse
import contextlib
import hashlib
import json
import math
import numpy as np
import scipy.linalg as la
import sympy as sp
from diagnose_vector_stage4 import module, E
from derive_vector_mass_extension import derive as derive_mass
from derive_odd_bounded_hilbert_series import image_dimension, relation_matrix

ROOT = Path(__file__).resolve().parents[1]
F, OPERATOR_MV = .093, .770
TOL = 1e-10
WIDTH_KEYS = ('universal', 'single_trace_spurion', 'single_trace_spurion_sign_alternative',
              'plus_double_trace_spurion')
HYPOTHESES = ('universal', 'single_trace_positive', 'single_trace_negative', 'double_trace')
SEED_ORIGINAL = np.array([.013, -.017, .019, .023, -.029, .031, .037, -.041, .043, .047, -.053])
BROAD_RHO_CHANNELS = {'rho0_pi0_gamma', 'rho0_eta_gamma', 'rhoplus_piplus_gamma'}
CONSTRAINED_KSTAR_CHANNELS = {'Kstar0_K0_gamma', 'Kstarplus_Kplus_gamma'}


def source_parameters(hypothesis):
    return ['F_V'] + ([] if hypothesis == 'universal' else ['z_V']) + (
        ['y_V'] if hypothesis == 'double_trace' else [])


def frozen_kernel():
    derived = derive_mass()
    mr, mk = sp.symbols('m_rho m_Kstar', positive=True)
    K = derived['calibrated'].subs({derived['r']: mr**2, derived['k']: mk**2})
    order = (mr, mk, derived['Delta'], derived['kap'], derived['xl'], derived['xs'])
    return sp.lambdify(order, K, 'numpy'), K


MASS_KERNEL, SYMBOLIC_MASS_MATRIX = frozen_kernel()


def action_sector(point, hypothesis, assumptions):
    mr, mk, delta, kappa = point[:4]
    K = np.asarray(MASS_KERNEL(mr, mk, delta, kappa,
                              assumptions['xl_gev2'], assumptions['xs_gev2']), dtype=float)
    squared = np.linalg.eigvalsh(K)
    if np.any(squared <= 0):
        raise ValueError('mass action is not positive')
    theta = .5*np.arctan2(2*K[0, 1], K[1, 1]-K[0, 0])
    sectors = module('vector_sector')
    generators = sectors.physical_generators(theta)
    FV = point[4]
    z = point[5] if hypothesis != 'universal' else 0.
    y = point[6] if hypothesis == 'double_trace' else 0.
    chi_plus = np.diag([assumptions['xl_gev2']]*2+[assumptions['xs_gev2']])
    residues = sectors.source_residues(generators, chi_plus, FV, z, y)
    sector = sectors.VectorSector((mr, *np.sqrt(squared)), generators, residues)
    return sector, mk


def channel_rows(sector, mk, tree, physical, e):
    api = module('odd_bounded_amplitudes')
    inputs = module('vector_sector_amplitudes').VectorAmplitudeInputs
    rho, omega, phi = sector.generators
    masses = dict(zip(('rho', 'omega', 'phi'), sector.masses))
    chi = tree['pi0']['chi']
    result = {}
    def radiative(label, TV, TP, source_mass, physical_mass, mv, orientation):
        x = inputs(F, OPERATOR_MV, source_mass, TP, chi, sector)
        row = api.radiative_gradients(x, TV, mv, e)
        if orientation == 'V_to_P':
            factor = (mv**2-physical_mass**2)**3/(96*math.pi*mv**3)
        else:
            factor = (physical_mass**2-mv**2)**3/(32*math.pi*physical_mass**3)
        if factor <= 0:
            raise ValueError('closed radiative channel '+label)
        result[label] = dict(row=row, offset=0., phase_factor=factor,
                             source_mass_gev=source_mass, phase_mass_gev=physical_mass,
                             orientation=orientation)
    neutral = [('omega_pi0_gamma', omega, 'pi0', 'omega', 'V_to_P'),
               ('rho0_pi0_gamma', rho, 'pi0', 'rho', 'V_to_P'),
               ('rho0_eta_gamma', rho, 'eta', 'rho', 'V_to_P'),
               ('omega_eta_gamma', omega, 'eta', 'omega', 'V_to_P'),
               ('etap_rho0_gamma', rho, 'etap', 'rho', 'P_to_V'),
               ('etap_omega_gamma', omega, 'etap', 'omega', 'P_to_V'),
               ('phi_eta_gamma', phi, 'eta', 'phi', 'V_to_P')]
    for label, TV, pname, vname, orientation in neutral:
        radiative(label, TV, tree[pname]['TP'], tree[pname]['mp'],
                  physical[pname], masses[vname], orientation)
    for label, i, j, pname, mv in (
        ('Kstar0_K0_gamma', 1, 2, 'K_neutral', mk),
        ('Kstarplus_Kplus_gamma', 0, 2, 'K_charged', mk),
        ('rhoplus_piplus_gamma', 0, 1, 'pi_charged', masses['rho'])):
        mp = np.sqrt((chi[i, i]+chi[j, j])/2)
        radiative(label, E(j, i), E(i, j), mp, physical[pname], mv, 'V_to_P')
    alpha = e**2/(4*math.pi)
    for pname in ('pi0', 'eta', 'etap'):
        x = inputs(F, OPERATOR_MV, tree[pname]['mp'], tree[pname]['TP'], chi, sector)
        result[pname+'_gamma_gamma'] = dict(row=api.tff_gradients(x), offset=api.wzw_offset(x),
            phase_factor=math.pi*alpha**2*physical[pname]**3/4,
            source_mass_gev=tree[pname]['mp'], phase_mass_gev=physical[pname], orientation='P_to_gamma_gamma')
    return result


def data_units(observation):
    unit = observation.get('unit', 'fraction')
    scale = {'fraction': 1., 'GeV': 1., 'MeV': 1e-3, 'keV': 1e-6, 'eV': 1e-9}[unit]
    return observation['value']*scale, observation['sigma']*scale


def measured_predictions(channels, observations, totals, coefficients):
    values, jacobian, amplitudes = [], [], []
    for observation in observations:
        channel = channels[observation['channel']]
        row = channel['row']
        G = channel['offset']+row@coefficients
        k = channel['phase_factor']
        total = totals[observation['total_width_key']] if observation['observable'] == 'branching_fraction' else 1.
        values.append(k*G**2/total)
        jacobian.append(2*k*G*row/total)
        amplitudes.append(G)
    return np.asarray(values), np.asarray(jacobian), np.asarray(amplitudes)


def even_observables(sector, mk, assumptions):
    masses = np.array(sector.masses)
    ratios = (assumptions['electron_mass_gev']/masses)**2
    widths = 8*math.pi*assumptions['alpha_em']**2*np.array(sector.residues)**2/(3*masses)
    widths *= (1+2*ratios)*np.sqrt(1-4*ratios)
    return np.r_[masses[0], mk, masses[1:], widths]


def rank(matrix, reference=0.):
    matrix = np.asarray(matrix, dtype=float)
    if matrix.size == 0:
        return 0, np.array([])
    norms = np.linalg.norm(matrix, axis=0)
    cutoff = TOL*max(reference, np.max(norms, initial=0))
    keep = norms > cutoff
    if not np.any(keep):
        return 0, np.array([])
    balanced = matrix[:, keep]/norms[keep]
    singular = la.svdvals(balanced)
    return int(np.sum(singular > TOL*max(1., singular[0]))), singular


def information_diagnostic(jacobian, selected, nuisance=None):
    matrix = jacobian[:, selected]
    fixed_rank, _ = rank(matrix)
    singular = la.svdvals(matrix)
    summary = dict(rank=fixed_rank, whitened_singular_values=singular.tolist(),
                   rank_method='nonzero column normalization; relative SVD tolerance 1e-10; row whitening retained')
    if nuisance is not None:
        nrank, _ = rank(nuisance)
        joint_rank, _ = rank(np.column_stack((nuisance, matrix)))
        n = np.linalg.norm(nuisance, axis=0)
        u, sv, _ = la.svd(nuisance[:, n > 0]/n[n > 0], full_matrices=False)
        projected = matrix-u[:, :nrank]@(u[:, :nrank].T@matrix)
        marginal_rank = joint_rank-nrank
        summary.update(nuisance_rank=nrank, joint_rank=joint_rank,
                       marginal_odd_rank=marginal_rank,
                       nuisance_projected_singular_values=la.svdvals(projected).tolist())
    # Representative measured combinations: pivot columns span the local image.
    if fixed_rank:
        _, _, pivot = la.qr(matrix, pivoting=True, mode='economic')
        chosen = pivot[:fixed_rank]
        basis = matrix[:, chosen]
        coefficients = la.lstsq(basis, matrix)[0]
        residual = np.linalg.norm(basis@coefficients-matrix)/max(1., np.linalg.norm(matrix))
        if residual > 5e-9:
            raise AssertionError('measurement-column reconstruction failed')
        summary.update(representative_parent_indices=[int(selected[i]) for i in chosen],
                       column_reconstruction_relative_residual=float(residual),
                       identifiable_combination_coefficients=coefficients.tolist(),
                       combination_column_indices=selected.tolist(),
                       raw_null_generator_formula='for each nonpivot k: e_k minus sum_l coefficients[l,k] e_pivot_l; divide out inherited action relations for physical null directions')
    return summary


def symbolic_observable_derivation():
    coupling = sp.Symbol('g', real=True)
    amplitude = sp.Function('G')(coupling)
    k, total, a = sp.symbols('k Gamma_total A_j', positive=True)
    observable = k*amplitude**2/total
    odd_derivative = sp.diff(observable, coupling).subs(sp.diff(amplitude, coupling), a)
    total_derivative = sp.diff(observable, total)
    return observable, odd_derivative, total_derivative


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--measurements', default='data/odd_parity/odd_low_energy_measurements.json')
    parser.add_argument('--limit-scenarios', type=int, default=0,
                        help='debug only; final artifact requires all 64 scenarios')
    args = parser.parse_args()
    measurement_path = ROOT/args.measurements
    data = json.loads(measurement_path.read_text())
    observations = [item for item in data['observations'] if item.get('status', 'included') != 'deferred']
    if not observations:
        raise ValueError('no accepted primary observations')
    for source in data.get('sources', {}).values():
        if 'local_path' in source and 'sha256' in source:
            if hashlib.sha256((ROOT/source['local_path']).read_bytes()).hexdigest() != source['sha256']:
                raise ValueError('measurement source hash mismatch')
    calibration = json.loads((ROOT/'data/vector_sector/stage3_results.json').read_text())
    pdg = json.loads((ROOT/'data/vector_sector/pdg2024.json').read_text())
    alda = json.loads((ROOT/'data/vector_sector/alda_arxiv_v1.json').read_text())
    assumptions = calibration['assumptions']
    e = np.sqrt(4*math.pi*assumptions['alpha_em'])
    tree = module('alda_neutral_inputs').alda_neutral_inputs(mass_scheme='tree', F=F, fpi=F)
    bounded = module('odd_bounded_amplitudes')
    ids = list(bounded.parent_ids())
    embedding = bounded.original_embedding()
    seed = embedding@SEED_ORIGINAL
    measured_values, errors = np.array([data_units(item) for item in observations]).T
    if np.any(errors <= 0):
        raise ValueError('measurement errors must be positive')
    totals_metadata = data.get('total_widths', {})
    used_totals = sorted(set(item['total_width_key'] for item in observations
                            if item['observable'] == 'branching_fraction'))
    totals = {key: totals_metadata[key]['value_gev'] for key in used_totals}
    total_errors = {key: totals_metadata[key]['sigma_gev'] for key in used_totals}
    if any(value <= 0 for value in [*totals.values(), *total_errors.values()]):
        raise ValueError('positive total widths and errors required')
    leading = json.loads((ROOT/'data/vector_sector/odd_leading_action_count.json').read_text())
    higher = json.loads((ROOT/'data/vector_sector/odd_three_flavor_action_quotient.json').read_text())
    exact_relations = relation_matrix(ids, leading['relations']+higher['all_higher_relations'])
    R = np.asarray(exact_relations, dtype=float)
    filtrations = {}
    for scheme in ('all_derivatives_soft', 'hard_vector_derivatives_promoted'):
        filtrations[scheme] = {}
        for level in (0, 1, 2):
            selected = bounded.retained_columns(scheme, level)
            filtrations[scheme][str(level)] = dict(columns=selected.tolist(), raw_parent_count=len(selected),
                parent_quotient_dimension=int(image_dimension(list(selected), exact_relations)))
    results = []
    for mass_case in calibration['scenarios']:
        mr = pdg['masses'][mass_case['rho_input']]['value_mev']/1000
        mk = pdg['masses'][mass_case['kstar_input']]['value_mev']/1000
        branches = mass_case['extended_mass_branches']
        sorted_branches = sorted(range(len(branches)), key=lambda i: abs(branches[i]['theta_deg']))
        for width_case in mass_case['width_tests']:
            for branch_label, branch_index in zip(('near_ideal', 'other'), sorted_branches):
                branch = branches[branch_index]
                sources = next(item for item in width_case['fixed_angle_cases']
                               if item['case'] == f'extended_mass_branch_{branch_index}')
                for hypothesis, source_key in zip(HYPOTHESES, WIDTH_KEYS):
                    source = sources[source_key]
                    point = np.array([mr, mk, branch['delta_gev2'], branch['kappa'], source['FV_gev']]
                         + ([] if hypothesis == 'universal' else [source['z_gev_inverse']])
                         + ([source['y_gev_inverse']] if hypothesis == 'double_trace' else []))
                    sector, mk_model = action_sector(point, hypothesis, assumptions)
                    held_sector = module('vector_sector').VectorSector.from_stage3(
                        rho_input=mass_case['rho_input'], kstar_input=mass_case['kstar_input'],
                        phi_input=width_case['phi_input'], source_hypothesis=hypothesis, branch=branch_label)
                    np.testing.assert_allclose(sector.response(), held_sector.response(), rtol=2e-12, atol=2e-13)
                    channels = channel_rows(sector, mk_model, tree, alda['physical_masses_gev'], e)
                    rows = np.array([channels[item['channel']]['row'] for item in observations])
                    predictions, J, amplitudes = measured_predictions(channels, observations, totals, seed)
                    if np.min(np.abs(amplitudes)) <= 1e-12:
                        raise AssertionError('reference seed lies on a zero-amplitude locus')
                    all_rows = np.array([item['row'] for item in channels.values()])
                    relation_error = np.max(np.abs(all_rows@R.T))/max(1., np.linalg.norm(all_rows))
                    if relation_error > 1e-11:
                        raise AssertionError('LO scalar-source relation defect')
                    mass_errors = [pdg['masses'][key]['sigma_mev']/1000
                                   for key in (mass_case['rho_input'], mass_case['kstar_input'], 'omega', 'phi')]
                    electron_errors = [pdg['leptonic_widths'][key]['sigma_kev']*1e-6
                                       for key in ('rho', 'omega', 'phi')]
                    if width_case['phi_input'] == 'phi_ee_average':
                        electron_errors[2] = pdg['alternatives']['phi_ee_average']['sigma_kev']*1e-6
                    # Read the actual alternative label rather than pooling it.
                    if width_case['phi_input'] != 'phi_evaluation':
                        electron_errors[2] = pdg['alternatives']['phi_ee_average']['sigma_kev']*1e-6
                    even_errors = np.array(mass_errors+electron_errors)
                    parameter_names = ['m_rho', 'm_Kstar', 'Delta_V', 'kappa_V']+source_parameters(hypothesis)
                    scales = np.array([.001, .001, .01, .01, .01]+(
                        [] if hypothesis == 'universal' else [.1])+([.1] if hypothesis == 'double_trace' else []))
                    def local_predictions(p):
                        local_sector, local_mk = action_sector(p, hypothesis, assumptions)
                        local_channels = channel_rows(local_sector, local_mk, tree, alda['physical_masses_gev'], e)
                        widths = measured_predictions(local_channels, observations, totals, seed)[0]
                        return np.r_[widths/errors, even_observables(local_sector, local_mk, assumptions)/even_errors]
                    nuisance_even = []
                    convergence = 0.
                    for i, scale in enumerate(scales):
                        columns = []
                        for fraction in (1e-3, 5e-4):
                            step = fraction*scale
                            hi, lo = point.copy(), point.copy()
                            hi[i] += step; lo[i] -= step
                            columns.append((local_predictions(hi)-local_predictions(lo))*scale/(2*step))
                        convergence = max(convergence, np.linalg.norm(columns[0]-columns[1])
                                          /max(1., np.linalg.norm(columns[1])))
                        nuisance_even.append(columns[1])
                    if convergence > 1e-6:
                        raise AssertionError('even/source derivative step convergence')
                    ndata, ntotal = len(observations), len(used_totals)
                    full_J = np.vstack((J/errors[:, None], np.zeros((7+ntotal, len(ids)))))
                    nuisance = np.zeros((ndata+7+ntotal, len(point)+ntotal))
                    nuisance[:ndata+7, :len(point)] = np.array(nuisance_even).T
                    for k, key in enumerate(used_totals):
                        for i, observation in enumerate(observations):
                            if observation.get('total_width_key') == key:
                                nuisance[i, len(point)+k] = -predictions[i]*total_errors[key]/totals[key]/errors[i]
                        nuisance[ndata+7+k, len(point)+k] = 1.
                    constraint_rank, _ = rank(nuisance[ndata:, :])
                    if constraint_rank != nuisance.shape[1]:
                        raise AssertionError('mass/electronic/total-width constraints do not identify their local nuisance columns')
                    nested = {}
                    for scheme, levels in filtrations.items():
                        nested[scheme] = {}
                        for level, info in levels.items():
                            selected = np.array(info['columns'])
                            diagnostic = information_diagnostic(full_J, selected, nuisance)
                            diagnostic.update(raw_parent_count=info['raw_parent_count'],
                                parent_quotient_dimension=info['parent_quotient_dimension'],
                                unmeasured_action_dimensions=info['parent_quotient_dimension']-diagnostic['marginal_odd_rank'])
                            nested[scheme][level] = diagnostic
                    original_diagnostic = information_diagnostic(full_J@embedding, np.arange(11), nuisance)
                    controls = {}
                    for label, excluded in (
                        ('without_broad_rho', BROAD_RHO_CHANNELS),
                        ('without_constrained_Kstar', CONSTRAINED_KSTAR_CHANNELS),
                        ('narrow_nonconstrained_summary_control', BROAD_RHO_CHANNELS|CONSTRAINED_KSTAR_CHANNELS)):
                        chosen_rows = [i for i, item in enumerate(observations) if item['channel'] not in excluded]
                        control_totals = set(observations[i]['total_width_key'] for i in chosen_rows
                                             if observations[i]['observable'] == 'branching_fraction')
                        chosen_rows += list(range(ndata, ndata+7))
                        chosen_rows += [ndata+7+i for i, key in enumerate(used_totals) if key in control_totals]
                        control_nuisance = nuisance[chosen_rows, :]
                        control_J = full_J[chosen_rows, :]
                        control_nested = {}
                        for scheme, levels in filtrations.items():
                            control_nested[scheme] = {}
                            for level, info in levels.items():
                                result = information_diagnostic(control_J, np.array(info['columns']), control_nuisance)
                                control_nested[scheme][level] = {key: result[key] for key in
                                    ('rank', 'nuisance_rank', 'joint_rank', 'marginal_odd_rank')}
                        controls[label] = dict(measurement_labels=[item['channel'] for item in observations
                                                                  if item['channel'] not in excluded],
                            nested=control_nested, covariance_status='still conditional diagonal summaries; exclusion does not supply missing covariance')
                    # Width differentiation is checked with a simultaneous
                    # finite displacement of all original couplings.
                    direction = embedding@np.cos(np.arange(11))
                    step = 1e-7
                    numerical = (measured_predictions(channels, observations, totals, seed+step*direction)[0]
                               - measured_predictions(channels, observations, totals, seed-step*direction)[0])/(2*step)
                    np.testing.assert_allclose(numerical, J@direction, rtol=1e-8, atol=1e-10)
                    fixed_rank, _ = rank(rows)
                    measured_rank, _ = rank(J/errors[:, None])
                    if fixed_rank != measured_rank:
                        raise AssertionError('nonzero width/amplitude local rank correspondence')
                    results.append(dict(rho_input=mass_case['rho_input'], kstar_input=mass_case['kstar_input'],
                        phi_input=width_case['phi_input'], branch=branch_label, source_hypothesis=hypothesis,
                        signed_residues_gev=list(sector.residues), measurement_labels=[item['channel'] for item in observations],
                        reference_amplitudes=amplitudes.tolist(), reference_predictions=predictions.tolist(),
                        fixed_measured_amplitude_rank=fixed_rank, original_11=original_diagnostic, nested=nested,
                        data_adequacy_controls=controls, locally_identified_nuisance_constraint_rank=constraint_rank,
                        measurement_amplitude_rows=rows.tolist(), whitened_odd_measurement_rows=(J/errors[:, None]).tolist(),
                        nuisance_parameter_names=parameter_names+['total_width:'+key for key in used_totals],
                        even_parameter_column_scales=scales.tolist(), whitened_nuisance_rows=nuisance.tolist(),
                        LO_relation_relative_residual=float(relation_error), nuisance_step_relative_change=float(convergence)))
                    if args.limit_scenarios and len(results) >= args.limit_scenarios:
                        break
                if args.limit_scenarios and len(results) >= args.limit_scenarios: break
            if args.limit_scenarios and len(results) >= args.limit_scenarios: break
        if args.limit_scenarios and len(results) >= args.limit_scenarios: break
    if not args.limit_scenarios:
        assert len(results) == 64
    paths = [str(measurement_path.relative_to(ROOT)), 'scripts/derive_odd_low_energy_jacobians.py',
             'scripts/implementation/odd_bounded_amplitudes.py', 'data/vector_sector/odd_cubic_parent_kernels.json',
             'data/vector_sector/stage3_results.json', 'data/vector_sector/pdg2024.json',
             'data/vector_sector/alda_arxiv_v1.json', 'data/vector_sector/odd_leading_action_count.json',
             'data/vector_sector/odd_three_flavor_action_quotient.json']
    output = ROOT/'data/vector_sector/odd_low_energy_jacobians.json'
    if args.limit_scenarios:
        output = ROOT/'data/vector_sector/odd_low_energy_jacobian_debug.json'
    artifact = dict(scope='fixed-scenario local measured-width Jacobians; no physical odd fit',
        complete_64_scenario_run=not bool(args.limit_scenarios), operator_order=ids, filtrations=filtrations,
        reference_original_coefficients=SEED_ORIGINAL.tolist(), reference_status='deterministic nonzero point; no OPE, fit or prior',
        measured_values=measured_values.tolist(), measurement_errors=errors.tolist(), measurement_snapshot=args.measurements,
        vertex_scalar_mass_scheme='LO Alda eigenmasses and charged source-EOM masses',
        phase_space_scalar_mass_scheme='physical Alda input-table masses; labelled kinematic replacement',
        domain='real photons and on-shell light vectors; no hard TFF or three-pion input',
        nuisance_policy='original mass/electronic/total-width observations once; diagonal summary errors conditional on absent cross covariance',
        alda_uncertainties='fixed central generators/quark spurion; posterior covariance and truncation effects unpropagated',
        zero_radiative_amplitude_warning='at zero odd couplings radiative width derivatives vanish; this singular locus is not a global identifiability result',
        scenarios=results, sha256={path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in paths})
    output.write_text(json.dumps(artifact, indent=2, allow_nan=False)+'\n')
    if not args.limit_scenarios:
        transcript = ROOT/'data/vector_sector/odd_low_energy_jacobians.md'
        with transcript.open('w') as handle, contextlib.redirect_stdout(handle):
            print('# Fixed-scenario low-energy measurement Jacobians\n')
            for label, value in zip(('Measured BR', 'Odd-coupling derivative', 'Total-width derivative'),
                                    symbolic_observable_derivation()):
                print(label+':\n$$\n'+sp.latex(value)+'\n$$\n')
            print('The raw columns, whitened measured Jacobians, nuisance Jacobians and all separate scenarios are in the JSON artifact.\n')
            for scheme in filtrations:
                for level in ('0', '1', '2'):
                    ranks = sorted(set(item['nested'][scheme][level]['marginal_odd_rank'] for item in results))
                    print(f'{scheme}, ceiling {level}: marginal odd ranks {ranks}.\n')
        print(transcript.relative_to(ROOT))
    print(output.relative_to(ROOT))


if __name__ == '__main__':
    main()
