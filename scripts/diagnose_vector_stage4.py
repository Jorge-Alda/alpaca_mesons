"""Conditional action-kernel ranks; inherited OPE rows are NOT rematched.

Run: .venv/bin/python scripts/diagnose_vector_stage4.py
No data fit, uncertainty pooling, or physical-angle inference is performed.
"""
from pathlib import Path
import hashlib
import importlib
import json
import sys
import types
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
# Avoid eager import of unrelated implementation modules.
PACKAGE = '_vector_stage4_impl'
pkg = types.ModuleType(PACKAGE)
pkg.__path__ = [str(ROOT / 'scripts/implementation')]
sys.modules[PACKAGE] = pkg

def module(name):
    return importlib.import_module(PACKAGE + '.' + name)

odd = module('odd_parity_rcht')
ORDER = [f'c{i}' for i in range(1, 8)] + [f'd{i}' for i in range(1, 5)]
TOL = 1e-11
OPE = np.zeros((5, 11))
OPE[0, [0, 2]] = [1, 4]
OPE[1, [0, 1, 4]] = [1, -1, 1]
OPE[2, [4, 5]] = [1, -1]
OPE[3, [7, 8]] = [1, 8]
OPE[4, 9] = 1
NULLS = np.eye(11)[:, [6, 10]].T.tolist() + [[0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0]]
BASIS = [odd.OddParityCouplings(odd.VJPCouplings(*r[:7]), odd.VVPCouplings(*r[7:]))
         for r in np.eye(11)]
ZERO = odd.OddParityCouplings(odd.VJPCouplings(*([0]*7)), odd.VVPCouplings(*([0]*4)))

def sensitivity(function):
    offset = function(ZERO)
    values = np.asarray([function(c)-offset for c in BASIS], dtype=complex)
    assert np.max(np.abs(values.imag)) < 1e-12
    row = values.real
    for vector in NULLS:
        assert abs(row @ vector) < 1e-11 * max(1, np.linalg.norm(row))
    return row

def diagnostic(rows):
    matrix = np.vstack(rows)
    norms = np.linalg.norm(matrix, axis=1)
    normalized = matrix[norms > 0] / norms[norms > 0, None]
    singular = np.linalg.svd(normalized, compute_uv=False)
    return dict(rank=int(np.count_nonzero(singular > TOL)),
                singular_values=singular.tolist(), row_norms=norms.tolist(),
                normalization='each nonzero sensitivity row divided by its Euclidean norm',
                absolute_singular_value_tolerance=TOL)

def E(i, j):
    result = np.zeros((3, 3)); result[i, j] = 1
    return result

def rows_for(sector, pseudoscalars, masses, *, F=.093, MV=.770, e=.3028):
    api = module('vector_sector_amplitudes')
    rho, omega, phi = sector.generators
    chi = pseudoscalars['pi0']['chi']
    channels = [('omega pi0', omega, 'pi0', 'omega'), ('rho0 pi0', rho, 'pi0', 'rho'),
                ('rho0 eta', rho, 'eta', 'rho'), ('omega eta', omega, 'eta', 'omega'),
                ('etap rho0', rho, 'etap', 'rho'), ('etap omega', omega, 'etap', 'omega'),
                ('phi eta', phi, 'eta', 'phi')]
    rows, names = [], []
    def radiative(label, TV, TP, mp, mv):
        x = api.VectorAmplitudeInputs(F, MV, mp, TP, chi, sector)
        rows.append(sensitivity(lambda c: sum(api.radiative_components(c, x, TV, mv, e))))
        names.append(label)
    for label, TV, pname, vname in channels:
        p = pseudoscalars[pname]
        radiative(label, TV, p['TP'], p['mp'], masses[vname])
    for label, TV, TP, mp, mv in (
        ('Kstar0 K0', E(2,1), E(1,2), masses['K0'], masses['Kstar0']),
        ('Kstar+ K+', E(2,0), E(0,2), masses['K+'], masses['Kstar+']),
        ('rho+ pi+', E(1,0), E(0,1), masses['pi+'], masses['rho+'])):
        radiative(label, TV, TP, mp, mv)
    real, virtual = [], []
    for pname in ('pi0','eta','etap'):
        p = pseudoscalars[pname]
        x = api.VectorAmplitudeInputs(F, MV, p['mp'], p['TP'], p['chi'], sector)
        real.append(sensitivity(lambda c: sum(api.tff_components(c, x, 0, 0))))
        for s in (-1., -5.):
            virtual.append(sensitivity(lambda c: sum(api.tff_components(c, x, s, 0))))
    return dict(channel_labels=names,
                radiative=diagnostic([*OPE, *rows]),
                radiative_and_real_tff=diagnostic([*OPE, *rows, *real]),
                radiative_and_real_and_virtual_tff=diagnostic([*OPE, *rows, *real, *virtual]),
                radiative_rows=np.asarray(rows).tolist(), real_tff_rows=np.asarray(real).tolist(),
                virtual_tff_rows=np.asarray(virtual).tolist())


def baseline_check():
    """Reproduce old illustrative 7/common, 8/split with independent chi."""
    sector_api = module('vector_sector')
    rho, omega, phi = sector_api.physical_generators(0)
    eta8 = np.diag([1,1,-2])/np.sqrt(12)
    eta1 = np.eye(3)/np.sqrt(3)
    chi = np.diag([.1396**2, .1396**2, 2*.4937**2-.1396**2])
    p = {name:dict(TP=TP, chi=chi, mp=mp) for name,TP,mp in (
        ('pi0',rho,.1396), ('eta',2*np.sqrt(2)/3*eta8-eta1/3,.5479),
        ('etap',eta8/3+2*np.sqrt(2)/3*eta1,.9578))}
    masses=dict(rho=.7753,omega=.7827,phi=1.0195,K0=.4976,
                **{'K+':.4937,'pi+':.1396,'rho+':.7753,'Kstar0':.8956,'Kstar+':.8922})
    result={}
    for name, internal in (('common',(.770,)*3),('split',(.7753,.7827,1.0195))):
        generators=(rho,omega,phi)
        sector=sector_api.VectorSector(internal,generators,
            sector_api.source_residues(generators,np.zeros((3,3)),.13),name)
        r=rows_for(sector,p,masses,F=.092)
        assert r['radiative']['rank']==(7 if name=='common' else 8)
        assert r['radiative_and_real_and_virtual_tff']['rank']==8
        result[name]=r
    return result


def main():
    sectors=module('vector_sector')
    alda=module('alda_neutral_inputs')
    calibration=json.loads((ROOT/'data/vector_sector/stage3_results.json').read_text())
    pdg=json.loads((ROOT/'data/vector_sector/pdg2024.json').read_text())
    pseudoscalars=alda.alda_neutral_inputs(mass_scheme='physical',F=.093,fpi=.093)
    mass_values={name:entry['value_mev']/1000 for name,entry in pdg['masses'].items()}
    # Physical kinematic masses from Alda's arXiv-v1 input table, not tree eigenmasses.
    snapshot=json.loads((ROOT/'data/vector_sector/alda_arxiv_v1.json').read_text())
    pm=snapshot['physical_masses_gev']
    fixed={'K0':pm['K_neutral'],'K+':pm['K_charged'],'pi+':pm['pi_charged']}
    outcomes=[]
    for scenario in calibration['scenarios']:
        rho_input,kstar_input=scenario['rho_input'],scenario['kstar_input']
        masses={**fixed,'rho':mass_values[rho_input], 'omega':mass_values['omega'],
                'phi':mass_values['phi'], 'rho+':mass_values[rho_input],
                'Kstar0':mass_values[kstar_input],'Kstar+':mass_values[kstar_input]}
        for width in scenario['width_tests']:
            for branch in ('near_ideal','other'):
                for hypothesis in ('universal','single_trace_positive','single_trace_negative','double_trace'):
                    sector=sectors.VectorSector.from_stage3(rho_input=rho_input,kstar_input=kstar_input,
                        phi_input=width['phi_input'],source_hypothesis=hypothesis,branch=branch)
                    outcomes.append(dict(rho_input=rho_input,kstar_input=kstar_input,
                        phi_input=width['phi_input'],branch=branch,preferred_if_otherwise_equal=branch=='near_ideal',
                        source_hypothesis=hypothesis, signed_residues_gev=list(sector.residues),
                        **rows_for(sector,pseudoscalars,masses)))
    alda_common=rows_for(sectors.VectorSector.universal(.770,.13),pseudoscalars,
        {**fixed,'rho':mass_values['rho_neutral'],'omega':mass_values['omega'],
         'phi':mass_values['phi'],'rho+':mass_values['rho_charged'],
         'Kstar0':mass_values['kstar_neutral'],'Kstar+':mass_values['kstar_charged']})
    assert alda_common['radiative']['rank']==7
    assert alda_common['radiative_and_real_and_virtual_tff']['rank']==8
    paths=['scripts/diagnose_vector_stage4.py','scripts/implementation/vector_sector.py',
           'scripts/implementation/vector_sector_amplitudes.py',
           'scripts/implementation/alda_neutral_inputs.py','data/vector_sector/stage3_results.json',
           'data/vector_sector/pdg2024.json','data/vector_sector/alda_arxiv_v1.json']
    result=dict(coupling_order=ORDER,scope='action sensitivity diagnostic; no fitted odd-parity couplings',
        inherited_ope='five baseline homogeneous sensitivity rows; NOT a rematched short-distance prescription',
        persistent_null_directions=dict(c7=NULLS[0],d4=NULLS[1],common_c2_c5_c6_shift=NULLS[2]),
        inputs=dict(F_gev=.093,operator_MV_gev=.770,virtualities_gev2=[[0,0],[-1,0],[-5,0]],
                    charged_pseudoscalar_masses_gev=fixed,
                    charged_pseudoscalar_mass_status='Alda arxiv-v1 input table; fixed central kinematic inputs',
                    calibrated_external_vector_mass_policy='rho+ uses selected rho_input; charged and neutral Kstar both use selected kstar_input, exact-isospin vector action',
                    vector_spurion_policy='even vector mass/source block uses isospin chi_plus from stage3; odd mass spurion and neutral P eigenvectors retain full Alda isospin; partial isospin truncation, rho mixing omitted',
                    alda_common_control_external_mass_policy='measured charged/neutral external masses retained as explicitly labeled old reference',
                    neutral_pseudoscalar_scheme='physical masses with Alda arxiv-v1 eigenvectors and independent chi',
                    uncertainties='fixed central inputs; no pooled alternatives or model/input uncertainty inference'),
        sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        baseline_consistency=baseline_check(),alda_universal_common_pole_control=alda_common,
        scenarios=outcomes)
    output=ROOT/'data/vector_sector/stage4_rank_results.json'
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(output.relative_to(ROOT))

if __name__=='__main__':
    main()
