"""Conditional original-operator fit of widths, vector Dalitz and spacelike TFFs.

Run: PYTHONPATH=scripts .venv/bin/python scripts/fit_joint_dalitz_tff.py
GeV units; fixed preliminary even sector, tree poles, diagonal published errors.
Complete CLEO/BaBar ranges diagnose this truncation; they do not establish its
validity at hard virtuality. Canonical coefficients are negated for the legacy
amplitude backend, including anomaly interference. Shared total widths and the
published BaBar 2.3 percent normalization are profiled with Gaussian priors.
"""
from pathlib import Path
import hashlib
import json
import math
import unittest
import numpy as np
from scipy.linalg import svd, svdvals
from scipy.optimize import least_squares
from scipy.stats import chi2 as chi2_distribution
import fit_vector_dalitz as dalitz
import calibrate_odd_low_energy as calibration
from implementation.vector_sector_amplitudes import VectorAmplitudeInputs, tff_components, radiative_components
from implementation.radiative_odd_parity import _couplings_from_mapping
from implementation.odd_parity_data import load_cleo_tff, load_babar_pi0_table11

ROOT = Path(__file__).resolve().parents[1]
NAMES = [f'c{i}' for i in range(1,8)]+[f'd{i}' for i in range(1,5)]


def affine_amplitude(function):
    """Evaluate exact affine dependence on canonical original coefficients."""
    zero = dict.fromkeys(NAMES,0.)
    offset = complex(function(_couplings_from_mapping(zero)))
    row = []
    for name in NAMES:
        # Legacy coefficients have the opposite orientation to canonical ones.
        row.append(complex(function(_couplings_from_mapping(dict(zero,**{name:-1.}))))-offset)
    return offset,np.asarray(row)


def prepare():
    """Load pinned measurements and construct exact amplitude coefficient rows."""
    parameters,pion_inputs = dalitz.preliminary_inputs()
    even_assumptions=json.loads((ROOT/'data/vector_sector/stage3_results.json').read_text())['assumptions']
    alda=json.loads((ROOT/'data/vector_sector/alda_2026nrd_arxiv_v1.json').read_text())
    tree=calibration.J.module('alda_neutral_inputs').alda_2026nrd_neutral_inputs(mass_scheme='tree',F=calibration.J.F,fpi=calibration.J.F)
    inputs={name:VectorAmplitudeInputs(calibration.J.F,calibration.J.OPERATOR_MV,item['mp'],item['TP'],item['chi'],pion_inputs.sector) for name,item in tree.items()}
    measurements=json.loads((ROOT/'data/odd_parity/odd_low_energy_measurements.json').read_text())
    observations=measurements['observations']
    channels=calibration.FastChannels(tree,alda['physical_masses_gev'],math.sqrt(4*math.pi*even_assumptions['alpha_em']))(pion_inputs.sector,parameters['m_Kstar_GeV'])
    embedding=calibration.model_maps()['original_11']
    width_rows=np.array([-channels[item['channel']]['row']@embedding for item in observations])
    width_offsets=np.array([channels[item['channel']]['offset'] for item in observations])
    width_phase=np.array([channels[item['channel']]['phase_factor'] for item in observations])
    total_names=list(measurements['total_widths'])
    total_values=np.array([measurements['total_widths'][name]['value_gev'] for name in total_names])
    total_sigmas=np.array([measurements['total_widths'][name]['sigma_gev'] for name in total_names])
    total_index=np.array([total_names.index(item['total_width_key']) if item['observable']=='branching_fraction' else -1 for item in observations])
    width_values=np.array([calibration.J.data_units(item)[0] for item in observations])
    width_errors=np.array([calibration.J.data_units(item)[1] for item in observations])
    blocks=[]
    experiments=json.loads((ROOT/'data/odd_parity/vector_dalitz_measurements.json').read_text())['experiments']
    for experiment in experiments:
        rows=[]
        for item in experiment['bins']:
            mass=(item['mass_low_gev']+item['mass_high_gev'])/2
            rows.append(affine_amplitude(lambda couplings:sum(radiative_components(couplings,pion_inputs,pion_inputs.sector.generators[1],pion_inputs.sector.masses[1],1.,mass*mass)))[1])
        denominator=affine_amplitude(lambda couplings:sum(radiative_components(couplings,pion_inputs,pion_inputs.sector.generators[1],pion_inputs.sector.masses[1],1.)))[1]
        blocks.append(dict(kind='dalitz',label=experiment['bibtex_key'],rows=np.array(rows),denominator=denominator,values=np.array([item['value'] for item in experiment['bins']]),errors=np.array([item.get('quoted_sigma',math.hypot(item['stat_sigma'],item['syst_sigma'])) for item in experiment['bins']])))
    cleo=load_cleo_tff().rows
    for meson in ('pi0',):
        rows=cleo[cleo['meson']==meson]
        blocks.append(tff_block('CLEO:1997fho_'+meson,inputs[meson],rows['q2_GeV2'],rows['q2F_GeV'],np.hypot(rows['stat_GeV'],rows['syst_GeV'])))
    babar=load_babar_pi0_table11().rows
    blocks.append(tff_block('BaBar:2009rrj_pi0',inputs['pi0'],babar['q2bar_GeV2'],babar['q2F_MeV']/1000,(babar['q2F_err_low_MeV']+babar['q2F_err_high_MeV'])/2000))
    return locals()


def tff_block(label,inputs,virtualities,values,errors):
    amplitudes=[affine_amplitude(lambda couplings:sum(tff_components(couplings,inputs,-float(virtuality),0.))) for virtuality in virtualities]
    return dict(kind='tff',label=label,q2=np.array(virtualities),offsets=np.array([item[0] for item in amplitudes]),rows=np.array([item[1] for item in amplitudes]),values=np.array(values),errors=np.array(errors))


def block_prediction(block,coefficients,normalization=0.):
    if block['kind']=='dalitz':
        return abs((block['rows']@coefficients)/(block['denominator']@coefficients))**2
    prediction=block['q2']*abs(block['offsets']+block['rows']@coefficients)
    if block['label'].startswith('BaBar:'):
        prediction=prediction*np.exp(np.log1p(.023)*normalization)
    return prediction


def fit(prepared,kinds):
    """Refit every identifiable original combination, with fixed even inputs."""
    blocks=[block for block in prepared['blocks'] if block['kind'] in kinds]
    amplitude_rows=np.vstack([prepared['width_rows'],*[block['rows'].real for block in blocks]])
    _,singular,right=svd(amplitude_rows,full_matrices=False)
    rank=int(np.sum(singular>singular[0]*1e-10));basis=right[:rank].T
    initial=np.array([prepared['parameters'][name] for name in NAMES])
    coordinates=basis.T@initial
    total_count=len(prepared['total_names'])
    has_babar=any(block['label'].startswith('BaBar:') for block in blocks)
    def residual_parts(point):
        coefficients=basis@point[:rank]
        total_coordinates=point[rank:rank+total_count]
        totals=prepared['total_values']*np.exp(np.log1p(prepared['total_sigmas']/prepared['total_values'])*total_coordinates)
        denominators=np.array([totals[index] if index>=0 else 1. for index in prepared['total_index']])
        predicted=prepared['width_phase']*(prepared['width_offsets']+prepared['width_rows']@coefficients)**2/denominators
        parts={'low_energy':(predicted-prepared['width_values'])/prepared['width_errors'],'total_width_priors':total_coordinates}
        for block in blocks:
            parts[block['label']]=(block_prediction(block,coefficients,point[-1] if has_babar else 0.)-block['values'])/block['errors']
        if has_babar:parts['babar_normalization_prior']=point[-1:]
        return parts
    def residual(point):return np.concatenate(list(residual_parts(point).values()))
    starts=[np.r_[coordinates,np.zeros(total_count+has_babar)]]
    for shift in (-.1,.1):
        changed=initial.copy();changed[0]+=shift;changed[4]+=shift
        starts.append(np.r_[basis.T@changed,np.zeros(total_count+has_babar)])
    fits=[least_squares(residual,start,x_scale='jac',max_nfev=2500,ftol=1e-10,xtol=1e-10,gtol=1e-9) for start in starts]
    best=min(fits,key=lambda item:item.fun@item.fun)
    coefficient=basis@best.x[:rank]
    from summarize_odd_preliminary_fit import main as preliminary
    _,old_null=preliminary(write_output=False)
    derivative=[];normalizations=[]
    for block in blocks:
        matrix=np.column_stack([(block_prediction(block,coefficient+1e-6*direction)-block_prediction(block,coefficient-1e-6*direction))/(2e-6*block['errors']) for direction in old_null.T])
        derivative.append(matrix)
        normalizations.append(block_prediction(block,coefficient)/block['errors'])
    null_matrix=np.vstack(derivative) if derivative else np.zeros((0,4))
    projected=[]
    for matrix,normalization in zip(derivative,normalizations):
        projected.append(matrix-np.outer(normalization,normalization@matrix)/(normalization@normalization))
    parts=residual_parts(best.x);dof=len(best.fun)-len(best.x)
    # c1-c2+c5 multiplies the other photon virtuality in the VJP polynomial.
    combination=np.zeros(11);combination[[0,1,4]]=[1.,-1.,1.]
    combination_gradient=np.r_[combination@basis,np.zeros(total_count+has_babar)]
    _,jacobian_singular,jacobian_right=svd(best.jac,full_matrices=False)
    covariance=(jacobian_right.T/jacobian_singular**2)@jacobian_right
    combination_best=float(combination@coefficient)
    combination_sigma=float(np.sqrt(combination_gradient@covariance@combination_gradient))
    eliminated=int(np.argmax(abs(combination_gradient)))
    kept=np.array([index for index in range(len(best.x)) if index!=eliminated])
    profile=[]
    if rank==8:
        for standardized in (-4.,-2.,-1.,0.,1.,2.,4.):
            fixed=combination_best+standardized*combination_sigma
            def expand(reduced):
                point=best.x.copy();point[kept]=reduced
                point[eliminated]=(fixed-combination_gradient[kept]@reduced)/combination_gradient[eliminated]
                return point
            profiled=least_squares(lambda reduced:residual(expand(reduced)),best.x[kept],x_scale='jac',max_nfev=1000,ftol=1e-10,xtol=1e-10,gtol=1e-9)
            profile.append(dict(c1_minus_c2_plus_c5=fixed,delta_chi2=float(profiled.fun@profiled.fun-best.fun@best.fun),success=bool(profiled.success)))
    return dict(c1_minus_c2_plus_c5=combination_best,local_combination_sigma=combination_sigma if rank==8 else None,combination_identifiable=rank==8,combination_profile=profile,blocks=[block['label'] for block in blocks],odd_rank=rank,chi2=float(best.fun@best.fun),nominal_dof=dof,p_value=float(chi2_distribution.sf(best.fun@best.fun,dof)),block_chi2={name:float(values@values) for name,values in parts.items()},parameters=dict(zip(NAMES,coefficient.tolist())),local_covariance=covariance.tolist(),odd_coordinate_basis=basis.tolist(),covariance_coordinates='orthonormal odd-coordinate basis followed by standardized total widths and optional BaBar normalization',old_null_singular_values=svdvals(null_matrix).tolist(),normalization_projected_old_null_singular_values=svdvals(np.vstack(projected)).tolist() if projected else [],success=bool(best.success),optimality=float(best.optimality),multistart_chi2=[float(item.fun@item.fun) for item in fits],total_width_standardized_shifts=dict(zip(prepared['total_names'],best.x[rank:rank+total_count].tolist())),babar_normalization_standardized_shift=float(best.x[-1]) if has_babar else None,residuals={name:value.tolist() for name,value in parts.items()})


def main():
    prepared=prepare()
    output=dict(assumptions=['fixed preliminary single_trace_negative even sector','original eleven odd operators only','tree-level multi-vector spectral photon response','Alda:2026nrd tree pseudoscalar source masses and matrices','diagonal measurement covariances and independent experiment blocks','full published hard spacelike ranges are conditional model diagnostics','profile shared total widths and BaBar 2.3 percent normalization','bin-center Dalitz spectra','canonical coefficients negated for legacy amplitude backend; common tensor phase invisible in data'],results={label:fit(prepared,kinds) for label,kinds in [('widths',()),('widths_dalitz',('dalitz',)),('widths_tff',('tff',)),('widths_dalitz_tff',('dalitz','tff'))]})
    output['measurements']={block['label']:dict(count=len(block['values']),q2_min=float(block['q2'].min()) if 'q2' in block else None,q2_max=float(block['q2'].max()) if 'q2' in block else None) for block in prepared['blocks']}
    paths=['data/vector_sector/odd_preliminary_fit.json','data/odd_parity/odd_low_energy_measurements.json','data/odd_parity/vector_dalitz_measurements.json','data/odd_parity/cleo_tff.csv','data/odd_parity/babar_pi0_table11.csv']
    output['sha256']={path:hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in paths}
    path=ROOT/'data/vector_sector/joint_dalitz_tff_conditional_fit.json';path.write_text(json.dumps(output,indent=2)+'\n');print(path.relative_to(ROOT))


class JointFitTests(unittest.TestCase):
    def test_affine_canonical_adapter(self):
        parameters,inputs=dalitz.preliminary_inputs()
        function=lambda couplings:sum(tff_components(couplings,inputs,-2.,0.))
        offset,row=affine_amplitude(function)
        coefficients=np.array([parameters[name] for name in NAMES])
        legacy={name:-parameters[name] for name in NAMES}
        self.assertAlmostEqual(abs(offset+row@coefficients-function(_couplings_from_mapping(legacy))),0.,places=12)
    def test_tff_old_null_selection(self):
        _,inputs=dalitz.preliminary_inputs()
        _,row=affine_amplitude(lambda couplings:sum(tff_components(couplings,inputs,-2.,0.)))
        self.assertGreater(abs(row[0]+row[4]),1.)
        self.assertAlmostEqual(abs(row[1]+row[4]+row[5]),0.,places=12)
        self.assertEqual(abs(row[6])+abs(row[10]),0.)

if __name__=='__main__':main()
