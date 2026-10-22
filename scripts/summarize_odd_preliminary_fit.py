"""Export central single-trace negative fit parameters; no new fit scan.

Regenerate the conditional G_V parameter file before exporting the combined fit::

    .venv/bin/python scripts/fit_vector_gv.py data/vector_sector/gv_measurements.json data/vector_sector/gv_conditional_fit.json
    .venv/bin/python scripts/summarize_odd_preliminary_fit.py

Both output files contain only metadata.bibtex_references and parameters.
The standalone G_V file supplies parameters.G_V_GeV to the combined export.
Numerical channel-image projection uses scipy least squares/SVD. GeV units;
all c_i,d_i are dimensionless with the original operator M_V factors retained.
"""
from pathlib import Path
import json
import math
import unittest
import sys
import numpy as np
import scipy.linalg as la
sys.path.insert(0,str(Path(__file__).resolve().parent))
import calibrate_odd_low_energy as C

ROOT=C.ROOT
NAMES=[f'c{i}' for i in range(1,8)]+[f'd{i}' for i in range(1,5)]
CENTRAL=dict(rho_input='rho_charged',kstar_input='kstar_charged',phi_input='phi_evaluation',branch='near_ideal')
HYPOTHESES=('universal','single_trace_positive','single_trace_negative','double_trace')

def read(name):
    return json.loads((ROOT/'data/vector_sector'/name).read_text())

def project_original(rows, embedding, parent_coefficients):
    """Minimum Euclidean norm representative, at fixed even-sector parameters."""
    matrix=rows@embedding
    target=rows@parent_coefficients
    coefficients=la.lstsq(matrix,target,cond=1e-10)[0]
    null=la.null_space(matrix,rcond=1e-10)
    defect=float(la.norm(matrix@coefficients-target)/max(1.,la.norm(target)))
    return coefficients,null,defect

def main(*, write_output=True):
    files=['odd_low_energy_calibration.json','odd_low_energy_profiles.json','odd_low_energy_continuation.json','odd_low_energy_profile_reanchoring.json','odd_low_energy_calibration_verification.json','odd_low_energy_jacobians.json','stage3_results.json','alda_2026nrd_arxiv_v1.json','pdg2024.json']
    archive={name:read(name) for name in files}
    profiles=archive['odd_low_energy_profiles.json']['results']
    overrides={}
    for name in ('odd_low_energy_continuation.json','odd_low_energy_profile_reanchoring.json'):
        for item in archive[name]['results']:
            overrides.setdefault(item['pair_index'],[]).append((name,item['best']))
    calibration=archive['odd_low_energy_calibration.json']['results']
    local={x['fit_index']:x for x in archive['odd_low_energy_calibration_verification.json']['local_prediction_uncertainties']}
    ass=archive['stage3_results.json']['assumptions'];alda=archive['alda_2026nrd_arxiv_v1.json']
    tree=C.J.module('alda_neutral_inputs').alda_2026nrd_neutral_inputs(mass_scheme='tree',F=C.J.F,fpi=C.J.F)
    fast=C.FastChannels(tree,alda['physical_masses_gev'],np.sqrt(4*np.pi*ass['alpha_em']))
    embedding=C.model_maps()['original_11']
    selections=[('central_single_trace_negative', dict(CENTRAL, source_hypothesis='single_trace_negative'))]
    output=[]
    observations=json.loads((ROOT/'data/odd_parity/odd_low_energy_measurements.json').read_text())['observations']
    for label,scenario in selections:
        control=label if label.startswith('without_') else 'primary'
        pair_index,pair=next((i,x) for i,x in enumerate(profiles) if x['scenario']==scenario and x['control']==control)
        ci,base=next((i,x) for i,x in enumerate(calibration) if x['scenario']==scenario and x['control']==control and x['model']=='original_11')
        candidates=[('odd_low_energy_calibration.json',base),('odd_low_energy_profiles.json',pair['best']),*overrides.get(pair_index,[])]
        source,best=min(candidates,key=lambda x:x[1]['chi2'])
        sector,mk=C.J.action_sector(np.array(best['even_parameters']),scenario['source_hypothesis'],ass)
        channels=fast(sector,mk);labels=base['measurement_labels'];rows=np.array([channels[k]['row'] for k in labels])
        coef,null,defect=project_original(rows,embedding,np.array(best['gauge_representative_coefficients']))
        # Width and interference invariance follows the documented epsilon adapter.
        canonical=-coef
        if defect>1e-8 or la.norm(rows@embedding@null)>1e-8:
            raise ValueError('original coefficient projection or nullspace failed')
        if not np.allclose(rows@embedding@coef,-rows@embedding@canonical,rtol=1e-12,atol=1e-12):
            raise ValueError('epsilon coefficient adapter failed')
        even_names=['m_rho_GeV','m_Kstar_GeV','Delta_V_GeV2','kappa_V','F_V_GeV']
        if scenario['source_hypothesis']!='universal':even_names.append('z_V_GeV_inverse')
        if scenario['source_hypothesis']=='double_trace':even_names.append('y_V_GeV_inverse')
        even_values=best['even_parameters']
        if len(even_names)!=len(even_values):raise ValueError('source parameter naming does not match fitted vector')
        parameters = dict(zip(even_names, even_values))
        parameters.update(zip(NAMES, canonical.tolist()))
        parameters.update({f'Gamma_total_{name}_GeV': value
                           for name, value in best['shared_total_widths_gev'].items()})
        light_source, strange_source = 0.036616, 0.943684
        off_diagonal = (math.sqrt(2)*parameters['Delta_V_GeV2']/3
                        + parameters['kappa_V']*(light_source+strange_source)/math.sqrt(2))
        diagonal_difference = (2*parameters['m_Kstar_GeV']**2-2*parameters['m_rho_GeV']**2
                               - parameters['Delta_V_GeV2']/3
                               + 2*parameters['kappa_V']*(strange_source-light_source))
        parameters['omega_phi_mixing_angle_rad'] = .5*math.atan2(2*off_diagonal, diagonal_difference)
        parameters['G_V_GeV'] = read('gv_conditional_fit.json')['parameters']['G_V_GeV']
        output.append(parameters)
    result = dict(metadata=dict(bibtex_references=[
        'Alda:2026nrd', 'Ruiz-Femenia:2003jdx', 'ParticleDataGroup:2024cfk',
        'Ecker:1988te', 'Ecker:1989yg']), parameters=output[0])
    if write_output:
        path=ROOT/'data/vector_sector/odd_preliminary_fit.json'
        path.write_text(json.dumps(result,indent=2)+'\n')
        print(path.relative_to(ROOT))
    return result, null

class ProjectionTests(unittest.TestCase):
    def test_representative_preserves_predictions_and_exposes_ambiguity(self):
        rows=np.array([[1.,2.,0.],[0.,0.,1.]])
        embedding=np.eye(3)
        representative,null,defect=project_original(rows,embedding,np.array([3.,4.,5.]))
        self.assertLess(defect,1e-12)
        self.assertEqual(null.shape,(3,1))
        np.testing.assert_allclose(rows@(representative+7*null[:,0]),[11.,5.],atol=1e-12)
        canonical=-representative
        np.testing.assert_allclose(-rows@canonical,[11.,5.],atol=1e-12)

if __name__=='__main__':main()
