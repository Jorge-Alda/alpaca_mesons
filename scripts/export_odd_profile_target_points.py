"""Export full primary profiled nuisance points, with separate coefficient gauges.
Run .venv/bin/python scripts/export_odd_profile_target_points.py.
"""
import json
from pathlib import Path
import numpy as np
import scipy.linalg as la
import calibrate_odd_low_energy as C
ROOT=C.ROOT

def main():
    profiles=json.loads((ROOT/'data/vector_sector/odd_low_energy_profiles.json').read_text())
    verification=json.loads((ROOT/'data/vector_sector/odd_low_energy_profiles_verification.json').read_text())
    assert profiles['complete_256_pair_run'] and verification['complete']
    rows=[r for r in profiles['results'] if r['control']=='primary'];assert len(rows)==64
    assumptions=json.loads((ROOT/'data/vector_sector/stage3_results.json').read_text())['assumptions']
    alda=json.loads((ROOT/'data/vector_sector/alda_arxiv_v1.json').read_text())
    tree=C.J.module('alda_neutral_inputs').alda_neutral_inputs(mass_scheme='tree',F=C.J.F,fpi=C.J.F)
    fast=C.FastChannels(tree,alda['physical_masses_gev'],np.sqrt(4*np.pi*assumptions['alpha_em']))
    maps=C.model_maps();points=[]
    calibration=json.loads((ROOT/'data/vector_sector/odd_low_energy_calibration.json').read_text())
    for row in rows:
        best=row['best'];sector,mk=C.J.action_sector(np.array(best['even_parameters']),row['scenario']['source_hypothesis'],assumptions)
        channels=fast(sector,mk)
        # Training rows have the same order as physical_amplitudes in the saved fit.
        old=next(r for r in calibration['results'] if r['scenario_index']==row['scenario_index'] and r['control']=='primary' and r['model']=='soft_leading')
        measured=np.array([channels[k]['row'] for k in old['measurement_labels']]);leading=np.array(best['gauge_representative_coefficients'])
        original=maps['original_11']@la.lstsq(measured@maps['original_11'],measured@leading,cond=1e-10)[0]
        assert np.linalg.norm(measured@(original-leading))/max(1.,np.linalg.norm(measured@leading))<1e-9
        for model,coefficients in [('original_11',original),('soft_leading',leading)]:
            points.append(dict(scenario_index=row['scenario_index'],model=model,even_parameters=best['even_parameters'],gauge_representative_coefficients=coefficients.tolist()))
    target=ROOT/'data/vector_sector/odd_profile_primary_target_points.json'
    target.write_text(json.dumps(points,indent=2,allow_nan=False)+'\n');print(target.relative_to(ROOT))

if __name__=='__main__':main()
