"""Reproduce arXiv-v1 neutral generator inputs and verify source values.
Run .venv/bin/python scripts/verify_alda_neutral_inputs.py.
"""
import importlib.util
import json
from pathlib import Path
import numpy as np


def main():
    path=Path(__file__).parent/'implementation/alda_neutral_inputs.py'
    spec=importlib.util.spec_from_file_location('alda_inputs',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    data=json.loads(module.SNAPSHOT.read_text())
    inputs=module.alda_neutral_inputs(mass_scheme='physical')
    matrices=[x['TP'] for x in inputs.values()]
    for name,x in inputs.items():
        np.testing.assert_allclose(np.sqrt(2)*np.diag(x['TP']),data['published_lambda_diagonals'][name],atol=3e-5,rtol=0)
        assert x['mp']==data['physical_masses_gev'][name]
    np.testing.assert_allclose([[np.trace(a@b) for b in matrices] for a in matrices],np.eye(3),atol=1e-14)
    scaled=module.alda_neutral_inputs(mass_scheme='tree',F=.1)
    for name,x in inputs.items():
        np.testing.assert_allclose(scaled[name]['TP'],x['TP']*.1/.093)
        assert scaled[name]['mp']==scaled[name]['tree_mass']
    assert abs(inputs['eta']['tree_mass']-.49410)<3e-5
    print(path)


if __name__=='__main__':main()
