"""Run focused vector and existing amplitude/data regressions.

Run .venv/bin/python scripts/verify_vector_stage4.py. This loads only the
implementation package path, avoiding its unrelated eager alpaca imports.
"""
import importlib
from pathlib import Path
import sys
import types
import unittest


def main():
    sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
    scripts=importlib.import_module('scripts')
    name='scripts.implementation'
    package=types.ModuleType(name)
    package.__path__=[str(Path(__file__).parent/'implementation')]
    sys.modules[name]=package
    scripts.implementation=package
    modules=('test_vector_sector','test_vector_sector_data','test_odd_parity_data',
             'test_pseudoscalar_tff_verification','test_radiative_sm_verification')
    suite=unittest.TestSuite(unittest.defaultTestLoader.loadTestsFromModule(
        importlib.import_module(name+'.'+module)) for module in modules)
    result=unittest.TextTestRunner(verbosity=1).run(suite)
    if not result.wasSuccessful():raise SystemExit(1)


if __name__=='__main__':main()
