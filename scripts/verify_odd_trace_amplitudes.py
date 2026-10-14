"""Focused new trace kernels plus established amplitude/data regressions."""
import importlib
from pathlib import Path
import sys
import types
import unittest

def main():
    root=Path(__file__).resolve().parents[1];sys.path.insert(0,str(root))
    scripts=importlib.import_module('scripts');name='scripts.implementation'
    package=types.ModuleType(name);package.__path__=[str(root/'scripts/implementation')]
    sys.modules[name]=package;scripts.implementation=package
    modules=('test_odd_trace_amplitudes','test_vector_sector','test_vector_sector_data',
             'test_odd_parity_data','test_pseudoscalar_tff_verification','test_radiative_sm_verification')
    suite=unittest.TestSuite(unittest.defaultTestLoader.loadTestsFromModule(importlib.import_module(name+'.'+m)) for m in modules)
    result=unittest.TextTestRunner(verbosity=1).run(suite)
    if not result.wasSuccessful():raise SystemExit(1)

if __name__=='__main__':main()
