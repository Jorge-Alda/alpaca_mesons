#!/usr/bin/env python3
"""Generate Pelaez:2019eqa P-wave phases for the external KT solver.

Run: .venv/bin/python scripts/prepare_vppp_pipi_phases.py
Tests: .venv/bin/python -m unittest discover -s scripts -p prepare_vppp_pipi_phases.py

Masses/invariants inside the formulas use GeV/GeV^2. Output columns are
s/m_pi^2 and delta_1^1 in radians. The fixed masses below explicitly use an
isospin-average kaon approximation; covariance is not supplied by the source.
"""
from pathlib import Path
import hashlib
import json
import unittest
import tempfile
import numpy as np

PION_MASS = 0.13957039
KAON_MASS = (0.493677 + 0.497611)/2
RHO_MASS = 0.7749
CONFORMAL_SCALE = 1.43**2
CONFORMAL_ALPHA = 0.3
PARAMETERS = {
    "I": {"B": [.96,.09,-.07,.58,1.39], "B_errors": [.01,.03,.08,.19,.38], "K": [.049,-.0051], "K_errors": [.002,.0003]},
    "II": {"B": [.96,.09,-.03,.64,1.1], "B_errors": [.01,.03,.08,.19,.32], "K": [.054,-.0062], "K_errors": [.001,.0001]},
    "III": {"B": [.97,.07,-.003,.67,1.0], "B_errors": [.01,.03,.08,.18,.34], "K": [.045,-.0037], "K_errors": [.002,.0002]},
}


def scattering_matrix(squared_energy, fit="I", *, use_product=False):
    """Return upper-lip physical S=1+2i sigma t, Eq. (20).

    The inelastic exponent is the source's delta_in, so its contribution
    to the physical scattering phase is half its real part.
    """
    squared_energy = np.asarray(squared_energy, dtype=float)
    above_threshold = squared_energy > 4*PION_MASS**2
    safe_energy = np.where(above_threshold, squared_energy, 4*PION_MASS**2 + 1e-15)
    pion_momentum = np.sqrt(safe_energy/4-PION_MASS**2)
    kaon_momentum_squared = safe_energy/4-KAON_MASS**2
    sigma_kaon = np.sqrt(1-4*KAON_MASS**2/safe_energy+0j)
    # At the kaon threshold the logarithmic term has a finite zero limit.
    with np.errstate(divide="ignore", invalid="ignore"):
        loop_function = 2/np.pi + sigma_kaon/np.pi*np.log((sigma_kaon-1)/(sigma_kaon+1))
    loop_function = np.where(sigma_kaon == 0, 2/np.pi+0j, loop_function)
    conformal_variable = (np.sqrt(safe_energy)-CONFORMAL_ALPHA*np.sqrt(CONFORMAL_SCALE-safe_energy))/(np.sqrt(safe_energy)+CONFORMAL_ALPHA*np.sqrt(CONFORMAL_SCALE-safe_energy))
    parameters = PARAMETERS[fit]
    polynomial = np.polynomial.polynomial.polyval(conformal_variable, parameters["B"])
    phi = np.sqrt(safe_energy)/(2*pion_momentum**3)*(RHO_MASS**2-safe_energy)*(2*PION_MASS**3/(RHO_MASS**2*np.sqrt(safe_energy))+polynomial)
    inelastic_exponent = loop_function*(parameters["K"][0]+parameters["K"][1]*safe_energy/KAON_MASS**2)*pion_momentum**3/(np.sqrt(safe_energy)*PION_MASS**2)*kaon_momentum_squared/KAON_MASS**2
    pion_sigma = 2*pion_momentum/np.sqrt(safe_energy)
    conformal_amplitude = 1/(pion_sigma*(phi-1j))
    inelastic_amplitude = np.expm1(1j*inelastic_exponent)/(2j*pion_sigma)
    full_amplitude = (conformal_amplitude + inelastic_amplitude
                      + 2j*pion_sigma*conformal_amplitude*inelastic_amplitude)
    result = 1+2j*pion_sigma*full_amplitude
    if use_product:
        result = (phi+1j)/(phi-1j)*np.exp(1j*inelastic_exponent)
    return np.where(above_threshold, result, 1+0j)


def phase_table(fit="I"):
    """Return dense KT table, stopping at its fixed continuation point 80."""
    pion_units = np.linspace(4,80,3801)
    phase = np.unwrap(np.angle(scattering_matrix(pion_units*PION_MASS**2,fit)))/2
    return np.column_stack((pion_units,phase))


def write_tables(output_directory):
    """Write central phase tables and auditable input/diagnostic manifest."""
    output_directory = Path(output_directory)
    output_directory.mkdir(parents=True,exist_ok=True)
    manifest = {"bibtex_keys":["Pelaez:2019eqa","ParticleDataGroup:2024cfk"], "erratum_url":"https://link.springer.com/article/10.1140/epjc/s10052-024-13589-6", "source":"arXiv:1907.13162v3, Eq. (20)-(24), P-wave coefficient tables tab:para11, tab:para112, tab:para113", "pion_mass_GeV":PION_MASS,"kaon_mass_GeV":KAON_MASS,"rho_mass_GeV":RHO_MASS,"rho_mass_error_GeV":.0012,"mass_assumption":"Charged pion mass and common (0.493677+0.497611)/2 GeV kaon mass (isospin-average approximation).", "conformal_alpha":CONFORMAL_ALPHA,"conformal_scale_GeV2":CONFORMAL_SCALE,"columns":["s/m_pi^2","delta_1^1 [rad]"],"covariance":None,"uncertainty_note":"Marginal parameter errors only; three global fits are alternatives, not covariance samples.","fits":{}}
    for fit, parameters in PARAMETERS.items():
        table = phase_table(fit)
        path = output_directory/f"pipi_pwave_pelaez_2019_{fit}.dat"
        np.savetxt(path,table,fmt="%.12e")
        manifest["fits"][fit] = {**parameters,"filename":path.name,"sha256":hashlib.sha256(path.read_bytes()).hexdigest(),"phase_at_80_rad":float(table[-1,1]),"derivative_at_80_per_pion_unit":float((3*table[-1,1]-4*table[-2,1]+table[-3,1])/.04),"minimum_phase_increment_rad":float(np.min(np.diff(table[:,1])))}
    manifest_path = output_directory/"pipi_pwave_pelaez_2019_manifest.json"
    manifest_path.write_text(json.dumps(manifest,indent=2)+"\n")
    return manifest


class PhaseTests(unittest.TestCase):
    def test_threshold_and_elastic_unitarity(self):
        for fit in PARAMETERS:
            self.assertEqual(scattering_matrix(4*PION_MASS**2,fit),1+0j)
            energies = np.linspace(4*PION_MASS**2+1e-8,4*KAON_MASS**2-1e-8,200)
            np.testing.assert_allclose(abs(scattering_matrix(energies,fit)),1,atol=5e-15)

    def test_physical_phase_and_continuity(self):
        for fit in PARAMETERS:
            table = phase_table(fit)
            self.assertTrue(np.all(np.diff(table[:,1])>0))
            self.assertLess(np.max(np.diff(table[:,1])),.01)
            self.assertGreater(table[-1,1],2)
            self.assertLess(table[-1,1],np.pi)
            matrix = scattering_matrix(table[:,0]*PION_MASS**2,fit)
            product_matrix = scattering_matrix(table[:,0]*PION_MASS**2,fit,use_product=True)
            np.testing.assert_allclose(matrix,product_matrix,atol=3e-15)
            np.testing.assert_allclose(np.exp(2j*table[:,1]),product_matrix/abs(product_matrix),atol=3e-15)
            threshold=4*KAON_MASS**2
            self.assertLess(abs(scattering_matrix(threshold-1e-9,fit)-scattering_matrix(threshold+1e-9,fit)),1e-6)

    def test_two_column_roundtrip_and_units(self):
        with tempfile.TemporaryDirectory() as directory:
            manifest=write_tables(directory)
            for fit in PARAMETERS:
                table=np.loadtxt(Path(directory)/manifest["fits"][fit]["filename"])
                self.assertEqual(table.shape,(3801,2))
                np.testing.assert_allclose(table,phase_table(fit),rtol=5e-13,atol=5e-13)
                self.assertAlmostEqual(table[-1,0]*PION_MASS**2,80*PION_MASS**2)


if __name__ == "__main__":
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-directory",default="data/odd_parity")
    arguments=parser.parse_args()
    results=write_tables(arguments.output_directory)
    for fit,values in results["fits"].items():
        print(f"{fit}: delta(80)={values['phase_at_80_rad']:.10f} rad, derivative={values['derivative_at_80_per_pion_unit']:.10g}")
