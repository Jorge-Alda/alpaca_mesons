"""Export a separately installed KT solver's one-subtraction basis to a width grid.

Run: .venv/bin/python scripts/export_vppp_kt_grid.py CONFIG.json OUTPUT.npz
The JSON supplies build_directory, phase_file, phase_source, pion_mass [GeV],
spectral_s [GeV^2], resonance_masses [GeV], on_pole_three_pion_widths [GeV],
iterations, cutoff [m_pi^2], pair_nodes, angle_nodes. The phase file must use
s/m_pi^2 and radians. No phase parametrization or physical widths are chosen.
The external GPL solver remains an optional, separately installed dependency.
"""

import argparse
import hashlib
import importlib
import json
from pathlib import Path
import sys
import subprocess
import unittest

import numpy as np


def dalitz_quadrature(invariant, pion_mass, pair_count, angle_count):
    """Return Gauss nodes, weights and t,u in GeV units for an open channel."""
    pair_nodes, pair_weights = np.polynomial.legendre.leggauss(pair_count)
    angles, angle_weights = np.polynomial.legendre.leggauss(angle_count)
    lower = 4*pion_mass**2
    upper = (np.sqrt(invariant)-pion_mass)**2
    pair = (upper+lower)/2+(upper-lower)*pair_nodes/2
    pair_weights *= (upper-lower)/2
    triangle = (invariant-pair-pion_mass**2)**2-4*pair*pion_mass**2
    kappa = np.sqrt(np.maximum((1-lower/pair)*triangle, 0))
    crossed_mean = (invariant+3*pion_mass**2-pair[:, None])/2
    pair_t = crossed_mean+kappa[:, None]*angles/2
    pair_u = crossed_mean-kappa[:, None]*angles/2
    return pair, pair_weights, angles, angle_weights, pair_t, pair_u


def export_grid(configuration):
    """Evaluate the external basis at dimensionless pair invariants."""
    build_directory = Path(configuration['build_directory']).resolve()
    phase_file = Path(configuration['phase_file']).resolve()
    revision = subprocess.check_output(['git', '-C', str(build_directory.parent), 'rev-parse', 'HEAD'], text=True).strip()
    phase = np.loadtxt(phase_file)
    if phase.ndim != 2 or phase.shape[1] != 2 or phase[0, 0] > 4 or phase[-1, 0] < 80:
        raise ValueError('Phase table must contain two columns and cover s/m_pi^2=4 through 80')
    sys.path.insert(0, str(build_directory))
    solver_module = importlib.import_module('_khuri_pinocchio')
    pion_mass = configuration['pion_mass']
    invariants = np.asarray(configuration['spectral_s'], float)
    masses = np.asarray(configuration['resonance_masses'], float)
    pair_count, angle_count = configuration['pair_nodes'], configuration['angle_nodes']
    pair_grid = np.zeros((len(invariants), pair_count))
    weights_grid = np.zeros_like(pair_grid)
    full_scalar = np.zeros((len(invariants), pair_count, angle_count, len(masses)), complex)
    normalization_values = []
    output_directory = build_directory/'alpaca_diagnostics'
    output_directory.mkdir(exist_ok=True)
    for index, invariant in enumerate(invariants):
        if invariant <= 9*pion_mass**2:
            continue
        solver = solver_module.IterationV3Pi(1., np.sqrt(invariant)/pion_mass,
            str(phase_file), .1, .05, .3, .2, configuration['cutoff'],
            configuration['iterations'], 1, str(output_directory/f'a_{index}.dat'),
            str(output_directory/f'b_{index}.dat'))
        def evaluate(pair):
            return solver(np.asarray(pair)/pion_mass**2, 1,
                solver_module.SubtractionConstant.a1, solver_module.Setting.above,
                configuration['iterations'])
        pair, weights, angles, angle_weights, pair_t, pair_u = dalitz_quadrature(
            invariant, pion_mass, pair_count, angle_count)
        pair_grid[index], weights_grid[index] = pair, weights
        scalar = evaluate(pair)[:, None]+evaluate(pair_t)+evaluate(pair_u)
        if not np.isfinite(scalar).all():
            raise ValueError('External KT basis contains nonfinite values')
        full_scalar[index] = scalar[:, :, None]
        normalization_values.append(complex(evaluate(0.)))
    angles, angle_weights = np.polynomial.legendre.leggauss(angle_count)
    return dict(spectral_s=invariants, resonance_masses=masses,
        on_pole_three_pion_widths=np.asarray(configuration['on_pole_three_pion_widths']),
        pion_mass=np.array(pion_mass), pair_s=pair_grid, pair_s_weights=weights_grid,
        z=angles, z_weights=angle_weights, full_scalar_unit=full_scalar,
        kt_grid_source=np.array(configuration['phase_source']+'; external _khuri_pinocchio; '+
            'phase SHA256 '+hashlib.sha256(phase_file.read_bytes()).hexdigest()),
        normalization_assumptions=np.array('Single-variable unit a1; full scalar is F(s)+F(t)+F(u); '+
            'exact isospin; no production a(q^2); fixed vector-current coupling'),
        external_basis_at_zero=np.asarray(normalization_values),
        external_solver_revision=np.array(revision),
        external_module_sha256=np.array(hashlib.sha256(Path(solver_module.__file__).read_bytes()).hexdigest()),
        configuration_json=np.array(json.dumps(configuration, sort_keys=True)))


class DalitzQuadratureTests(unittest.TestCase):
    """Check conservation and quadrature against exact integrals."""

    def test_conservation_and_measure(self):
        pair, weights, angles, angle_weights, pair_t, pair_u = dalitz_quadrature(.64, .14, 8, 6)
        np.testing.assert_allclose(pair[:, None]+pair_t+pair_u, .64+3*.14**2)
        self.assertAlmostEqual(weights.sum(), (.8-.14)**2-4*.14**2)
        self.assertAlmostEqual(np.dot(angle_weights, 1-angles**2), 4/3)
        np.testing.assert_allclose(pair_t, pair_u[:, ::-1])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('configuration', type=Path)
    parser.add_argument('output', type=Path)
    arguments = parser.parse_args()
    np.savez(arguments.output, **export_grid(json.loads(arguments.configuration.read_text())))
