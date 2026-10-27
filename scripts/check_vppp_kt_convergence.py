"""Check external KT numerical convergence; no physical pole widths are assumed.

Run .venv/bin/python scripts/check_vppp_kt_convergence.py OUTPUT_DIRECTORY.
Large grids remain in the supplied ignored external directory. Units are GeV,
with solver invariants divided by m_pi**2. Each solve is reused for all
iteration counts and Dalitz quadratures. Internal solver quadrature is fixed
by upstream; this checks the external phase-space quadrature only.
"""
import argparse
import hashlib
import importlib
import json
from pathlib import Path
import subprocess
import sys
import time
import unittest
import numpy as np
from export_vppp_kt_grid import dalitz_quadrature
from build_vppp_three_pion_width import integrate_unit_scalar

PION_MASS = .13957039
ENERGIES = [.7, .78266, .9, 1.019461, 1.1]

def integrate_solution(solver, module, invariant, iteration, pair_count, angle_count):
    """Return the unit-scalar integral, normalization and reflection error."""
    pair, weights, angles, angle_weights, pair_t, pair_u = dalitz_quadrature(
        invariant, PION_MASS, pair_count, angle_count)
    def evaluate(values):
        return solver(np.asarray(values)/PION_MASS**2, 1,
            module.SubtractionConstant.a1, module.Setting.above, iteration)
    scalar = evaluate(pair)[:, None]+evaluate(pair_t)+evaluate(pair_u)
    integral = integrate_unit_scalar([invariant], PION_MASS, pair[None],
        weights[None], angles, angle_weights, scalar[None, :, :, None])[0, 0]
    zero = complex(evaluate(0.))
    return dict(integral=float(integral), basis_zero_real=zero.real,
        basis_zero_imag=zero.imag, finite=bool(np.isfinite(scalar).all()),
        reflection_relative=float(np.max(np.abs(scalar-scalar[:, ::-1]))/np.max(np.abs(scalar))))

def run_checks(output_directory):
    """Solve bounded phase, cutoff, iteration and phase-space quadrature scans."""
    build_directory = Path('.external/khuri_treiman_solver/build').resolve()
    sys.path.insert(0, str(build_directory))
    module = importlib.import_module('_khuri_pinocchio')
    output_directory.mkdir(parents=True, exist_ok=True)
    records = []
    settings = [('I', 2.), ('I', 2.5), ('I', 3.), ('II', 2.5), ('III', 2.5)]
    for phase, cutoff_energy in settings:
        phase_file = Path(f'data/odd_parity/pipi_pwave_pelaez_2019_{phase}.dat').resolve()
        for energy in ENERGIES:
            start = time.monotonic()
            stem = output_directory/f'{phase}_{cutoff_energy}_{energy}'
            solver = module.IterationV3Pi(1., energy/PION_MASS, str(phase_file),
                .1, .05, .3, .2, (cutoff_energy/PION_MASS)**2, 6, 1,
                str(stem)+'_a.dat', str(stem)+'_b.dat')
            combinations = [(iteration, 24, 16) for iteration in [1, 2, 3, 4, 5, 6]]
            if phase == 'I' and cutoff_energy == 2.5:
                combinations += [(6, pair, angle) for pair, angle in [(12, 8), (48, 32), (96, 64)]]
            for iteration, pair_count, angle_count in combinations:
                result = integrate_solution(solver, module, energy**2, iteration, pair_count, angle_count)
                records.append(dict(phase=phase, cutoff_energy_GeV=cutoff_energy,
                    energy_GeV=energy, iteration=iteration, pair_nodes=pair_count,
                    angle_nodes=angle_count, **result))
            (output_directory/'metrics.json').write_text(json.dumps(dict(records=records), indent=2)+'\n')
            print(f'{phase} {cutoff_energy} {energy}: {time.monotonic()-start:.1f}s', flush=True)
    manifest = dict(pion_mass_GeV=PION_MASS,
        phase_sha256={phase: hashlib.sha256(Path(f'data/odd_parity/pipi_pwave_pelaez_2019_{phase}.dat').read_bytes()).hexdigest() for phase in ['I', 'II', 'III']},
        solver_revision=subprocess.check_output(['git', '-C', str(build_directory.parent), 'rev-parse', 'HEAD'], text=True).strip(),
        module_sha256=hashlib.sha256(Path(module.__file__).read_bytes()).hexdigest(),
        convention='exact isospin, unit a1, full scalar sum of three channels; no on-pole widths',
        internal_quadrature='upstream fixed; not varied', records=records)
    (output_directory/'metrics.json').write_text(json.dumps(manifest, indent=2)+'\n')
    (output_directory/'summary.json').write_text(json.dumps(summarize(records), indent=2)+'\n')

def compare_phase_points():
    """Interpolate the prepared phase tables at the Roy matching energies."""
    energies = np.asarray([.8, 1.15])
    results = {}
    for phase in ['I', 'II', 'III']:
        table = np.loadtxt(f'data/odd_parity/pipi_pwave_pelaez_2019_{phase}.dat')
        results[phase] = np.rad2deg(np.interp((energies/PION_MASS)**2,
            table[:, 0], table[:, 1])).tolist()
    return dict(energies_GeV=energies.tolist(), madrid_degrees=results,
        Colangelo_2018mtw_fitted_central_degrees=[110.4, 165.7],
        Colangelo_2018mtw_reference_estimates_degrees=[108.9, 166.5])

def summarize(records):
    """Compare common quadratures and normalize width shapes to each pole."""
    def select(phase='I', cutoff=2.5, iteration=6, pair=24, angle=16):
        return [record for record in records if record['phase'] == phase
            and record['cutoff_energy_GeV'] == cutoff and record['iteration'] == iteration
            and record['pair_nodes'] == pair and record['angle_nodes'] == angle]
    def relative(first, second):
        return [a['integral']/b['integral']-1 for a, b in zip(first, second)]
    baseline = select()
    shape_sets = {}
    for phase, cutoff in [('I', 2.), ('I', 2.5), ('I', 3.), ('II', 2.5), ('III', 2.5)]:
        selected = select(phase, cutoff)
        integrals = np.asarray([record['integral'] for record in selected])
        energies = np.asarray([record['energy_GeV'] for record in selected])
        shape_sets[f'{phase}_{cutoff}'] = {
            'omega_normalized': ((energies/.78266)**3*integrals/integrals[1]).tolist(),
            'phi_normalized': ((energies/1.019461)**3*integrals/integrals[3]).tolist()}
    return dict(phase_matching_points=compare_phase_points(), energies_GeV=ENERGIES,
        baseline_integrals=[record['integral'] for record in baseline],
        iteration_relative_to_six={str(iteration): relative(select(iteration=iteration), baseline)
            for iteration in [1, 2, 3, 4, 5]},
        quadrature_relative_to_96x64={f'{pair}x{angle}': relative(select(pair=pair, angle=angle),
            select(pair=96, angle=64)) for pair, angle in [(12, 8), (24, 16), (48, 32)]},
        cutoff_relative_to_2p5={str(cutoff): relative(select(cutoff=cutoff), baseline) for cutoff in [2., 3.]},
        phase_relative_to_I={phase: relative(select(phase=phase), baseline) for phase in ['II', 'III']},
        normalized_width_shapes=shape_sets,
        max_zero_normalization_error=max(abs(complex(record['basis_zero_real'], record['basis_zero_imag'])-1) for record in records),
        max_reflection_relative=max(record['reflection_relative'] for record in records),
        all_finite=all(record['finite'] for record in records))

def check_amplitude_iterations(output_directory):
    """Compare full complex scalars through eight iterations at both poles."""
    build_directory = Path('.external/khuri_treiman_solver/build').resolve()
    sys.path.insert(0, str(build_directory))
    module = importlib.import_module('_khuri_pinocchio')
    output_directory.mkdir(parents=True, exist_ok=True)
    records = []
    for energy in [.78266, 1.019461]:
        phase_file = Path('data/odd_parity/pipi_pwave_pelaez_2019_I.dat').resolve()
        stem = output_directory/f'amplitude_{energy}'
        solver = module.IterationV3Pi(1., energy/PION_MASS, str(phase_file),
            .1, .05, .3, .2, (2.5/PION_MASS)**2, 8, 1,
            str(stem)+'_a.dat', str(stem)+'_b.dat')
        pair, weights, angles, angle_weights, pair_t, pair_u = dalitz_quadrature(
            energy**2, PION_MASS, 48, 32)
        scalars = []
        for iteration in range(1, 9):
            def evaluate(values):
                return solver(np.asarray(values)/PION_MASS**2, 1,
                    module.SubtractionConstant.a1, module.Setting.above, iteration)
            scalars.append(evaluate(pair)[:, None]+evaluate(pair_t)+evaluate(pair_u))
        last = scalars[-1]
        for iteration, scalar in enumerate(scalars, 1):
            records.append(dict(energy_GeV=energy, iteration=iteration,
                maximum_relative_to_final_maximum=float(np.max(np.abs(scalar-last))/np.max(np.abs(last))),
                relative_l2=float(np.linalg.norm(scalar-last)/np.linalg.norm(last)),
                max_pointwise_relative=float(np.max(np.abs((scalar-last)/last)))))
        np.savez(output_directory/f'amplitude_iterations_{energy}.npz',
            full_scalar=np.asarray(scalars), pair_s=pair, z=angles)
        (output_directory/'amplitude_iterations.json').write_text(json.dumps(records, indent=2)+'\n')

class ConvergenceSummaryTests(unittest.TestCase):
    """Check pole normalization using an exactly known energy dependence."""
    def test_normalized_width_and_relative_changes(self):
        records = []
        for phase, cutoff in [('I', 2.), ('I', 2.5), ('I', 3.), ('II', 2.5), ('III', 2.5)]:
            for iteration in range(1, 7):
                for pair, angle in [(12, 8), (24, 16), (48, 32), (96, 64)]:
                    for energy in ENERGIES:
                        records.append(dict(phase=phase, cutoff_energy_GeV=cutoff,
                            iteration=iteration, pair_nodes=pair, angle_nodes=angle,
                            energy_GeV=energy, integral=energy**2,
                            basis_zero_real=1., basis_zero_imag=0.,
                            finite=True, reflection_relative=0.))
        summary = summarize(records)
        np.testing.assert_allclose(summary['normalized_width_shapes']['I_2.5']['omega_normalized'],
            (np.asarray(ENERGIES)/.78266)**5)
        self.assertEqual(summary['normalized_width_shapes']['I_2.5']['omega_normalized'][1], 1.)
        self.assertEqual(summary['normalized_width_shapes']['I_2.5']['phi_normalized'][3], 1.)
        np.testing.assert_array_equal(summary['iteration_relative_to_six']['5'], np.zeros(5))
        self.assertTrue(summary['all_finite'])

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output_directory', type=Path)
    parser.add_argument('--mode', choices=['scan', 'amplitude'], default='scan')
    arguments = parser.parse_args()
    if arguments.mode == 'amplitude':
        check_amplitude_iterations(arguments.output_directory)
    else:
        run_checks(arguments.output_directory)
