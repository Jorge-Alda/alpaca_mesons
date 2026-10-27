"""Conditional two-body running widths for supplied VPPP spectral inputs.

Run: .venv/bin/python scripts/build_vppp_partial_widths.py INPUT.npz OUTPUT.npz
Derive: .venv/bin/python scripts/build_vppp_partial_widths.py --derive-formulas
Tests: .venv/bin/python -m unittest discover -s scripts -p build_vppp_partial_widths.py

Natural units GeV. Explicit approximation: each channel's reduced vertex is
constant in s. V->P gamma then scales as k(s)^3; V->P1 P2 scales as
k(s)^3/s, with exact two-body threshold. On-pole partial widths, pole
masses, and daughter masses are caller inputs; no PDG central defaults.
Radiative convention: memory/rcht/odd-parity/radiative-width-kinematics.md,
with sources Ruiz-Femenia:2003jdx, Ecker:1989yg, ParticleDataGroup:2024cfk.
The two-pseudoscalar amplitude is explicitly assumed g epsilon.(p1-p2).
No barrier factor, off-shell vertex correction or three-body approximation
is introduced. Rescattered 3pi widths remain a required EXTERNAL grid.

Input NPZ: spectral_s(n,), resonance_masses(r,), channel_types(c,) strings
'Pgamma' or 'PP', channel_resonance_indices(c,), daughter_masses(c,2),
on_pole_partial_widths(c,), external_three_pion_widths(n,r),
partial_width_source (scalar string), three_pion_width_source (scalar
string), omitted_channels (scalar string; explicitly state 'none' if so).
For Pgamma the second daughter mass must be zero. Other spectral-template
inputs can be included and are copied unchanged. The output resonance_widths
is the SUM OF SUPPLIED CHANNELS, not a certified physical total width.
"""

import argparse
import hashlib
import json
from pathlib import Path
import unittest

import numpy as np
import sympy as sp


def derive_formulas():
    """Derive polarization averages, phase space and normalized running shapes."""
    invariant, momentum, pole_mass, pole_momentum = sp.symbols('s k M k_on', positive=True)
    radiative_vertex, pseudoscalar_vertex = sp.symbols('G g', real=True)
    # Explicit rest-frame polarization contraction. Daughters lie along z;
    # the Levi-Civita contraction is a determinant, squared before averaging.
    vector_polarizations = [sp.eye(4)[:, index] for index in [1, 2, 3]]
    photon_polarizations = vector_polarizations[:2]
    four_momentum = sp.Matrix([sp.sqrt(invariant), 0, 0, 0])
    daughter_energy = sp.Symbol('E_P', real=True)
    daughter_momentum = sp.Matrix([daughter_energy, 0, 0, momentum])
    radiative_average = sp.simplify(sum(
        (radiative_vertex*sp.det(sp.Matrix.hstack(vector, photon,
            four_momentum, daughter_momentum)))**2
        for vector in vector_polarizations for photon in photon_polarizations)/3)
    metric = sp.diag(1, -1, -1, -1)
    energy_difference = sp.Symbol('E_1_minus_E_2', real=True)
    daughter_difference = sp.Matrix([energy_difference, 0, 0, 2*momentum])
    pseudoscalar_average = sp.simplify(sum(
        (pseudoscalar_vertex*(vector.T*metric*daughter_difference)[0])**2
        for vector in vector_polarizations)/3)
    phase_space = momentum/(8*sp.pi*invariant)
    radiative_width = sp.factor(phase_space*radiative_average)
    pseudoscalar_width = sp.factor(phase_space*pseudoscalar_average)
    radiative_ratio = sp.factor(radiative_width/radiative_width.subs({invariant: pole_mass**2, momentum: pole_momentum}))
    pseudoscalar_ratio = sp.factor(pseudoscalar_width/pseudoscalar_width.subs({invariant: pole_mass**2, momentum: pole_momentum}))
    first_mass, second_mass = sp.symbols('m_1 m_2', nonnegative=True)
    energy = (invariant+first_mass**2-second_mass**2)/(2*sp.sqrt(invariant))
    momentum_squared = sp.factor(energy**2-first_mass**2)
    return [('Radiative spin-averaged amplitude squared', radiative_average),
        ('Two-pseudoscalar spin-averaged amplitude squared', pseudoscalar_average),
        ('Two-body phase-space factor', phase_space),
        ('Radiative partial width', radiative_width),
        ('Two-pseudoscalar partial width', pseudoscalar_width),
        ('Radiative width divided by on-pole width', radiative_ratio),
        ('Two-pseudoscalar width divided by on-pole width', pseudoscalar_ratio),
        ('Two-body momentum squared from energy conservation', momentum_squared)]


def two_body_momentum(spectral_s, first_mass, second_mass):
    """Return rest-frame k in GeV, exactly zero at and below threshold."""
    spectral_s = np.asarray(spectral_s, float)
    momentum = np.zeros_like(spectral_s)
    open_channel = spectral_s > (first_mass+second_mass)**2
    allowed_s = spectral_s[open_channel]
    momentum[open_channel] = np.sqrt((allowed_s-(first_mass+second_mass)**2)
        *(allowed_s-(first_mass-second_mass)**2))/(2*np.sqrt(allowed_s))
    return momentum


def normalized_partial_width(spectral_s, pole_mass, first_mass, second_mass,
                             on_pole_width, channel_type):
    """Normalize the stated constant-vertex two-body shape at the pole."""
    spectral_s = np.asarray(spectral_s, float)
    if channel_type not in ['Pgamma', 'PP']:
        raise ValueError('channel_type must be Pgamma or PP')
    if channel_type == 'Pgamma' and second_mass != 0:
        raise ValueError('Pgamma requires a massless second daughter')
    if first_mass < 0 or second_mass < 0 or pole_mass <= first_mass+second_mass:
        raise ValueError('The supplied on-pole decay must be open with nonnegative daughter masses')
    if np.any(spectral_s < 0) or not np.isfinite(spectral_s).all() or on_pole_width < 0:
        raise ValueError('Finite nonnegative s and on-pole width are required')
    pole_momentum = float(two_body_momentum(np.array(pole_mass**2), first_mass, second_mass))
    momentum = two_body_momentum(spectral_s, first_mass, second_mass)
    width = on_pole_width*(momentum/pole_momentum)**3
    if channel_type == 'PP':
        width = np.divide(width*pole_mass**2, spectral_s, out=np.zeros_like(width), where=spectral_s>0)
    return width


def sum_channel_widths(channel_arrays):
    """Sum supplied (channel,s,resonance) grids without claiming completeness."""
    channel_arrays = np.asarray(channel_arrays, float)
    if channel_arrays.ndim != 3 or np.any(channel_arrays < 0) or not np.isfinite(channel_arrays).all():
        raise ValueError('Channel arrays must be finite nonnegative (channel,s,resonance) grids')
    return np.sum(channel_arrays, axis=0)


def build_width_grid(inputs):
    """Combine explicit two-body parameters with the required external 3pi grid."""
    required = ['spectral_s', 'resonance_masses', 'channel_types', 'channel_resonance_indices',
        'daughter_masses', 'on_pole_partial_widths', 'external_three_pion_widths',
        'partial_width_source', 'three_pion_width_source', 'omitted_channels']
    missing = [name for name in required if name not in inputs]
    if missing:
        raise ValueError('Missing explicit width inputs: '+', '.join(missing))
    spectral_s, masses = np.asarray(inputs['spectral_s'], float), np.asarray(inputs['resonance_masses'], float)
    types = np.asarray(inputs['channel_types'])
    indices = np.asarray(inputs['channel_resonance_indices'])
    daughters = np.asarray(inputs['daughter_masses'], float)
    on_pole = np.asarray(inputs['on_pole_partial_widths'], float)
    external = np.asarray(inputs['external_three_pion_widths'], float)
    if spectral_s.ndim != 1 or masses.ndim != 1 or types.ndim != 1 or indices.shape != types.shape or on_pole.shape != types.shape or daughters.shape != (types.size, 2):
        raise ValueError('Inconsistent channel parameter dimensions')
    if not np.issubdtype(indices.dtype, np.integer) or np.any(indices < 0) or np.any(indices >= masses.size):
        raise ValueError('Channel resonance indices must identify supplied resonance masses')
    if external.shape != (spectral_s.size, masses.size):
        raise ValueError('External rescattered 3pi grid must have shape (s,resonance)')
    if not all(np.isfinite(array).all() for array in [spectral_s, masses, daughters, on_pole]):
        raise ValueError('Width parameters must be finite')
    for name in ['partial_width_source', 'three_pion_width_source', 'omitted_channels']:
        value = np.asarray(inputs[name])
        if value.ndim != 0 or value.dtype.kind not in 'US' or not str(value.item()).strip():
            raise ValueError(name+' must be an explicit nonempty scalar string')
    channels = np.zeros((types.size+1, spectral_s.size, masses.size))
    channels[-1] = external
    for index, kind in enumerate(types):
        resonance_index = indices[index]
        channels[index, :, resonance_index] = normalized_partial_width(spectral_s,
            masses[resonance_index], *daughters[index], on_pole[index], str(kind))
    output = {name: np.asarray(inputs[name]) for name in inputs}
    output['supplied_partial_width_channels'] = channels
    output['resonance_widths'] = sum_channel_widths(channels)
    output['width_source'] = np.array('Constant-reduced-vertex two-body channels plus EXTERNAL rescattered 3pi: '
        +str(inputs['partial_width_source'])+'; '+str(inputs['three_pion_width_source'])
        +'; omitted channels: '+str(inputs['omitted_channels']))
    return output


class PartialWidthTests(unittest.TestCase):
    """Nonphysical algebraic fixtures verify thresholds and normalization."""

    def test_pole_and_threshold(self):
        for kind, first, second in [('Pgamma', .2, 0.), ('PP', .2, .3)]:
            threshold = (first+second)**2
            grid = np.array([0., threshold/2, threshold, 1.])
            width = normalized_partial_width(grid, 1., first, second, .01, kind)
            np.testing.assert_array_equal(width[:3], [0., 0., 0.])
            self.assertAlmostEqual(width[-1], .01)

    def test_radiative_and_pp_shapes(self):
        grid = np.array([1., 4.])
        # Massless daughters: k=sqrt(s)/2, hence s^(3/2) and sqrt(s).
        np.testing.assert_allclose(normalized_partial_width(grid, 1., 0., 0., .01, 'Pgamma'), [.01, .08])
        np.testing.assert_allclose(normalized_partial_width(grid, 1., 0., 0., .01, 'PP'), [.01, .02])

    def test_sum_and_missing_three_pion(self):
        np.testing.assert_array_equal(sum_channel_widths(np.array([[[1., 2.]], [[3., 4.]]])), [[4., 6.]])
        with self.assertRaisesRegex(ValueError, 'external_three_pion_widths'):
            build_width_grid({})

    def test_builder_preserves_external_grid_and_channel_assignment(self):
        inputs = dict(spectral_s=np.array([1., 4.]), resonance_masses=np.array([1., 2.]),
            channel_types=np.array(['Pgamma']), channel_resonance_indices=np.array([0]),
            daughter_masses=np.array([[0., 0.]]), on_pole_partial_widths=np.array([.01]),
            external_three_pion_widths=np.array([[.02, .03], [.04, .05]]),
            partial_width_source=np.array('nonphysical algebraic fixture'),
            three_pion_width_source=np.array('nonphysical algebraic fixture'),
            omitted_channels=np.array('all other channels; fixture only'))
        output = build_width_grid(inputs)
        np.testing.assert_allclose(output['resonance_widths'], [[.03, .03], [.12, .05]])
        np.testing.assert_array_equal(output['supplied_partial_width_channels'][-1], inputs['external_three_pion_widths'])

    def test_symbolic_width_normalization(self):
        formulas = dict(derive_formulas())
        invariant, momentum, pole_mass, pole_momentum = sp.symbols('s k M k_on', positive=True)
        self.assertEqual(formulas['Radiative width divided by on-pole width'], momentum**3/pole_momentum**3)
        self.assertEqual(formulas['Two-pseudoscalar width divided by on-pole width'], pole_mass**2*momentum**3/(invariant*pole_momentum**3))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path, nargs='?')
    parser.add_argument('output', type=Path, nargs='?')
    parser.add_argument('--derive-formulas', action='store_true')
    arguments = parser.parse_args()
    if arguments.derive_formulas:
        print('Constant reduced vertices; GeV units; spin average over three initial vector states.')
        for label, expression in derive_formulas():
            print(label+':\n'+sp.latex(expression)+'\n')
        return
    if arguments.input is None or arguments.output is None or arguments.output.suffix != '.npz':
        parser.error('explicit input and output NPZ paths are required')
    if arguments.input.resolve() == arguments.output.resolve():
        parser.error('preserve the input file by selecting a separate output')
    with np.load(arguments.input, allow_pickle=False) as inputs:
        output = build_width_grid(inputs)
    np.savez(arguments.output, **output)
    manifest = {'assumption': 'Constant reduced vertices; sum of supplied channels only.',
        'bibtex_keys': ['Ruiz-Femenia:2003jdx', 'Ecker:1989yg', 'ParticleDataGroup:2024cfk'],
        'input_sha256': hashlib.sha256(arguments.input.read_bytes()).hexdigest(),
        'output_sha256': hashlib.sha256(arguments.output.read_bytes()).hexdigest(),
        'width_source': str(output['width_source'].item()),
        'three_pion_width': 'Required caller grid; no three-body width calculated here.'}
    arguments.output.with_suffix('.provenance.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(arguments.output)
    print(arguments.output.with_suffix('.provenance.json'))


if __name__ == '__main__':
    main()
