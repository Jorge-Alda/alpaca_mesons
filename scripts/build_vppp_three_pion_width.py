"""Integrate an EXTERNAL unit-subtraction full KT scalar into a 3pi width shape.

Run: .venv/bin/python scripts/build_vppp_three_pion_width.py INPUT.npz OUTPUT.npz
Derive: .venv/bin/python scripts/build_vppp_three_pion_width.py --derive-formulas
Tests: .venv/bin/python -m unittest discover -s scripts -p build_vppp_three_pion_width.py

Hoferichter:2019mqg Eq. (5), repeated by Hoferichter:2025lcz Eq. (9).
Exact isospin; pion_mass in GeV. No KT solution, interpolation, phase model,
resonance line shape or subtraction function a(q^2) is generated here.
The caller supplies the FULL scalar F(s,t,u;q^2)/a(q^2), including all three
single-variable terms and crossed-channel rescattering, on the quadrature
points. In this unit-subtraction convention it is dimensionless. Inputs
must exclude the photon normalization a(q^2) and its resonance poles.
The fixed vector-current coupling is the explicit normalization assumption:
Gamma(s)/Gamma(M^2)=(s/M^2)^(3/2)*sigma_unit(s)/sigma_unit(M^2).
Only the ratio is used: alpha^2 cancels. sigma_unit/alpha^2 has units GeV^4
for dimensionless unit-subtraction scalar and is not a measured cross section.

Required NPZ: spectral_s(nq,) positive increasing GeV^2;
resonance_masses(r,) GeV; on_pole_three_pion_widths(r,) GeV;
pion_mass scalar GeV; pair_s(nq,ns) GeV^2;
pair_s_weights(nq,ns) GeV^2; z(nz,), z_weights(nz,);
full_scalar_unit(nq,ns,nz,r) complex;
kt_grid_source and normalization_assumptions scalar nonempty strings.
Every pole M^2 must occur explicitly in spectral_s: no pole interpolation.
For open q^2, pair_s lies in [4m_pi^2,(sqrt(q^2)-m_pi)^2].
At closed q^2 the integration entries are ignored and width is exactly zero.
Other fields are copied, including downstream spectral quadrature/tail inputs.
Output external_three_pion_widths(nq,r) feeds build_vppp_partial_widths.py.
Grid convergence remains the caller's responsibility; repeat with refined
EXTERNAL KT/quadrature grids. No synthetic physical amplitudes are supplied.
"""

import argparse
import hashlib
import json
from pathlib import Path
import unittest

import numpy as np
import sympy as sp


def derive_formulas():
    """Derive spectral-tensor width/cross-section ratio and angular Jacobian."""
    invariant, charge, current_coupling, spectral_density = sp.symbols('q_squared e g_current rho', positive=True)
    # Evaluate the spectral-tensor contractions explicitly in the q rest frame.
    metric = sp.diag(1, -1, -1, -1)
    four_momentum = sp.Matrix([sp.sqrt(invariant), 0, 0, 0])
    electron = sp.Matrix([sp.sqrt(invariant)/2, 0, 0, sp.sqrt(invariant)/2])
    positron = sp.Matrix([sp.sqrt(invariant)/2, 0, 0, -sp.sqrt(invariant)/2])
    current_tensor = (four_momentum*four_momentum.T-invariant*metric)*spectral_density
    lepton_tensor = (electron*positron.T+positron*electron.T
        -metric*(electron.T*metric*positron)[0])
    lepton_contraction = sp.simplify(sp.trace(lepton_tensor*metric*current_tensor*metric))
    cross_section = sp.simplify(charge**4/(2*invariant**3)*lepton_contraction)
    # Three-vector-polarization average gives s*rho; decay flux 1/(2sqrt(s)).
    vector_average = (-metric+four_momentum*four_momentum.T/invariant)/3
    vector_contraction = sp.simplify(sp.trace(vector_average*metric*current_tensor*metric))
    width = sp.simplify(current_coupling**2/(2*sp.sqrt(invariant))*vector_contraction)
    ratio = sp.simplify(width/cross_section)
    pole_mass = sp.Symbol('M', positive=True)
    sigma, sigma_pole = sp.symbols('sigma_unit sigma_unit_on', positive=True)
    normalized = sp.simplify(ratio/ratio.subs(invariant, pole_mass**2)*sigma/sigma_pole)
    pair_s, pion_mass, angle, kappa = sp.symbols('pair_s m_pi z kappa', real=True)
    pair_t = (3*pion_mass**2+invariant-pair_s+kappa*angle)/2
    pair_u = (3*pion_mass**2+invariant-pair_s-kappa*angle)/2
    return [('Integrated current tensor contracted with lepton tensor', lepton_contraction),
        ('Cross section from current spectral tensor', cross_section),
        ('Width for a fixed vector-current coupling', width),
        ('Width divided by cross section', ratio),
        ('Normalized running three-pion width', normalized),
        ('Pair invariant t', pair_t), ('Pair invariant u', pair_u),
        ('Jacobian dt/dz', sp.diff(pair_t, angle))]


def integrate_unit_scalar(spectral_s, pion_mass, pair_s, pair_weights, z, z_weights, scalar):
    """Integrate the cited cross-section formula divided by alpha^2."""
    spectral_s, pair_s = np.asarray(spectral_s, float), np.asarray(pair_s, float)
    pair_weights, z, z_weights = np.asarray(pair_weights, float), np.asarray(z, float), np.asarray(z_weights, float)
    scalar = np.asarray(scalar, complex)
    if spectral_s.ndim != 1 or np.any(spectral_s <= 0) or np.any(np.diff(spectral_s) <= 0):
        raise ValueError('spectral_s must be positive and increasing')
    if pair_s.ndim != 2 or pair_s.shape[0] != spectral_s.size or pair_weights.shape != pair_s.shape:
        raise ValueError('pair_s and pair_s_weights must have matching (nq,ns) dimensions')
    if z.ndim != 1 or z_weights.shape != z.shape or np.any(np.abs(z)>1):
        raise ValueError('Angular nodes and weights must have matching shape inside [-1,1]')
    if scalar.ndim != 4 or scalar.shape[:3] != (spectral_s.size, pair_s.shape[1], z.size):
        raise ValueError('Full scalar grid must have shape (nq,ns,nz,r)')
    if np.any(pair_weights < 0) or np.any(z_weights < 0) or pion_mass <= 0:
        raise ValueError('Quadrature weights must be nonnegative and pion mass positive')
    if not all(np.isfinite(array).all() for array in [spectral_s, pair_s, pair_weights, z, z_weights, scalar]):
        raise ValueError('KT quadrature inputs must be finite')
    result = np.zeros((spectral_s.size, scalar.shape[-1]))
    for index, invariant in enumerate(spectral_s):
        if invariant <= 9*pion_mass**2:
            continue
        pair = pair_s[index]
        if np.any(pair < 4*pion_mass**2-1e-12) or np.any(pair > (np.sqrt(invariant)-pion_mass)**2+1e-12):
            raise ValueError('Pair-invariant quadrature nodes outside the physical Dalitz region')
        triangle = np.maximum((invariant-pair-pion_mass**2)**2-4*pair*pion_mass**2, 0.)
        kappa = np.sqrt(np.maximum(1-4*pion_mass**2/pair, 0.)*triangle)
        # The supplied scalar is evaluated at t,u defined in the module contract;
        # integrate dt=(kappa/2) dz and sin^2(theta)=1-z^2.
        measure = pair_weights[index, :, None]*z_weights[None, :]*kappa[:, None]/2
        kinematics = (pair[:, None]-4*pion_mass**2)*triangle[:, None]*(1-z[None, :]**2)/(768*np.pi*invariant**3)
        result[index] = np.sum((measure*kinematics)[:, :, None]*np.abs(scalar[index])**2, axis=(0, 1))
    return result


def build_three_pion_width(inputs):
    """Normalize externally rescattered scalar integrals to explicit pole widths."""
    names = ['spectral_s', 'resonance_masses', 'on_pole_three_pion_widths', 'pion_mass',
        'pair_s', 'pair_s_weights', 'z', 'z_weights', 'full_scalar_unit', 'kt_grid_source', 'normalization_assumptions']
    missing = [name for name in names if name not in inputs]
    if missing:
        raise ValueError('Missing external KT/width inputs: '+', '.join(missing))
    for name in ['kt_grid_source', 'normalization_assumptions']:
        value = np.asarray(inputs[name])
        if value.ndim != 0 or value.dtype.kind not in 'US' or not str(value.item()).strip():
            raise ValueError(name+' must be an explicit nonempty scalar string')
    spectral_s = np.asarray(inputs['spectral_s'], float)
    masses, pole_widths = np.asarray(inputs['resonance_masses'], float), np.asarray(inputs['on_pole_three_pion_widths'], float)
    pion_mass = float(np.asarray(inputs['pion_mass']))
    integrated = integrate_unit_scalar(spectral_s, pion_mass, inputs['pair_s'], inputs['pair_s_weights'],
        inputs['z'], inputs['z_weights'], inputs['full_scalar_unit'])
    if masses.ndim != 1 or pole_widths.shape != masses.shape or integrated.shape[1] != masses.size or np.any(pole_widths < 0) or not np.isfinite(pole_widths).all():
        raise ValueError('Mass, pole-width and scalar resonance dimensions do not agree')
    widths = np.zeros_like(integrated)
    for index, mass in enumerate(masses):
        matches = np.flatnonzero(np.isclose(spectral_s, mass**2, rtol=1e-12, atol=1e-14))
        if matches.size != 1 or integrated[matches[0], index] <= 0:
            raise ValueError('Each open on-pole mass squared must occur once with a positive KT integral')
        widths[:, index] = pole_widths[index]*(spectral_s/mass**2)**1.5*integrated[:, index]/integrated[matches[0], index]
    output = {name: np.asarray(inputs[name]) for name in inputs}
    output.update(unit_cross_section_divided_by_alpha_squared=integrated,
        external_three_pion_widths=widths,
        three_pion_width_source=np.array('External full unit-subtraction KT scalar: '+str(inputs['kt_grid_source'])
            +'; fixed vector-current coupling; '+str(inputs['normalization_assumptions'])))
    return output


class ThreePionWidthTests(unittest.TestCase):
    """Algebraic fixtures test the integral without supplying physical KT amplitudes."""

    def test_spectral_tensor_power(self):
        expressions = dict(derive_formulas())
        invariant, mass = sp.symbols('q_squared M', positive=True)
        sigma, sigma_pole = sp.symbols('sigma_unit sigma_unit_on', positive=True)
        self.assertEqual(sp.simplify(expressions['Normalized running three-pion width']
            -(invariant/mass**2)**sp.Rational(3,2)*sigma/sigma_pole), 0)

    def test_zero_scalar_has_zero_integral(self):
        result = integrate_unit_scalar(np.array([1.]), .1, np.array([[.2, .5]]),
            np.array([[.1, .1]]), np.array([-.5, .5]), np.ones(2), np.zeros((1, 2, 2, 1), complex))
        np.testing.assert_array_equal(result, [[0.]])

    def test_missing_external_amplitude_rejected(self):
        with self.assertRaisesRegex(ValueError, 'full_scalar_unit'):
            build_three_pion_width({})

    def test_closed_threshold_zero_and_scalar_rescaling(self):
        # Arbitrary algebraic array; it is not a KT solution or a physical amplitude.
        inputs = dict(spectral_s=np.array([.09, 1.]), resonance_masses=np.array([1.]),
            on_pole_three_pion_widths=np.array([.01]), pion_mass=np.array(.1),
            pair_s=np.array([[.04, .04], [.2, .5]]), pair_s_weights=np.full((2, 2), .1),
            z=np.array([-.5, .5]), z_weights=np.ones(2), full_scalar_unit=np.ones((2, 2, 2, 1), complex),
            kt_grid_source=np.array('nonphysical algebraic fixture'),
            normalization_assumptions=np.array('nonphysical algebraic fixture'))
        first = build_three_pion_width(inputs)
        self.assertAlmostEqual(first['external_three_pion_widths'][0, 0], 0.)
        self.assertAlmostEqual(first['external_three_pion_widths'][1, 0], .01)
        inputs['full_scalar_unit'] *= 3j
        second = build_three_pion_width(inputs)
        np.testing.assert_allclose(first['external_three_pion_widths'], second['external_three_pion_widths'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', nargs='?', type=Path)
    parser.add_argument('output', nargs='?', type=Path)
    parser.add_argument('--derive-formulas', action='store_true')
    arguments = parser.parse_args()
    if arguments.derive_formulas:
        print('Exact isospin; fixed vector-current coupling; unit-subtraction full KT scalar; natural units GeV.')
        for label, expression in derive_formulas():
            print(label+':\n'+sp.latex(expression)+'\n')
        return
    if arguments.input is None or arguments.output is None or arguments.output.suffix != '.npz':
        parser.error('explicit input/output NPZ paths required')
    if arguments.input.resolve() == arguments.output.resolve():
        parser.error('preserve original amplitude inputs with a separate output path')
    with np.load(arguments.input, allow_pickle=False) as inputs:
        output = build_three_pion_width(inputs)
    np.savez(arguments.output, **output)
    manifest = {'bibtex_keys': ['Hoferichter:2019mqg', 'Hoferichter:2025lcz'],
        'input_sha256': hashlib.sha256(arguments.input.read_bytes()).hexdigest(),
        'output_sha256': hashlib.sha256(arguments.output.read_bytes()).hexdigest(),
        'source': str(output['three_pion_width_source']),
        'scope': 'External KT grid integrated only; no KT solution or physical amplitude generated.'}
    arguments.output.with_suffix('.provenance.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(arguments.output)
    print(arguments.output.with_suffix('.provenance.json'))


if __name__ == '__main__':
    main()
