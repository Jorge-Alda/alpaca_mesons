"""Generate VPPP resonance/conformal spectral templates from explicit inputs.

Run: .venv/bin/python scripts/build_vppp_spectral_templates.py INPUT.npz OUTPUT.npz
Derive: .venv/bin/python scripts/build_vppp_spectral_templates.py --derive-formulas
Tests: .venv/bin/python -m unittest discover -s scripts -p build_vppp_spectral_templates.py

Hoferichter:2019mqg Eqs. (15),(16), repeated by Hoferichter:2025lcz:
Im B_V(s)=Im[1/(M_V^2-s-i sqrt(s) Gamma_V(s))] per unit residue;
Im C_n(s)=Im[z(s+i0)^n-z(0)^n] per unit conformal coefficient.
Natural units GeV: spectral_s/weights in GeV^2, resonance_masses/widths in
GeV. Resonance templates have units GeV^-2; conformal templates are
dimensionless, multiplying coefficients with units GeV^-3. This is the
real-residue resonance part; rho-omega mixing and complex residues require
additional explicitly calculated spectral functions and are not generated.

Required NPZ fields: spectral_s(n,), spectral_weights(n,),
resonance_masses(r,), resonance_widths(n,r), conformal_powers(c,) positive
integers, tail_description (scalar string), width_source (scalar string).
Optional conformal_s_inel_GeV2=1 and conformal_s1_GeV2=-1 use source defaults.
The finite grid/quadrature and all tail handling belong to the supplied
prescription; no tail is extrapolated or guessed. Endpoint truncation must
be explicitly stated in tail_description. Arrays for the existing profile
calculator can be included and are copied unchanged. This program supplies
only spectral templates, never masses, widths, observations or covariance.
"""

import argparse
import hashlib
import json
from pathlib import Path
import unittest

import numpy as np
import sympy as sp

if __package__:
    from .profile_vppp_dispersive_basis import derive_zero_width_identity
else:
    from profile_vppp_dispersive_basis import derive_zero_width_identity



def derive_formulas():
    """Derive absorptive resonance and upper-lip phase using native SymPy."""
    spectral_s, mass, width = sp.symbols('s M Gamma', positive=True)
    resonance = 1/(mass**2-spectral_s-sp.I*sp.sqrt(spectral_s)*width)
    absorptive = sp.factor(sp.im(sp.expand_complex(resonance)))
    threshold_distance, map_scale = sp.symbols('s_minus_s_inel s_inel_minus_s1', positive=True)
    upper_lip_root = -sp.I*sp.sqrt(threshold_distance)
    mapped = (sp.sqrt(map_scale)-upper_lip_root)/(sp.sqrt(map_scale)+upper_lip_root)
    real_part = sp.simplify(sp.re(sp.expand_complex(mapped)))
    imaginary_part = sp.simplify(sp.im(sp.expand_complex(mapped)))
    angle = sp.Symbol('theta', positive=True)
    substitutions = {threshold_distance: map_scale*sp.tan(angle)**2}
    # 0<theta<pi/2 fixes sqrt(tan(theta)^2)=tan(theta).
    real_angle = sp.trigsimp(real_part.subs(substitutions))
    imaginary_angle = sp.trigsimp(imaginary_part.subs(substitutions).refine(sp.Q.positive(sp.tan(angle))))
    phase_form = sp.cos(2*angle)+sp.I*sp.sin(2*angle)
    power = sp.Symbol('n', integer=True, positive=True)
    powered_phase = sp.expand_complex(sp.exp(2*sp.I*power*angle))
    imaginary_power = sp.simplify(sp.im(powered_phase))
    physical_phase = imaginary_power.subs(angle, sp.atan(sp.sqrt(threshold_distance/map_scale)))
    return [('Resonance input', resonance), ('Derived resonance imaginary part', absorptive),
        ('Upper-lip square root above threshold', upper_lip_root),
        ('Upper-lip conformal map', mapped), ('Derived real part of z', real_part),
        ('Derived imaginary part of z', imaginary_part),
        ('Substitute s-s_inel=(s_inel-s1) tan(theta)^2, real part', real_angle),
        ('Same substitution, imaginary part (0<theta<pi/2)', imaginary_angle),
        ('Conformal map on unit circle', phase_form),
        ('Integer power from exponential representation', powered_phase),
        ('Derived Im z^n above threshold', physical_phase)]


def conformal_map_upper_lip(spectral_s, s_inel=1., s1=-1.):
    """Evaluate the exact upper-lip square-root branch without an epsilon."""
    spectral_s = np.asarray(spectral_s, float)
    if not np.isfinite(s_inel) or not np.isfinite(s1) or s_inel <= 0 or s1 >= s_inel:
        raise ValueError('Conformal inputs require finite s_inel>0 and s1<s_inel')
    root_at_s1 = np.sqrt(s_inel-s1)
    root = np.sqrt(np.maximum(s_inel-spectral_s, 0.)).astype(complex)
    root -= 1j*np.sqrt(np.maximum(spectral_s-s_inel, 0.))
    return (root_at_s1-root)/(root_at_s1+root)


def build_templates(inputs):
    """Validate supplied grids and evaluate the two source spectral formulas."""
    required = ['spectral_s', 'spectral_weights', 'resonance_masses',
                'resonance_widths', 'conformal_powers', 'tail_description', 'width_source']
    missing = [name for name in required if name not in inputs]
    if missing:
        raise ValueError('Missing explicit spectral inputs: '+', '.join(missing))
    spectral_s = np.asarray(inputs['spectral_s'], float)
    weights = np.asarray(inputs['spectral_weights'], float)
    masses = np.asarray(inputs['resonance_masses'], float)
    widths = np.asarray(inputs['resonance_widths'], float)
    powers = np.asarray(inputs['conformal_powers'])
    if spectral_s.ndim != 1 or spectral_s.size < 2 or np.any(spectral_s <= 0) or np.any(np.diff(spectral_s) <= 0):
        raise ValueError('spectral_s must be strictly increasing and positive')
    if weights.shape != spectral_s.shape or np.any(weights < 0):
        raise ValueError('Explicit nonnegative quadrature weights must match spectral_s')
    if masses.ndim != 1 or np.any(masses <= 0) or widths.shape != (spectral_s.size, masses.size) or np.any(widths < 0):
        raise ValueError('Mass/width arrays must have dimensions (r,) and (n,r), with masses positive and widths nonnegative')
    if any(not np.isfinite(array).all() for array in [spectral_s, weights, masses, widths]):
        raise ValueError('Spectral inputs must be finite')
    if np.any(np.max(widths, axis=0) == 0):
        raise ValueError('A zero-width resonance requires a delta distribution, not a sampled finite-width grid')
    if powers.ndim != 1 or not np.issubdtype(powers.dtype, np.number) or np.any(powers < 1) or np.any(powers != np.floor(powers)):
        raise ValueError('conformal_powers must contain positive integers')
    if len(np.unique(powers)) != len(powers):
        raise ValueError('conformal_powers must be distinct')
    for name in ['tail_description', 'width_source']:
        value = np.asarray(inputs[name])
        if value.ndim != 0 or value.dtype.kind not in 'US' or not str(value.item()).strip():
            raise ValueError(name+' must be an explicit nonempty scalar string')
    s_inel = float(np.asarray(inputs.get('conformal_s_inel_GeV2', 1.)))
    s1 = float(np.asarray(inputs.get('conformal_s1_GeV2', -1.)))
    absorptive = np.sqrt(spectral_s)[:, None]*widths
    denominator = (masses[None, :]**2-spectral_s[:, None])**2+absorptive**2
    if np.any(denominator == 0):
        raise ValueError('Zero-width pole on spectral grid cannot be represented by ordinary samples')
    resonance_spectral = absorptive/denominator
    mapped = conformal_map_upper_lip(spectral_s, s_inel, s1)
    mapped_zero = conformal_map_upper_lip(0., s_inel, s1)
    conformal_spectral = (mapped[:, None]**powers[None, :]-mapped_zero**powers[None, :]).imag
    # Preserve caller-supplied profile fields without generating missing physical inputs.
    output = {name: np.asarray(inputs[name]) for name in inputs}
    output.update(resonance_spectral=resonance_spectral, conformal_spectral=conformal_spectral,
        conformal_s_inel_GeV2=np.array(s_inel), conformal_s1_GeV2=np.array(s1))
    return output


class SpectralTemplateTests(unittest.TestCase):
    """Algebraic tests; numerical fixtures are not physical resonance inputs."""

    def test_symbolic_spectral_derivations(self):
        formulas = dict(derive_formulas())
        spectral_s, mass, width = sp.symbols('s M Gamma', positive=True)
        self.assertEqual(sp.simplify(formulas['Derived resonance imaginary part']
            - sp.sqrt(spectral_s)*width/((mass**2-spectral_s)**2+spectral_s*width**2)), 0)
        angle = sp.Symbol('theta', positive=True)
        self.assertEqual(sp.trigsimp(formulas['Substitute s-s_inel=(s_inel-s1) tan(theta)^2, real part']-sp.cos(2*angle)), 0)
        self.assertEqual(sp.trigsimp(formulas['Same substitution, imaginary part (0<theta<pi/2)']-sp.sin(2*angle)), 0)

    def test_upper_lip_against_direct_complex_evaluation(self):
        invariant = np.array([.25, .75, 1.5, 3.])
        complex_root = np.sqrt(1.-(invariant+1e-12j))
        direct = (np.sqrt(2)-complex_root)/(np.sqrt(2)+complex_root)
        exact = conformal_map_upper_lip(invariant)
        np.testing.assert_allclose(exact, direct, rtol=1e-10, atol=1e-10)
        self.assertTrue(np.all(exact[invariant>1].imag > 0))
        self.assertTrue(np.all(exact[invariant<1].imag == 0))
        self.assertEqual(conformal_map_upper_lip(1.), 1.)

    def test_resonance_formula_and_conformal_threshold(self):
        inputs = dict(spectral_s=np.array([.5, 1., 2.]), spectral_weights=np.ones(3),
            resonance_masses=np.array([1.1]), resonance_widths=np.full((3, 1), .1),
            conformal_powers=np.array([1, 2, 3]), tail_description=np.array('algebraic fixture only'),
            width_source=np.array('algebraic fixture only'))
        result = build_templates(inputs)
        direct = (1/(1.1**2-inputs['spectral_s']-1j*np.sqrt(inputs['spectral_s'])*.1)).imag
        np.testing.assert_allclose(result['resonance_spectral'][:, 0], direct)
        np.testing.assert_array_equal(result['conformal_spectral'][:2], np.zeros((2, 3)))

    def test_exact_zero_width_limit_matches_existing_matching(self):
        _, difference = derive_zero_width_identity()
        self.assertEqual(difference, 0)

    def test_no_implicit_physical_inputs(self):
        with self.assertRaisesRegex(ValueError, 'Missing explicit spectral inputs'):
            build_templates({})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path, nargs='?')
    parser.add_argument('output', type=Path, nargs='?')
    parser.add_argument('--derive-formulas', '--derive-branches', action='store_true', help='derive source spectral formulas and print SymPy LaTeX')
    arguments = parser.parse_args()
    if arguments.derive_formulas:
        print('GeV units; upper lip above s_inel; real positive width; integer n>0.\n'
              'Below s_inel the conformal map and z(0) are real, so Im[z^n-z(0)^n]=0.\n'
              'Above threshold theta=atan(sqrt((s-s_inel)/(s_inel-s1))) lies in (0,pi/2).')
        for label, expression in derive_formulas():
            print(label+':\n'+sp.latex(expression)+'\n')
        return
    if arguments.input is None or arguments.output is None:
        parser.error('explicit input and output NPZ paths are required unless --derive-formulas is used')
    if arguments.output.suffix != '.npz':
        parser.error('output filename must end in .npz')
    if arguments.input.resolve() == arguments.output.resolve():
        parser.error('use separate input and output files to preserve source arrays')
    with np.load(arguments.input, allow_pickle=False) as inputs:
        output = build_templates(inputs)
    np.savez(arguments.output, **output)
    manifest = {'bibtex_keys': ['Hoferichter:2019mqg', 'Hoferichter:2025lcz'],
        'input_sha256': hashlib.sha256(arguments.input.read_bytes()).hexdigest(),
        'output_sha256': hashlib.sha256(arguments.output.read_bytes()).hexdigest(),
        'tail_description': str(output['tail_description'].item()),
        'width_source': str(output['width_source'].item()),
        'conformal_s_inel_GeV2': float(output['conformal_s_inel_GeV2']),
        'conformal_s1_GeV2': float(output['conformal_s1_GeV2']),
        'scope': 'Spectral templates from explicit caller inputs; no physical prediction or fit inferred.'}
    arguments.output.with_suffix('.provenance.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(arguments.output)
    print(arguments.output.with_suffix('.provenance.json'))


if __name__ == '__main__':
    main()
