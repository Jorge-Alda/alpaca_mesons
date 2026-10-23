"""Conditional omega-pion Dalitz fit along delta c1 = delta c5.

Run .venv/bin/python scripts/fit_vector_dalitz.py. All masses are GeV.
This preserves the preliminary fit's even parameters and measured odd image.
The input bins contain extracted normalized |F|^2. Diagonal statistical plus
systematic errors and bin-center evaluation are explicit approximations.
Canonical odd coefficients have a common opposite sign to the inherited
reference tensor backend; that common sign cancels in normalized |F|^2.
"""
from pathlib import Path
import json
import math
import unittest
import numpy as np
from scipy.integrate import quad
from scipy.optimize import least_squares
from scipy.linalg import svdvals
from implementation.vector_sector_amplitudes import VectorAmplitudeInputs, radiative_components
from implementation.radiative_odd_parity import _couplings_from_mapping
import derive_odd_low_energy_jacobians as calibration

ROOT = Path(__file__).resolve().parents[1]


def preliminary_inputs():
    """Reconstruct the preliminary even sector and Alda:2026nrd tree pion."""
    parameters = json.loads((ROOT/'data/vector_sector/odd_preliminary_fit.json').read_text())['parameters']
    assumptions = json.loads((ROOT/'data/vector_sector/stage3_results.json').read_text())['assumptions']
    point = np.array([parameters[name] for name in ('m_rho_GeV','m_Kstar_GeV','Delta_V_GeV2','kappa_V','F_V_GeV','z_V_GeV_inverse')])
    sector, _ = calibration.action_sector(point, 'single_trace_negative', assumptions)
    from implementation.vector_sector import VectorSector
    sector = VectorSector(sector.masses,sector.generators,sector.residues,sector.hypothesis)
    tree = calibration.module('alda_neutral_inputs').alda_2026nrd_neutral_inputs(mass_scheme='tree',F=calibration.F,fpi=calibration.F)['pi0']
    inputs = VectorAmplitudeInputs(calibration.F, calibration.OPERATOR_MV, tree['mp'],tree['TP'],tree['chi'],sector)
    return parameters, inputs


def normalized_form_factor_squared(parameters, inputs, mass):
    """Dimensionless |G(mass²)/G(0)|² with the full retained vector response."""
    couplings = _couplings_from_mapping(parameters)
    generator = inputs.sector.generators[1]
    vector_mass = inputs.sector.masses[1]
    amplitude = sum(radiative_components(couplings, inputs,generator,vector_mass,1.,mass*mass))
    real_amplitude = sum(radiative_components(couplings,inputs,generator,vector_mass,1.))
    return float(abs(amplitude/real_amplitude)**2)


def qed_weight(mass, vector_mass, pion_mass, lepton_mass):
    """Pointlike dGamma/dmass weight, with constants cancelling in bin ratios."""
    if mass <= 2*lepton_mass or mass >= vector_mass-pion_mass:
        return 0.
    squared = mass*mass
    triangle = (vector_mass**2+pion_mass**2-squared)**2-4*vector_mass**2*pion_mass**2
    return math.sqrt(1-4*lepton_mass**2/squared)*(1+2*lepton_mass**2/squared)*triangle**1.5/mass


def predict_bin(parameters,inputs,bin_data,lepton_mass,mode):
    """Evaluate at center or average with a pointlike physical-phase-space rate."""
    low,high = bin_data['mass_low_gev'],bin_data['mass_high_gev']
    if mode == 'center':
        return normalized_form_factor_squared(parameters,inputs,(low+high)/2)
    physical = json.loads((ROOT/'data/vector_sector/alda_2026nrd_arxiv_v1.json').read_text())['physical_masses_gev']
    weight = lambda mass: qed_weight(mass,inputs.sector.masses[1],physical['pi0'],lepton_mass)
    denominator = quad(weight,low,high,epsabs=1e-12)[0]
    return quad(lambda mass: weight(mass)*normalized_form_factor_squared(parameters,inputs,mass),low,high,epsabs=1e-12)[0]/denominator


def fit(experiments,mode='center'):
    """Fit the common c1,c5 shift at fixed remaining preliminary coordinates."""
    parameters,inputs = preliminary_inputs()
    bins = [(experiment,bin_data) for experiment in experiments for bin_data in experiment['bins']]
    def predictions(shift):
        shifted = dict(parameters,c1=parameters['c1']+shift,c5=parameters['c5']+shift)
        return np.array([predict_bin(shifted,inputs,b,experiment['lepton_mass_gev'],mode) for experiment,b in bins])
    observed = np.array([b['value'] for _,b in bins])
    # A2 supplies only the quoted combined uncertainty in this schema's
    # stat_sigma column; it must not be interpreted as a statistical error.
    errors = np.array([b.get('quoted_sigma',math.hypot(b['stat_sigma'],b['syst_sigma'])) for _,b in bins])
    residual = lambda coordinates: (predictions(coordinates[0])-observed)/errors
    fits = [least_squares(residual,[start],xtol=1e-12,gtol=1e-12,ftol=1e-12) for start in (-.1,0.,.1)]
    best = min(fits,key=lambda result: np.dot(result.fun,result.fun))
    shift = float(best.x[0])
    jacobian = best.jac
    result = dict(mode=mode,common_c1_c5_shift=shift,chi2=float(best.fun@best.fun),nominal_degrees_of_freedom=len(bins)-1,
                  local_shift_sigma=float(1/np.linalg.norm(jacobian)),weighted_singular_values=svdvals(jacobian).tolist(),
                  predictions=predictions(shift).tolist(),residuals=best.fun.tolist(),success=bool(best.success))
    result['c1_minus_c5'] = parameters['c1']-parameters['c5']
    result['c1_plus_c5'] = parameters['c1']+parameters['c5']+2*shift
    # Four exact preliminary null coordinates, in original coefficient notation.
    directions = ({'c1':1.,'c5':1.},{'c2':1.,'c5':1.,'c6':1.},{'c7':1.},{'d4':1.})
    accepted = dict(parameters,c1=parameters['c1']+shift,c5=parameters['c5']+shift)
    derivative_rows = []
    for direction in directions:
        changed = []
        for sign in (1.,-1.):
            shifted = {name:value+sign*1e-5*direction.get(name,0.) for name,value in accepted.items()}
            changed.append(np.array([predict_bin(shifted,inputs,b,experiment['lepton_mass_gev'],mode) for experiment,b in bins]))
        derivative_rows.append((changed[0]-changed[1])/(2e-5*errors))
    null_jacobian = np.array(derivative_rows).T
    from summarize_odd_preliminary_fit import main as preliminary_summary, NAMES
    _, orthonormal_null = preliminary_summary(write_output=False)
    named_null = np.array([[direction.get(name,0.) for direction in directions] for name in NAMES])
    conversion = np.linalg.lstsq(named_null,orthonormal_null,rcond=None)[0]
    defect = float(np.linalg.norm(named_null@conversion-orthonormal_null))
    if defect > 1e-8:
        raise ValueError('named candidate directions fail to span the measured null space')
    null_jacobian = null_jacobian@conversion
    normalization_columns = np.array([[prediction/errors[index] if experiment['bibtex_key']==key else 0. for key in [item['bibtex_key'] for item in experiments]] for index,((experiment,_),prediction) in enumerate(zip(bins,predictions(shift)))])
    projected = null_jacobian-normalization_columns@np.linalg.lstsq(normalization_columns,null_jacobian,rcond=None)[0]
    result['null_coordinates'] = ['delta_c1_equals_delta_c5','delta_c2_equals_delta_c5_equals_delta_c6','c7','d4']
    result['named_null_span_defect'] = defect
    result['null_jacobian_basis'] = 'orthonormal SVD basis from summarize_odd_preliminary_fit.main(write_output=False)'
    result['orthonormal_null_basis'] = orthonormal_null.tolist()
    result['null_jacobian_singular_values'] = svdvals(null_jacobian).tolist()
    result['normalization_projected_null_jacobian_singular_values'] = svdvals(projected).tolist()
    return result


def main():
    data = json.loads((ROOT/'data/odd_parity/vector_dalitz_measurements.json').read_text())
    results = []
    for mode in ('center','qed_weighted_average'):
        for experiments in ([data['experiments'][0]],[data['experiments'][1]],data['experiments']):
            result = fit(experiments,mode)
            result['experiments'] = [experiment['bibtex_key'] for experiment in experiments]
            results.append(result)
    output = dict(assumptions=['fixed preliminary even sector and measured odd image','zero-width tree vector poles','diagonal sum of statistical and systematic variances','no experimental systematic covariance provided','only the common c1,c5 null coordinate is fitted','bin-center result compared with pointlike-rate-weighted bin average'],results=results)
    path = ROOT/'data/vector_sector/vector_dalitz_conditional_fit.json'
    path.write_text(json.dumps(output,indent=2)+'\n')
    print(path.relative_to(ROOT))


class DalitzTests(unittest.TestCase):
    def test_finite_virtuality_polynomials_against_tensor_contractions(self):
        import sympy as sp
        from derive_radiative_sm_amplitude import derive,c,d,fq,fh,fk,m,s1,s2
        direct,exchange = derive()
        expected_direct = fq*((c[0]+c[1]-c[4])*m**2+(-c[0]+c[1]+c[4]-2*c[5])*s1+(c[0]-c[1]+c[4])*s2)+4*c[2]*fh-4*c[3]*fk
        expected_exchange = fq*((d[0]-d[2])*m**2+d[2]*(s1+s2))+4*d[1]*fh
        self.assertEqual(sp.simplify(direct-expected_direct),0)
        self.assertEqual(sp.simplify(exchange-expected_exchange),0)

    def test_real_photon_limit_and_null_coordinates(self):
        parameters,inputs = preliminary_inputs()
        self.assertAlmostEqual(normalized_form_factor_squared(parameters,inputs,0.),1.)
        shifted = dict(parameters)
        for name in ('c2','c5','c6'):
            shifted[name] += .032
        self.assertAlmostEqual(normalized_form_factor_squared(shifted,inputs,.5),normalized_form_factor_squared(parameters,inputs,.5),places=10)
        changed = dict(parameters,c1=parameters['c1']+.01,c5=parameters['c5']+.01)
        self.assertGreater(abs(normalized_form_factor_squared(changed,inputs,.5)-normalized_form_factor_squared(parameters,inputs,.5)),.01)

    def test_pointlike_thresholds(self):
        self.assertEqual(qed_weight(.2,.783,.135,.106),0.)
        self.assertEqual(qed_weight(.7,.783,.135,.106),0.)
        self.assertGreater(qed_weight(.4,.783,.135,.106),0.)


if __name__ == '__main__':
    main()
