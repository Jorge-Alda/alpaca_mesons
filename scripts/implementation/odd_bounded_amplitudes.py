"""Normalized neutral/charged SM cubic gradients of the bounded odd action.

GeV, (+---), source epsilon_lower0123=+1; coefficients multiply the common
source tensor E. Raw parent coefficients are dimensionless and redundant.
Use source-EOM-compatible LO pseudoscalar masses for the 228-action quotient;
phenomenological physical-mass substitutions require a separate EOM audit.
The compact input artifact contains the exact prior SymPy tensor derivation.
"""
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
import hashlib
import json
import math
import numpy as np
from .radiative_sm import Q
from .vector_sector_amplitudes import VectorAmplitudeInputs

ROOT = Path(__file__).resolve().parents[2]


@lru_cache(None)
def parent_catalog():
    path = ROOT/'data/vector_sector/odd_cubic_parent_kernels.json'
    saved = json.loads(path.read_text())
    for name, digest in saved['sha256'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest() != digest:
            raise ValueError('stale bounded cubic parent artifact: '+name)
    records = []
    for item in saved['parents']:
        record = dict(item)
        record['polynomial'] = np.array([float(Fraction(v))
                                        for v in item['kinematic_coefficients']])
        record['words'] = [(float(Fraction(term['coefficient'])), term['words'])
                           for term in item['flavor_words']]
        records.append(record)
    return tuple(records)


def parent_ids():
    return tuple(item['id'] for item in parent_catalog())


def retained_columns(scheme, ceiling):
    if scheme not in ('all_derivatives_soft', 'hard_vector_derivatives_promoted'):
        raise ValueError('unknown parent-counting scheme')
    if ceiling not in (0, 1, 2):
        raise ValueError('parent ceiling must be 0, 1 or 2')
    return np.array([i for i, item in enumerate(parent_catalog())
        if (0 if scheme == 'hard_vector_derivatives_promoted'
               and item['mode'] in ('Ydd', 'fmdd') else item['level'])
            + item['traces']-1 <= ceiling], dtype=int)


def _flavor(item, matrices):
    value = 0j
    for coefficient, words in item['words']:
        term = complex(coefficient)
        for word in words:
            product = matrices[word[0]]
            for name in word[1:]:
                product = product @ matrices[name]
            term *= np.trace(product)
        value += term
    return value


def _kernel(item, x, A, B, s, t):
    if item['mode'] in ('fminus', 'fmdd'):
        raise ValueError('composite fminus needs the radiative SM expansion')
    U = -math.sqrt(2)*x.TP/x.F
    Y = math.sqrt(2)*(x.TP@x.chi+x.chi@x.TP)/x.F
    # The exact flavor words use U as a generic slot in the Y blocks.
    source = Y if item['mode'] in ('Y', 'Ydd') else U
    matrices = dict(A=A, B=B, U=source, Y=Y, X=2*x.chi)
    polynomial = item['polynomial'] @ np.array([1., s, t, (x.mp**2-s-t)/2])
    return item['phase']*item['bose_factor']*_flavor(item, matrices)*polynomial


def _real_row(values):
    values = np.asarray(values, dtype=complex)
    if not np.all(np.isfinite(values)):
        raise ValueError('nonfinite bounded amplitude')
    if np.max(np.abs(values.imag), initial=0) > 1e-12*max(1., np.max(np.abs(values.real), initial=0)):
        raise ValueError('declared real tree channels have complex gradients')
    result = values.real
    result.setflags(write=False)
    return result


def radiative_gradients(x, TV, Mvector, e):
    """Raw dG/dg_j, GeV^-1; real photon, external canonical 1/Mvector."""
    if not isinstance(x, VectorAmplitudeInputs):
        raise TypeError('expected VectorAmplitudeInputs')
    if not np.isfinite(Mvector) or Mvector <= 0 or not np.isfinite(e) or e <= 0:
        raise ValueError('Mvector and e must be positive finite')
    TV = np.asarray(TV, dtype=complex)
    if TV.shape != (3, 3):
        raise ValueError('TV must be 3 by 3')
    W0 = x.sector.response(0.)
    values = []
    for item in parent_catalog():
        normalization = x.MV**(-item['operator_MV_power'])
        if item['mode'] in ('fminus', 'fmdd'):
            matrices = dict(A=TV, B=x.TP@Q-Q@x.TP, X=2*x.chi)
            polynomial = item['polynomial'] @ np.array(
                [1., Mvector**2, 0., (x.mp**2-Mvector**2)/2])
            value = (item['phase']*math.sqrt(2)*e/(x.F*Mvector)
                     * normalization*_flavor(item, matrices)*polynomial)
        elif item['kind'] == 'VJP':
            value = -2*e*normalization*_kernel(item, x, TV, Q, Mvector**2, 0.)/Mvector
        else:
            value = math.sqrt(2)*e*normalization*_kernel(item, x, TV, W0, Mvector**2, 0.)/Mvector
        values.append(value)
    return _real_row(values)


def tff_gradients(x, s1=0., s2=0.):
    """Raw dF/dg_j, GeV^-1; two photon placements, without WZW offset."""
    if not isinstance(x, VectorAmplitudeInputs):
        raise TypeError('expected VectorAmplitudeInputs')
    if not np.allclose(x.TP@Q-Q@x.TP, 0., atol=1e-14, rtol=0):
        raise ValueError('neutral electromagnetic-compatible TP required')
    W1, W2 = x.sector.response(s1), x.sector.response(s2)
    values = []
    for item in parent_catalog():
        normalization = x.MV**(-item['operator_MV_power'])
        if item['mode'] in ('fminus', 'fmdd'):
            value = 0.
        elif item['kind'] == 'VJP':
            # Responses are inserted after selecting each local polynomial.
            value = -2*math.sqrt(2)*normalization*(_kernel(item, x, W1, Q, s1, s2)
                                                  + _kernel(item, x, W2, Q, s2, s1))
        else:
            value = 2*normalization*_kernel(item, x, W1, W2, s1, s2)
        values.append(value)
    return _real_row(values)


def wzw_offset(x):
    """Inherited leading nonet WZW normalization, GeV^-1, common E phase."""
    return float((math.sqrt(2)*x.NC*np.trace(2*Q@Q@x.TP)/(8*math.pi**2*x.F)).real)


@lru_cache(None)
def original_embedding():
    """241-by-11 original c1..c7,d1..d4 embedding, no OPE elimination."""
    saved = json.loads((ROOT/'data/vector_sector/odd_leading_action_count.json').read_text())
    ids = parent_ids()
    result = np.zeros((len(ids), 11))
    # Leading single-trace families are identified by their saved raw words.
    def single(prefix):
        candidates = [item for item in parent_catalog() if item['id'].startswith(prefix)]
        item = next(item for item in candidates if item['traces'] == 1)
        return item['id'].split(':')[-2]
    jf, vf = single('VJP:F'), single('VVP:S:F')
    for j, coordinates in zip((0, 1, 4, 5, 6), saved['original_VJP_derivatives_to_generic_basis']):
        for k, value in enumerate(coordinates):
            result[ids.index(f'VJP:{jf}:L{k}'), j] += 2*float(Fraction(value))
    result[ids.index(f'VJP:Y:{jf}'), 2] = 2.
    result[ids.index('VJP:fminus'), 3] = 2.
    for j, coordinates in zip((7, 9, 10), saved['original_VVP_derivatives_to_symmetric_basis']):
        for k, value in enumerate(coordinates):
            result[ids.index(f'VVP:S:{vf}:L{k}'), j] += float(Fraction(value))
    result[ids.index(f'VVP:Y:{vf}'), 8] = 2.
    result.setflags(write=False)
    return result
