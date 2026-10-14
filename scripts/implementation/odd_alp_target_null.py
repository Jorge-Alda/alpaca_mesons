"""Conditional neutral ALP resonance gradients; GeV and common source E.

The photon coefficient factors out e**2; Vgamma retains e. These are local
bounded-action projections, not widths or independent anomaly matching.
"""
import math
import json
from pathlib import Path
from functools import lru_cache
import numpy as np
import sympy as sp
from . import odd_bounded_amplitudes as API

@lru_cache(None)
def action_relations():
    root=Path(__file__).resolve().parents[2]
    leading=json.loads((root/'data/vector_sector/odd_leading_action_count.json').read_text())
    higher=json.loads((root/'data/vector_sector/odd_three_flavor_action_quotient.json').read_text())
    relations=leading['relations']+higher['all_higher_relations']
    ids=API.parent_ids()
    return np.asarray([[float(sp.Rational(relation.get(identifier,'0'))) for identifier in ids] for relation in relations])

def relation_certificate(rows):
    matrix=np.asarray(list(rows.values()))
    relations=action_relations()
    residual=matrix@relations.T
    scales=np.maximum(1.,np.linalg.norm(matrix,axis=1)[:,None]*np.linalg.norm(relations,axis=1)[None,:])
    defect=float(np.max(np.abs(residual)/scales))
    return dict(relation_count=len(relations),maximum_normalized_relation_defect=defect,valid=bool(defect<1e-10))


def sources(chi, mass_squared, cAP, cGA, m0_squared, fa=1.):
    """LO diagonal physical sources, excluding neutral poles; fa in GeV.

    cAP and cGA are independent UV invariant inputs. Setting fa=1 GeV is
    an overall normalization for ranks, never a fitted decay constant.
    """
    chi = np.asarray(chi, float)
    if chi.shape != (3, 3) or not np.allclose(chi, np.diag(np.diag(chi))):
        raise ValueError('diagonal 3x3 quark spurion required')
    if fa <= 0:
        raise ValueError('positive fa required')
    g = m0_squared/3
    K = np.diag(np.diag(chi)-mass_squared)+g*np.ones((3, 3))
    if np.linalg.cond(K) > 1e10:
        raise ValueError('neutral pole or ill-conditioned resolvent')
    S = np.linalg.solve(K, -np.diag(chi)*np.asarray(cAP)-g*cGA)
    U = np.diag(-2*S/fa)
    Y = np.diag(4*(np.asarray(cAP)+S)*np.diag(chi)/fa)
    Z = -2*g*(cGA+sum(S))/fa
    if not np.allclose(Y+2*mass_squared*U-2*Z*np.eye(3), 0., atol=1e-12):
        raise ValueError('source EOM defect')
    return U, Y, Z


def kernel(item, chi, U, Y, A, B, w, s, t, MV=.770):
    if item['mode'] in ('fminus', 'fmdd'):
        return 0.
    slot = Y if item['mode'] in ('Y', 'Ydd') else U
    value = API._flavor(item, dict(A=A, B=B, U=slot, Y=Y, X=2*chi))
    polynomial = item['polynomial'] @ np.array([1., s, t, (w-s-t)/2])
    return (item['phase']*item['bose_factor']*value*polynomial
            * MV**(-item['operator_MV_power']))


def target_rows(chi, U, Y, w, sector, e, MV=.770):
    """Each labeled on-shell neutral target, 241 parent coefficients.

    Includes signed photon response and both photon assignments. VV is
    the VVP resonance sector only. Odd one-insertion vacuum 3P is zero.
    Daughter kinematics are on-shell algebraic projections; no threshold
    or unstable-vector spectral conclusion is asserted.
    """
    W = sector.response(0.)
    names = ('rho', 'omega', 'phi')
    rows = {}
    def row(evaluator):
        return API._real_row([evaluator(item) for item in API.parent_catalog()])
    for i, (A, m1) in enumerate(zip(sector.generators, sector.masses)):
        for j in range(i, 3):
            B, m2 = sector.generators[j], sector.masses[j]
            rows['VV:'+names[i]+':'+names[j]] = row(lambda p:
                kernel(p, chi, U, Y, A, B, w, m1*m1, m2*m2, MV)/(m1*m2)
                if p['kind'] == 'VVP' else 0.)
        rows['Vgamma:'+names[i]] = row(lambda p:
            (-2*e*kernel(p, chi, U, Y, A, API.Q, w, m1*m1, 0., MV)
             if p['kind'] == 'VJP' else math.sqrt(2)*e*kernel(
                 p, chi, U, Y, A, W, w, m1*m1, 0., MV))/m1)
    rows['gammagamma'] = row(lambda p:
        -4*math.sqrt(2)*kernel(p, chi, U, Y, W, API.Q, w, 0., 0., MV)
        if p['kind'] == 'VJP' else 2*kernel(p, chi, U, Y, W, W, w, 0., 0., MV))
    rows['3P'] = np.zeros(len(API.parent_ids()))
    return rows
