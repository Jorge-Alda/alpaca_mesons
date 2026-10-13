"""Derived tree-level neutral-P TFF with one common vector mass.

Derivation: scripts/derive_radiative_tff_amplitude.py and indexed memory record.
The signed convention has a positive pi0 WZW contribution. Data measure |F|.
Inputs explicitly supply FQ=Tr({Q,Q}TP), FH=Tr({Q,Q}{TP,chi}) in GeV^2,
where chi=2B0 mq. Physical eta/eta' masses do not define FH.
"""
from dataclasses import dataclass
import math
from typing import Mapping, Any
import numpy as np
from .odd_parity_rcht import OddParityCouplings
from .radiative_odd_parity import _couplings_from_mapping, _finite_real, _positive_real

@dataclass(frozen=True)
class TFFInputs:
    F: float
    FV: float
    MV: float
    mp: float
    FQ: float
    FH: float
    NC: int = 3

    def __post_init__(self):
        for name in ('F','FV','MV'):
            object.__setattr__(self,name,_positive_real(getattr(self,name),name))
        for name in ('mp','FQ','FH'):
            object.__setattr__(self,name,_finite_real(getattr(self,name),name))
        if self.mp<0:
            raise ValueError('mp must be nonnegative')
        if isinstance(self.NC,bool) or not isinstance(self.NC,(int,np.integer)) or self.NC<=0:
            raise ValueError('NC must be a positive integer')

def tff_components(couplings: OddParityCouplings, inputs: TFFInputs, s1=0., s2=0.):
    """Return signed (WZW,VJP,VVP), each in GeV^-1; virtualities in GeV^2.

    Tree amplitudes at a vector pole are undefined. Widths/loops are absent.
    """
    if not isinstance(couplings,OddParityCouplings) or not isinstance(inputs,TFFInputs):
        raise TypeError('couplings and inputs must be OddParityCouplings and TFFInputs')
    s1,s2=_finite_real(s1,'s1'),_finite_real(s2,'s2')
    x=inputs; c,d=couplings.vjp,couplings.vvp
    D1,D2=x.MV**2-s1,x.MV**2-s2
    if D1==0 or D2==0:
        raise ValueError('tree TFF is undefined at a vector pole')
    def C(s,t):
        return x.FQ*((c.c1+c.c2-c.c5)*x.mp**2+(-c.c1+c.c2+c.c5-2*c.c6)*s
                     +(c.c1-c.c2+c.c5)*t)+4*c.c3*x.FH
    H=x.FQ*((d.d1-d.d3)*x.mp**2+d.d3*(s1+s2))+4*d.d2*x.FH
    values=(math.sqrt(2)*x.NC*x.FQ/(8*math.pi**2*x.F),
            4*x.FV/(x.F*x.MV)*(C(s1,s2)/D1+C(s2,s1)/D2),
            -4*math.sqrt(2)*x.FV**2/x.F*H/(D1*D2))
    return tuple(_finite_real(v,'TFF component') for v in values)

def pseudoscalar_tff(couplings: OddParityCouplings, inputs: TFFInputs, q2):
    """Signed spacelike F(-Q^2,0), GeV^-1; Q^2 >=0, GeV^2."""
    values=np.asarray(q2,dtype=float)
    if np.iscomplexobj(q2) or not np.all(np.isfinite(values)) or np.any(values<0):
        raise ValueError('q2 must be finite, real and nonnegative')
    result=np.array([sum(tff_components(couplings,inputs,-q,0.)) for q in values.flat]).reshape(values.shape)
    return float(result) if result.ndim==0 else result

def q2_abs_tff_predictor(parameters: Mapping[str,float], config: Mapping[str,Any]):
    """Fit predictor for the published Q^2|F| observable in GeV."""
    if set(config)!={'inputs','q2'} or not isinstance(config['inputs'],TFFInputs):
        raise ValueError('TFF config must contain inputs=TFFInputs and q2')
    q2=np.asarray(config['q2'],dtype=float)
    if q2.ndim!=1 or q2.size==0:
        raise ValueError('q2 must be a nonempty one-dimensional array')
    prediction=q2*np.abs(pseudoscalar_tff(_couplings_from_mapping(parameters),config['inputs'],q2))
    prediction.setflags(write=False)
    return prediction
