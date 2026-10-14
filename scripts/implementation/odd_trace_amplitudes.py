"""Normalized scalar SM kernels for the leading single/double-trace sector.

No chi_plus-dressed derivative extension or ALP/source-functional predictor
is supplied here. Two axial-source-only antisymmetric VVP directions are
not represented by this scalar projection; no EOM/OPE relation is imposed.
"""
from dataclasses import dataclass
import math
import numpy as np
from .odd_parity_rcht import OddParityCouplings, VJPCouplings, VVPCouplings, _source_matrix
from .radiative_odd_parity import _finite_real, _positive_real
from .radiative_sm import Q, _ac, _comm
from .vector_sector_amplitudes import VectorAmplitudeInputs

ZC=VJPCouplings(*([0.]*7)); ZD=VVPCouplings(*([0.]*4))

@dataclass(frozen=True)
class OddTraceCouplings:
    """Four JP trace families; three symmetric VV families; one scalar VV odd pair.

    q_JP=2*(F0,F1,F2,F3); q_VV=(2F0,F1+F2,2F3).
    VJP c4 belongs only to the original commutator. c7/d4 are retained
    as scalar-blind slots to preserve the original baseline correspondence.
    """
    vjp: tuple
    vvp: tuple
    d_minus: float=0.

    def __post_init__(self):
        vjp,vvp=tuple(self.vjp),tuple(self.vvp)
        if len(vjp)!=4 or not all(isinstance(c,VJPCouplings) for c in vjp):
            raise ValueError('vjp requires four VJPCouplings trace families')
        if len(vvp)!=3 or not all(isinstance(d,VVPCouplings) for d in vvp):
            raise ValueError('vvp requires three VVPCouplings symmetric trace families')
        if any(c.c4!=0 for c in vjp[1:]):
            raise ValueError('no independent double-trace c4: its C-odd projection vanishes')
        object.__setattr__(self,'vjp',vjp);object.__setattr__(self,'vvp',vvp)
        object.__setattr__(self,'d_minus',_finite_real(self.d_minus,'d_minus'))

    @classmethod
    def from_baseline(cls,couplings):
        if not isinstance(couplings,OddParityCouplings):raise TypeError('expected OddParityCouplings')
        return cls((couplings.vjp,ZC,ZC,ZC),(couplings.vvp,ZD,ZD))

def trace_flavor(A,B,Z):
    """F0,F1,F2,F3; generic matrices, including charged external conjugations."""
    A,B,Z=(_source_matrix(v,n) for v,n in ((A,'A'),(B,'B'),(Z,'Z')))
    return np.array([np.trace(_ac(A,B)@Z)/2,
        np.trace(A)*np.trace(B@Z),np.trace(B)*np.trace(A@Z),np.trace(A@B)*np.trace(Z)])

def _inputs(c,x):
    if not isinstance(c,OddTraceCouplings) or not isinstance(x,VectorAmplitudeInputs):
        raise TypeError('expected OddTraceCouplings and VectorAmplitudeInputs')

def C(c,x,A,B,s,t):
    q=2*trace_flavor(A,B,x.TP); h=2*trace_flavor(A,B,_ac(x.TP,x.chi))
    return sum(qi*((v.c1+v.c2-v.c5)*x.mp**2+(-v.c1+v.c2+v.c5-2*v.c6)*s
                   +(v.c1-v.c2+v.c5)*t)+4*v.c3*hi
               for v,qi,hi in zip(c.vjp,q,h))

def H(c,x,A,B,s,t):
    fq=trace_flavor(A,B,x.TP);fh=trace_flavor(A,B,_ac(x.TP,x.chi))
    q=np.array([2*fq[0],fq[1]+fq[2],2*fq[3]])
    h=np.array([2*fh[0],fh[1]+fh[2],2*fh[3]])
    symmetric=sum(qi*((d.d1-d.d3)*x.mp**2+d.d3*(s+t))+4*d.d2*hi
                  for d,qi,hi in zip(c.vvp,q,h))
    return symmetric-c.d_minus*(fq[1]-fq[2])*(s-t)

def tff_components(c,x,s1=0.,s2=0.):
    _inputs(c,x);s1=_finite_real(s1,'s1');s2=_finite_real(s2,'s2')
    if not np.allclose(_comm(x.TP,Q),0,atol=1e-14,rtol=0):
        raise ValueError('neutral TFF requires [TP,Q]=0')
    w1,w2=x.sector.response(s1),x.sector.response(s2)
    result=(math.sqrt(2)*x.NC*np.trace(_ac(Q,Q)@x.TP)/(8*math.pi**2*x.F),
        4/(x.F*x.MV)*(C(c,x,w1,Q,s1,s2)+C(c,x,w2,Q,s2,s1)),
        -4*math.sqrt(2)/x.F*H(c,x,w1,w2,s1,s2))
    if not np.all(np.isfinite(result)):raise ValueError('TFF prediction must be finite')
    return result

def radiative_components(c,x,TV,Mvector,e):
    _inputs(c,x);TV=_source_matrix(TV,'TV')
    Mvector=_positive_real(Mvector,'Mvector');e=_positive_real(e,'e')
    fk=np.trace(TV@_comm(_comm(x.TP,Q),x.chi))
    result=(2*math.sqrt(2)*e/(x.F*Mvector*x.MV)*(C(c,x,TV,Q,Mvector**2,0)-4*c.vjp[0].c4*fk),
            -4*e/(x.F*Mvector)*H(c,x,TV,x.sector.response(0),Mvector**2,0))
    if not np.all(np.isfinite(result)):raise ValueError('radiative prediction must be finite')
    return result

def two_photon_width(c,x,alpha):
    alpha=_positive_real(alpha,'alpha')
    return _finite_real(math.pi*alpha**2*x.mp**3*abs(sum(tff_components(c,x)))**2/4,'two photon width')

def radiative_width(c,x,TV,Mvector,e,orientation='V_to_P'):
    Mvector=_positive_real(Mvector,'Mvector');amplitude=abs(sum(radiative_components(c,x,TV,Mvector,e)))**2
    if orientation=='V_to_P' and Mvector>x.mp:
        value=amplitude*(Mvector**2-x.mp**2)**3/(96*math.pi*Mvector**3)
    elif orientation=='P_to_V' and x.mp>Mvector:
        value=amplitude*(x.mp**2-Mvector**2)**3/(32*math.pi*x.mp**3)
    else:raise ValueError('invalid orientation or closed decay channel')
    return _finite_real(value,'radiative width')
