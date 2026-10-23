"""Neutral spectral response inserted in the established single-trace kernels.

chi here is chi=2B0 mq, whereas vector source chi_plus=2 chi in vacuum.
MV remains the reference operator scale, independently of pole masses.
"""
from dataclasses import dataclass
import math
import numpy as np
from .odd_parity_rcht import OddParityCouplings, _source_matrix
from .radiative_odd_parity import _finite_real, _positive_real
from .radiative_odd_parity import _couplings_from_mapping
from .radiative_sm import Q, _ac, _comm
from .vector_sector import VectorSector

@dataclass(frozen=True)
class VectorAmplitudeInputs:
    F: float
    MV: float
    mp: float
    TP: np.ndarray
    chi: np.ndarray
    sector: VectorSector
    NC: int = 3

    def __post_init__(self):
        for n in ('F','MV'): object.__setattr__(self,n,_positive_real(getattr(self,n),n))
        object.__setattr__(self,'mp',_finite_real(self.mp,'mp'))
        if self.mp<0: raise ValueError('mp must be nonnegative')
        for n in ('TP','chi'): object.__setattr__(self,n,_source_matrix(getattr(self,n),n))
        if not np.allclose(self.chi,np.diag(np.diag(self.chi))) or np.any(self.chi.imag):
            raise ValueError('chi must be real diagonal')
        if not isinstance(self.sector,VectorSector): raise TypeError('sector must be VectorSector')
        if isinstance(self.NC,bool) or not isinstance(self.NC,(int,np.integer)) or self.NC<=0:
            raise ValueError('NC must be positive integer')


def flavor(x,A,B):
    return np.trace(_ac(A,B)@x.TP), np.trace(_ac(A,B)@_ac(x.TP,x.chi))


def _C(c,x,A,B,s,t):
    fq,fh=flavor(x,A,B)
    return fq*((c.c1+c.c2-c.c5)*x.mp**2+(-c.c1+c.c2+c.c5-2*c.c6)*s
               +(c.c1-c.c2+c.c5)*t)+4*c.c3*fh


def _H(d,x,A,B,s,t):
    fq,fh=flavor(x,A,B)
    return fq*((d.d1-d.d3)*x.mp**2+d.d3*(s+t))+4*d.d2*fh


def tff_components(couplings,x,s1=0.,s2=0.):
    """Signed WZW,VJP,VVP, GeV^-1; neutral TP required."""
    if not isinstance(couplings,OddParityCouplings) or not isinstance(x,VectorAmplitudeInputs):
        raise TypeError('expected OddParityCouplings and VectorAmplitudeInputs')
    if not np.allclose(_comm(x.TP,Q),0,atol=1e-14,rtol=0):
        raise ValueError('TFF requires [TP,Q]=0')
    s1,s2=_finite_real(s1,'s1'),_finite_real(s2,'s2')
    w1,w2=x.sector.response(s1),x.sector.response(s2)
    values=(math.sqrt(2)*x.NC*flavor(x,Q,Q)[0]/(8*math.pi**2*x.F),
            4/(x.F*x.MV)*(_C(couplings.vjp,x,w1,Q,s1,s2)+_C(couplings.vjp,x,w2,Q,s2,s1)),
            -4*math.sqrt(2)/x.F*_H(couplings.vvp,x,w1,w2,s1,s2))
    if not np.all(np.isfinite(values)): raise ValueError('TFF components must be finite')
    return values


def radiative_components(couplings,x,TV,Mvector,e,photon_mass_squared=0.):
    """Signed direct/exchange G in GeV^-1 at photon invariant mass squared.

    The external tensor reduction uses 1/Mvector. The virtual photon is
    contracted with a conserved electromagnetic current. Tree poles and
    source residues follow the supplied sector, without an inserted width.
    """
    if not isinstance(couplings,OddParityCouplings) or not isinstance(x,VectorAmplitudeInputs):
        raise TypeError('expected OddParityCouplings and VectorAmplitudeInputs')
    TV=_source_matrix(TV,'TV'); Mvector=_positive_real(Mvector,'Mvector');e=_positive_real(e,'e')
    fk=np.trace(TV@_comm(_comm(x.TP,Q),x.chi))
    direct=2*math.sqrt(2)*e/(x.F*Mvector*x.MV)*(_C(couplings.vjp,x,TV,Q,Mvector**2,photon_mass_squared)-4*couplings.vjp.c4*fk)
    exchange=-4*e/(x.F*Mvector)*_H(couplings.vvp,x,TV,x.sector.response(photon_mass_squared),Mvector**2,photon_mass_squared)
    if not np.all(np.isfinite((direct,exchange))): raise ValueError('radiative components must be finite')
    return complex(direct),complex(exchange)


def two_photon_width(couplings,x,alpha):
    """P->gamma gamma width in GeV, including identical-photon factor."""
    alpha=_positive_real(alpha,'alpha')
    form_factor=sum(tff_components(couplings,x,0,0))
    return _finite_real(math.pi*alpha**2*x.mp**3*abs(form_factor)**2/4,'two-photon width')


def q2_abs_tff_predictor(parameters,config):
    """Fit-ready Q^2 |F(-Q^2,0)| in GeV, with explicit flavor inputs."""
    if set(config)!={'inputs','q2'} or not isinstance(config['inputs'],VectorAmplitudeInputs):
        raise ValueError('TFF config requires inputs=VectorAmplitudeInputs and q2')
    if np.iscomplexobj(config['q2']):
        raise ValueError('q2 must be finite, real and nonnegative')
    q2=np.asarray(config['q2'],dtype=float)
    if q2.ndim!=1 or q2.size==0 or not np.all(np.isfinite(q2)) or np.any(q2<0):
        raise ValueError('q2 must be a finite nonempty 1D nonnegative array')
    couplings=_couplings_from_mapping(parameters)
    result=q2*np.abs(np.array([sum(tff_components(couplings,config['inputs'],-q,0)) for q in q2]))
    if not np.all(np.isfinite(result)): raise ValueError('prediction must be finite')
    result.setflags(write=False)
    return result


def radiative_width(couplings,x,TV,Mvector,e,orientation='V_to_P'):
    """Crossed radiative partial width in GeV with explicit spin averages."""
    Mvector=_positive_real(Mvector,'Mvector')
    G=sum(radiative_components(couplings,x,TV,Mvector,e))
    if orientation=='V_to_P':
        if Mvector<=x.mp: raise ValueError('V_to_P requires Mvector > mp')
        width=abs(G)**2*(Mvector**2-x.mp**2)**3/(96*math.pi*Mvector**3)
    elif orientation=='P_to_V':
        if x.mp<=Mvector: raise ValueError('P_to_V requires mp > Mvector')
        width=abs(G)**2*(x.mp**2-Mvector**2)**3/(32*math.pi*x.mp**3)
    else:
        raise ValueError('orientation must be V_to_P or P_to_V')
    return _finite_real(width,'radiative width')


def radiative_width_predictor(parameters,config):
    """One-entry immutable prediction for a declared crossed channel."""
    if set(config)!={'inputs','TV','Mvector','e','orientation'} or not isinstance(config['inputs'],VectorAmplitudeInputs):
        raise ValueError('radiative config requires inputs, TV, Mvector, e, orientation')
    result=np.array([radiative_width(_couplings_from_mapping(parameters),config['inputs'],
                    config['TV'],config['Mvector'],config['e'],config['orientation'])])
    result.setflags(write=False)
    return result
