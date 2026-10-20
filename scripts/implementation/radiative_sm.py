"""Action-derived SM V-P-photon amplitude, including c4 and all neutral exchange.

Tree level, single-trace RChT basis and common FV; explicit flavor matrices
are coefficients of Phi_R. Alda:2026nrd conversion is Phi_R=F Phi_A/(sqrt(2) fpi).
No physical eta/eta' default or fitted coupling is supplied here.
"""
from dataclasses import dataclass
import math
import numpy as np
from .odd_parity_rcht import OddParityCouplings, ideal_vector_generators, _source_matrix
from .radiative_odd_parity import _finite_real, _positive_real, vector_to_pseudoscalar_gamma_width

Q = np.diag([2/3,-1/3,-1/3])
Q.setflags(write=False)

@dataclass(frozen=True)
class RadiativeSMInputs:
    e: float
    F: float
    FV: float
    MV: float
    Mvector: float
    mp: float
    TV: np.ndarray
    TP: np.ndarray
    chi: np.ndarray
    internal_masses: tuple[float,float,float] | None = None

    def __post_init__(self):
        for name in ('e','F','FV','MV','Mvector'):
            object.__setattr__(self,name,_positive_real(getattr(self,name),name))
        object.__setattr__(self,'mp',_finite_real(self.mp,'mp'))
        if self.mp<0:
            raise ValueError('mp must be nonnegative')
        for name in ('TV','TP','chi'):
            object.__setattr__(self,name,_source_matrix(getattr(self,name),name))
        if not np.array_equal(self.chi,np.diag(np.diag(self.chi))) or np.any(np.imag(self.chi)):
            raise ValueError('chi must be a real diagonal mass-squared spurion')
        masses=(self.MV,)*3 if self.internal_masses is None else tuple(self.internal_masses)
        if len(masses)!=3:
            raise ValueError('internal_masses requires rho0, omega, phi masses')
        object.__setattr__(self,'internal_masses',tuple(_positive_real(v,'internal mass') for v in masses))

def _ac(a,b): return a@b+b@a
def _comm(a,b): return a@b-b@a

def radiative_sm_flavor(inputs: RadiativeSMInputs):
    """Return direct FQ,FH,FK and (overlap,FQ_alpha,FH_alpha,M_alpha) entries."""
    if not isinstance(inputs,RadiativeSMInputs):
        raise TypeError('inputs must be RadiativeSMInputs')
    x=inputs
    mass_wave=_ac(x.TP,x.chi)
    aq=_ac(x.TV,Q)
    FQ=np.trace(aq@x.TP)
    FH=np.trace(aq@mass_wave)
    FK=np.trace(x.TV@_comm(_comm(x.TP,Q),x.chi))
    generators=ideal_vector_generators()
    exchange=[]
    for name,mass in zip(('rho0','omega','phi'),x.internal_masses):
        T=generators[name]
        a=_ac(x.TV,T)
        exchange.append((np.trace(T@Q),np.trace(a@x.TP),np.trace(a@mass_wave),mass))
    return FQ,FH,FK,tuple(exchange)

def radiative_sm_components(couplings: OddParityCouplings, inputs: RadiativeSMInputs):
    """Return signed direct and exchange G in GeV^-1, Eq27 phase convention."""
    if not isinstance(couplings,OddParityCouplings):
        raise TypeError('couplings must be OddParityCouplings')
    FQ,FH,FK,exchange=radiative_sm_flavor(inputs)
    x=inputs; c,d=couplings.vjp,couplings.vvp
    def C():
        return FQ*((c.c1+c.c2-c.c5)*x.mp**2+(-c.c1+c.c2+c.c5-2*c.c6)*x.Mvector**2)+4*c.c3*FH-4*c.c4*FK
    direct=2*math.sqrt(2)*x.e/(x.F*x.Mvector*x.MV)*C()
    vector_sum=sum(w*(fq*((d.d1-d.d3)*x.mp**2+d.d3*x.Mvector**2)+4*d.d2*fh)/mass**2
                   for w,fq,fh,mass in exchange)
    exchanged=-4*x.e*x.FV/(x.F*x.Mvector)*vector_sum
    for value in (direct,exchanged):
        if not math.isfinite(value.real) or not math.isfinite(value.imag):
            raise ValueError('radiative amplitude must be finite')
    return complex(direct),complex(exchanged)

def G_radiative_sm(couplings: OddParityCouplings, inputs: RadiativeSMInputs):
    return sum(radiative_sm_components(couplings,inputs))

def radiative_sm_width(couplings: OddParityCouplings, inputs: RadiativeSMInputs, orientation='V_to_P'):
    """Open V->Pgamma or P->Vgamma partial width in GeV, spin averages explicit."""
    x=inputs; G=G_radiative_sm(couplings,inputs)
    if orientation=='V_to_P':
        return vector_to_pseudoscalar_gamma_width(G,x.Mvector,x.mp)
    if orientation!='P_to_V':
        raise ValueError('orientation must be V_to_P or P_to_V')
    if x.mp<=x.Mvector:
        raise ValueError('P_to_V requires mp > Mvector')
    result=abs(G)**2*(x.mp**2-x.Mvector**2)**3/(32*math.pi*x.mp**3)
    return _finite_real(result,'P_to_V partial width')
