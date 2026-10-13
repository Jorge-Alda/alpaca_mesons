"""Declared tree tensor-vector sector; masses/residues GeV, chi_plus GeV^2.

Canonical neutral generators obey Tr(T_i T_j)=delta_ij. No loop widths,
kinetic breaking, or automatic uncertainty averaging are introduced.
"""
from dataclasses import dataclass
from pathlib import Path
import json
import hashlib
import numpy as np
from .odd_parity_rcht import ideal_vector_generators, _source_matrix
from .radiative_odd_parity import _finite_real, _positive_real
from .radiative_sm import Q

NAMES = ('rho0', 'omega', 'phi')


def physical_generators(theta):
    theta = _finite_real(theta, 'theta')
    g = ideal_vector_generators()
    return (g['rho0'], np.cos(theta)*g['omega']-np.sin(theta)*g['phi'],
            np.sin(theta)*g['omega']+np.cos(theta)*g['phi'])


def source_residues(generators, chi_plus, FV, z=0., y=0.):
    """Exact declared single- and double-trace source action, signed GeV."""
    x = _source_matrix(chi_plus, 'chi_plus')
    if not np.allclose(x, np.diag(np.diag(x))) or np.any(x.imag):
        raise ValueError('chi_plus must be real diagonal')
    FV,z,y = (_finite_real(v,n) for v,n in ((FV,'FV'),(z,'z'),(y,'y')))
    source = FV*Q + z*(x@Q+Q@x) + y*np.trace(x@Q)*np.eye(3)
    traces=tuple(np.trace(_source_matrix(t,'generator')@source) for t in generators)
    if any(abs(t.imag)>1e-14 for t in traces):
        raise ValueError('source residues must be real; check generator phases')
    return tuple(float(t.real) for t in traces)


def neutral_mass_matrix(MV, em, delta, kappa, xl, xs):
    """Mass matrix from the declared stage-2 invariant action, in GeV^2."""
    MV = _positive_real(MV,'MV')
    em,delta,kappa,xl,xs = (_finite_real(v,n) for v,n in
                          ((em,'em'),(delta,'delta'),(kappa,'kappa'),(xl,'xl'),(xs,'xs')))
    return np.array([[MV**2-4*em*xl+2*delta/3+2*kappa*xl,
                      np.sqrt(2)*delta/3+kappa*(xl+xs)/np.sqrt(2)],
                     [np.sqrt(2)*delta/3+kappa*(xl+xs)/np.sqrt(2),
                      MV**2-4*em*xs+delta/3+kappa*xs]])


@dataclass(frozen=True)
class VectorSector:
    masses: tuple
    generators: tuple
    residues: tuple
    hypothesis: str = 'explicit'

    def __post_init__(self):
        if len(self.masses)!=3 or len(self.generators)!=3 or len(self.residues)!=3:
            raise ValueError('neutral sector requires three states')
        masses=tuple(_positive_real(m,'mass') for m in self.masses)
        generators=tuple(_source_matrix(t,'generator') for t in self.generators)
        for t in generators:
            if not np.allclose(t,t.conj().T) or not np.allclose(t,np.diag(np.diag(t))):
                raise ValueError('neutral generators must be Hermitian diagonal')
        gram=np.array([[np.trace(a@b) for b in generators] for a in generators])
        if not np.allclose(gram,np.eye(3),atol=1e-12,rtol=1e-12):
            raise ValueError('generators must be orthonormal')
        object.__setattr__(self,'masses',masses)
        object.__setattr__(self,'generators',generators)
        object.__setattr__(self,'residues',tuple(_finite_real(r,'residue') for r in self.residues))

    def response(self,s=0.):
        """W(s)=sum R_i T_i/(M_i^2-s), GeV^-1; pole undefined."""
        s=_finite_real(s,'s')
        denominators=np.array(self.masses)**2-s
        if np.any(denominators==0):
            raise ValueError('tree response is undefined at a vector pole')
        return sum((r*t/d for r,t,d in zip(self.residues,self.generators,denominators)),
                   np.zeros((3,3),complex))

    @classmethod
    def universal(cls,MV,FV,theta=0.):
        g=physical_generators(theta)
        return cls((MV,)*3,g,source_residues(g,np.zeros((3,3)),FV),'universal')

    @classmethod
    def from_stage3(cls, *, rho_input, kstar_input, phi_input,
                    source_hypothesis, branch='near_ideal', results_path=None, inputs_path=None):
        """Explicit separate scenario; nearest ideal is a declared tie-break policy.

        branch may be 'near_ideal', 'other', or the retained integer root index.
        Width signs are those of the explicitly selected source hypothesis.
        """
        root=Path(__file__).resolve().parents[2]
        data=json.loads(Path(results_path or root/'data/vector_sector/stage3_results.json').read_text())
        input_bytes=Path(inputs_path or root/'data/vector_sector/pdg2024.json').read_bytes()
        if hashlib.sha256(input_bytes).hexdigest()!=data['input_sha256']:
            raise ValueError('stage-3 inputs checksum mismatch')
        inputs=json.loads(input_bytes)
        cases=[s for s in data['scenarios'] if s['rho_input']==rho_input and s['kstar_input']==kstar_input]
        if len(cases)!=1: raise ValueError('unknown mass scenario')
        scenario=cases[0]; branches=scenario['extended_mass_branches']
        order=sorted(range(len(branches)),key=lambda i:abs(branches[i]['theta_deg']))
        index=order[0] if branch=='near_ideal' else order[-1] if branch=='other' else branch
        if isinstance(index,bool) or not isinstance(index,int) or index not in range(len(branches)):
            raise ValueError('invalid mass branch')
        widths=[w for w in scenario['width_tests'] if w['phi_input']==phi_input]
        if len(widths)!=1: raise ValueError('unknown phi width scenario')
        case=next(c for c in widths[0]['fixed_angle_cases'] if c['case']==f'extended_mass_branch_{index}')
        if not np.isclose(case['theta_deg'],branches[index]['theta_deg'],rtol=0,atol=1e-12):
            raise ValueError('stage-3 source and mass angles differ')
        mapping={'universal':'universal','single_trace_positive':'single_trace_spurion',
                 'single_trace_negative':'single_trace_spurion_sign_alternative',
                 'double_trace':'plus_double_trace_spurion'}
        if source_hypothesis not in mapping: raise ValueError('explicit source hypothesis required')
        p=case[mapping[source_hypothesis]]
        g=physical_generators(np.radians(branches[index]['theta_deg']))
        a=data['assumptions']; chi_plus=np.diag([a['xl_gev2']]*2+[a['xs_gev2']])
        residues=source_residues(g,chi_plus,p['FV_gev'],p.get('z_gev_inverse',0.),p.get('y_gev_inverse',0.))
        masses=tuple(inputs['masses'][key]['value_mev']/1000 for key in (rho_input,'omega','phi'))
        return cls(masses,g,residues,f'{rho_input}/{kstar_input}/{phi_input}/branch_{index}/{source_hypothesis}')
