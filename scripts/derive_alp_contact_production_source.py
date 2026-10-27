"""Derive unprojected original-operator ALP V_i gamma sources.

Natural units; metric +---. All incoming k+q+p=0. Work in k rest frame
with q=(q0,0,0,Q), transverse photon epsilon=(0,0,1,0).
Tensor-normalized sources retain the Fourier phases of the existing
ordered-numerator derivation. Epsilon orientation and source prefactors
are fixed by reproducing its wedge projection, not by importing old LECs.
"""
import itertools
import sympy as sp

metric_sign=(1,-1,-1,-1)
s,y,alp_mass,Q=sp.symbols('s y m_a Q',positive=True,real=True)
energy=(alp_mass**2-s-y)/(2*sp.sqrt(s))
k=sp.Matrix([sp.sqrt(s),0,0,0]); q=sp.Matrix([energy,0,0,Q]); p=-k-q
photon=sp.Matrix([0,0,1,0]); vector_polarization=sp.Matrix([0,1,0,0])
pairs=tuple(itertools.combinations(range(4),2))
components=sp.symbols('V_01 V_02 V_03 V_12 V_13 V_23')
vector=sp.zeros(4)
for value,(mu,nu) in zip(components,pairs):
    vector[mu,nu]=value;vector[nu,mu]=-value


def wedge(first,second):
    return first*second.T-second*first.T


def contractions(first_momentum,second_momentum,first_tensor,second_tensor):
    """Published VJP contractions before restricting either tensor to a wedge."""
    axial=sp.I*p
    result=[sp.S.Zero]*7
    for indices in itertools.permutations(range(4)):
        mu,nu,rho,sigma=indices; sign=sp.LeviCivita(*indices)
        for alpha in range(4):
            result[0]+=sign*first_tensor[mu,nu]*second_tensor[rho,alpha]*sp.I*metric_sign[alpha]*p[alpha]*axial[sigma]
            result[1]+=sign*first_tensor[mu,alpha]*second_tensor[rho,sigma]*sp.I*metric_sign[alpha]*p[alpha]*axial[nu]
            result[4]+=sign*sp.I*metric_sign[alpha]*first_momentum[alpha]*first_tensor[mu,nu]*second_tensor[rho,alpha]*axial[sigma]
            result[5]+=sign*sp.I*metric_sign[alpha]*first_momentum[alpha]*first_tensor[mu,alpha]*second_tensor[rho,sigma]*axial[nu]
            result[6]+=sign*sp.I*first_momentum[sigma]*first_tensor[mu,nu]*second_tensor[rho,alpha]*metric_sign[alpha]*axial[alpha]
        result[2]+=sign*first_tensor[mu,nu]*second_tensor[rho,sigma]
    return [sp.factor(sp.rem(sp.expand(value),Q**2-energy**2+y,Q)) for value in result]


first=contractions(k,q,sp.I*vector,sp.I*wedge(q,photon))
second=contractions(q,k,sp.I*wedge(q,photon),sp.I*vector)
couplings=sp.symbols('c_1:8',real=True);pair_couplings=sp.symbols('d_1:5',real=True)
vector_axial_trace,vector_mass_trace=sp.symbols('Tr_TiQUa Tr_TiQYa')
pair_axial_trace,pair_mass_trace=sp.symbols('Tr_TiTjUa Tr_TiTjYa')
resonance_mass,electric_charge,photon_response=sp.symbols('M_V e W_j')
vjp=-2/resonance_mass*sum(couplings[index]*(vector_mass_trace if index==2 else vector_axial_trace)*first[index] for index in range(7))
vvp=-2*(pair_couplings[0]*pair_axial_trace*(first[0]+second[0])+pair_couplings[1]*pair_mass_trace*(first[2]+second[2])+pair_couplings[2]*pair_axial_trace*(first[4]+second[4])+pair_couplings[3]*pair_axial_trace*(first[6]+second[6]))
production=sp.expand(-2*electric_charge*vjp+sp.sqrt(2)*electric_charge*photon_response*vvp)
production_components=sp.Matrix([sp.factor(sp.diff(production,value)) for value in components])
# L=V_mu_nu A^mu_nu includes both antisymmetric orders.
source_upper=sp.zeros(4)
for coefficient,(mu,nu) in zip(production_components,pairs):
    source_upper[mu,nu]=coefficient/(2*metric_sign[mu]*metric_sign[nu])
    source_upper[nu,mu]=-source_upper[mu,nu]
source_spin_one=source_upper.copy()
for mu,nu in pairs:
    if mu and nu:
        source_spin_one[mu,nu]=0;source_spin_one[nu,mu]=0
source_local=sp.simplify(source_upper-source_spin_one)
vector_wedge=wedge(k,vector_polarization)
wedge_substitution={symbol:vector_wedge[mu,nu] for symbol,(mu,nu) in zip(components,pairs)}
reference=sum(sp.LeviCivita(*indices)*vector_polarization[indices[0]]*photon[indices[1]]*k[indices[2]]*q[indices[3]] for indices in itertools.permutations(range(4)))
projected_numerator=sp.factor(production.subs(wedge_substitution)/reference)


def show(label,value):
    print(f'{label}:\n{sp.latex(value)}\n')

if __name__=='__main__':
    show('All incoming k+q+p',k+q+p)
    show('Independent antisymmetric components ordered 01,02,03,12,13,23',production_components)
    show('Wedge-projected numerator',projected_numerator)
    show('Production spin-one tensor source',source_spin_one)
    show('Production local tensor source',source_local)

# Contract with the full decay tensor source. In this frame the transverse
# photon polarization selects A01 and A13. No restriction is imposed on B13.
vector_mass=sp.symbols('M_i',positive=True)
decay_01,decay_13=sp.symbols('B_01 B_13')
pole_contraction=sp.factor(-4*source_upper[0,1]*decay_01/(s-vector_mass**2))
local_contraction=sp.factor(-4*source_upper[1,3]*decay_13/vector_mass**2)
full_contraction=pole_contraction+local_contraction
vjp_local_coefficients=sp.Matrix([sp.factor(sp.diff(vjp,value).diff(components[4])) for value in couplings])
vvp_local_coefficients=sp.Matrix([sp.factor(sp.diff(vvp,value).diff(components[4])) for value in pair_couplings])
if __name__=='__main__':
    show('VJP coefficient of V13 per c1..c7',vjp_local_coefficients)
    show('VVP coefficient of V13 per d1..d4',vvp_local_coefficients)
    show('Production numerator minus C01/Q',sp.factor(projected_numerator-production_components[0]/Q))
    show('Pole-projector contraction with full three-pion source',pole_contraction)
    show('Regular tensor-projector contraction with full three-pion source',local_contraction)

# Insert all five independently weighted VPPP endpoint tensors. Their
# epsilon two-form components are kept covariant rather than expanding
# individual pion four-vector components in this frame.
import derive_vector_three_pion_contacts as contacts
endpoint_coefficients=(sp.Matrix(1,5,contacts.couplings)*contacts.tensor_source_coefficients/contacts.resonance_mass).applyfunc(sp.factor)
dual_01=sp.Matrix(sp.symbols('epsilon_01_pplus_pminus epsilon_01_pplus_pzero epsilon_01_pminus_pzero'))
dual_13=sp.Matrix(sp.symbols('epsilon_13_pplus_pminus epsilon_13_pplus_pzero epsilon_13_pminus_pzero'))
# contacts derives lower tensor coefficients in V^{mu nu} B_{mu nu};
# raising 01 reverses its sign, whereas raising 13 preserves it.
endpoint_01=-(endpoint_coefficients*dual_01)[0]
endpoint_13=(endpoint_coefficients*dual_13)[0]
contact_pole_contraction=pole_contraction.subs(decay_01,endpoint_01)
contact_local_contraction=local_contraction.subs(decay_13,endpoint_13)
if __name__=='__main__':
    show('Five-coupling contact endpoint B upper 01',endpoint_01)
    show('Five-coupling contact endpoint B upper 13',endpoint_13)
    show('Full local-projector ALP contact contribution',contact_local_contraction)
