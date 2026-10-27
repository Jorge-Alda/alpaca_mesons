"""Derive the local a V_i -> pi+ pi- VPPP vertex before photon mixing.

Run .venv/bin/python scripts/derive_alp_vppp_isovector_tff.py [--test].
Natural units; epsilon^0123=+1; metric (+---); tree level; one ALP.
U_a,Y_a are supplied no-mixing sources, chi_minus=-i Y. No equal-light-
quark-spurion assumption is used. Charged pion masses are equal by charge
conjugation; the intermediate expressions keep their on-shell mass separate
from chi. The final LO wave uses m_pi^2=(chi_u+chi_d)/2 through a trace of
the charged generators and chi; neither chi_u nor chi_d is equated to it.
The finite Wick sum retains derivative slots; SymPy evaluates all traces,
substitutions, scalar products, angular projections and displayed algebra.
"""
import argparse
import itertools
import unittest
from functools import lru_cache
import sympy as sp


def operator_coefficients(vector, derivative_sources, mass_sources, scalar_products):
    """Expand the five literal operators over every pseudoscalar assignment.

    The result divides out i epsilon(p+,p-,-P,epsilon_i). The vector tensor
    here is q wedge epsilon_i, without 1/M_i. q=p++p--P=-k. Native SymPy
    matrix differentiation does not retain the separate momentum slots;
    only this finite assignment sum is written explicitly.
    """
    basis=tuple(sp.eye(4)[:,index] for index in range(4))
    q=sum(basis[:3],sp.zeros(4,1))
    result=sp.zeros(1,5)
    for first,second,third in itertools.permutations(range(3)):
        first_u,second_u,third_u=(derivative_sources[index] for index in (first,second,third))
        def epsilon_ratio(first_index,second_index):
            return sp.det(sp.Matrix.hstack(q,basis[3],basis[first_index],basis[second_index]))
        difference=(sp.Trace(vector*first_u*second_u*third_u)
                    -sp.Trace(vector*third_u*second_u*first_u))
        result[0,0]+=4*epsilon_ratio(first,third)*scalar_products[first,second]*difference
        result[0,1]+=4*epsilon_ratio(first,second)*scalar_products[first,third]*difference
        result[0,2]+=4*epsilon_ratio(second,third)*scalar_products[first,second]*difference
        result[0,3]+=2*epsilon_ratio(first,second)*(sp.Trace(vector*first_u*second_u*mass_sources[third])
                                                         +sp.Trace(first_u*second_u*vector*mass_sources[third]))
        result[0,4]+=2*epsilon_ratio(first,second)*sp.Trace(first_u*vector*second_u*mass_sources[third])
    return result


@lru_cache(maxsize=1)
def derive():
    """Construct vector vertex with abstract flavour traces, then mix photon."""
    s,t,u,y,ma2,mpi2=sp.symbols('s t u y m_a_squared m_pi_squared',real=True)
    fpi,fa,mv,mi=sp.symbols('f_pi f_a M_V M_i',positive=True)
    gs=sp.symbols('g_1:6',real=True)
    vector,first_pion,second_pion,derivative_alp,mass_alp,chi,response=(
        sp.MatrixSymbol(name,3,3) for name in
        ('T_i','lambda_+^dagger','lambda_-^dagger','U_a','Y_a','chi','W_gamma'))
    derivative=(-first_pion/fpi,-second_pion/fpi,derivative_alp)
    mass=((first_pion*chi+chi*first_pion)/fpi,
          (second_pion*chi+chi*second_pion)/fpi,mass_alp)
    dot=sp.Matrix([[mpi2,(s-2*mpi2)/2,(t-mpi2-ma2)/2],
                   [(s-2*mpi2)/2,mpi2,(u-mpi2-ma2)/2],
                   [(t-mpi2-ma2)/2,(u-mpi2-ma2)/2,ma2]])
    coefficients=operator_coefficients(vector,derivative,mass,dot)
    coefficients=coefficients.applyfunc(lambda value:sp.expand(value.subs(u,ma2+y+2*mpi2-s-t)))
    vector_coefficients=coefficients/(mv*mi)
    vector_wave=(vector_coefficients*sp.Matrix(gs))[0]
    # External-vector coefficient multiplies i epsilon(p+,p-,-P,eps_i).
    # Photon response F_gamma=i k wedge eps=-i q wedge eps cancels that i.
    # Crossing yields epsilon(p+,p-,-P,eps)=epsilon(eps,p+,p-,k).
    photon_coefficients=vector_coefficients.subs(vector,sp.sqrt(2)*mi*response).applyfunc(sp.expand)
    photon_wave=(photon_coefficients*sp.Matrix(gs))[0]
    angle=sp.Symbol('z',real=True)
    span=sp.Symbol('B_pi',real=True)
    projected=sp.integrate(sp.Rational(3,4)*(1-angle**2)*photon_wave.subs(
        t,mpi2+(ma2-s+y)/2-span*angle),(angle,-1,1))
    # Explicit entries serve only to test the abstract trace expression.
    chi_entries=sp.symbols('chi_u chi_d chi_s',real=True)
    vector_entries=sp.symbols('t_u t_d t_s',real=True)
    axial_entries=sp.symbols('c_u^A c_d^A c_s^A',real=True)
    pseudoscalar_entries=sp.symbols('c_u^P c_d^P c_s^P',real=True)
    response_entries=sp.symbols('w_u w_d w_s',real=True)
    first=sp.zeros(3);second=sp.zeros(3)
    first[1,0]=sp.sqrt(2);second[0,1]=sp.sqrt(2)
    chi_component=sp.diag(*chi_entries)
    derivative_component=2*sp.diag(*axial_entries)/fa
    mass_component=-2*(sp.diag(*pseudoscalar_entries)*chi_component
                       +chi_component*sp.diag(*pseudoscalar_entries))/fa
    component_replacements={vector:sp.diag(*vector_entries),first_pion:first,
                            second_pion:second,derivative_alp:derivative_component,
                            mass_alp:mass_component,chi:chi_component,
                            response:sp.diag(*response_entries)}
    component_vector=vector_wave.subs(component_replacements).doit().expand()
    def commutator(first,second):
        return first*second-second*first
    def anticommutator(first,second):
        return first*second+second*first
    # Cyclic traces and commuting neutral matrices reduce the derivative
    # sector to these two charged-generator words. Determine their scalar
    # coefficients from the already calculated Wick sum, rather than
    # transcribing a component polynomial.
    derivative_trace_first=sp.Trace(vector*derivative_alp*commutator(first_pion,second_pion))
    derivative_trace_second=sp.Trace(vector*(first_pion*derivative_alp*second_pion
                                                       -second_pion*derivative_alp*first_pion))
    scalar_first,scalar_second=sp.symbols('coefficient_first coefficient_second')
    compact_operator_terms=[]
    common_factor=2/(mv*mi*fpi**2)
    for index in range(3):
        ansatz=common_factor*(scalar_first*derivative_trace_first+scalar_second*derivative_trace_second)
        difference=(ansatz-vector_coefficients[0,index]).subs(component_replacements).doit().expand()
        polynomial=sp.Poly(difference,*vector_entries,*axial_entries)
        solution=sp.solve(polynomial.coeffs(),(scalar_first,scalar_second),dict=True)[0]
        compact_operator_terms.append((ansatz/common_factor).subs(solution))
    # Group the six mass-source assignments before trace expansion. This
    # retains the literal anticommutators with the full unequal-mass chi.
    mass_trace_fourth=sp.Trace(
        anticommutator(vector,commutator(first_pion,second_pion))*mass_alp
        +anticommutator(vector,commutator(first_pion,derivative_alp))*anticommutator(second_pion,chi)
        -anticommutator(vector,commutator(second_pion,derivative_alp))*anticommutator(first_pion,chi))
    mass_trace_fifth=sp.Trace(
        (first_pion*vector*second_pion-second_pion*vector*first_pion)*mass_alp
        +(first_pion*vector*derivative_alp-derivative_alp*vector*first_pion)*anticommutator(second_pion,chi)
        -(second_pion*vector*derivative_alp-derivative_alp*vector*second_pion)*anticommutator(first_pion,chi))
    compact_operator_terms.extend((mass_trace_fourth,mass_trace_fifth))
    compact_bracket=((gs[0]-gs[2])*compact_operator_terms[0]
                     +gs[1]*compact_operator_terms[1]+gs[3]*mass_trace_fourth+gs[4]*mass_trace_fifth)
    compact_vector=common_factor*compact_bracket
    compact_photon=2*sp.sqrt(2)/(mv*fpi**2)*compact_bracket.subs(vector,response)
    charged_pion_mass_squared=sp.Trace((first_pion*second_pion+second_pion*first_pion)*chi)/4
    lo_projected=projected.subs(mpi2,charged_pion_mass_squared)
    return locals()


class VPPPTests(unittest.TestCase):
    """Tests retain unequal chi_u,chi_d throughout the flavour comparison."""
    def test_abstract_trace_matches_component_assignments(self):
        d=derive()
        derivative=(-d['first']/d['fpi'],-d['second']/d['fpi'],d['derivative_component'])
        chi=d['chi_component']
        mass=((d['first']*chi+chi*d['first'])/d['fpi'],
              (d['second']*chi+chi*d['second'])/d['fpi'],d['mass_component'])
        component=operator_coefficients(sp.diag(*d['vector_entries']),derivative,mass,d['dot'])
        component=(component*sp.Matrix(d['gs']))[0]/(d['mv']*d['mi'])
        component=component.doit().subs(d['u'],d['ma2']+d['y']+2*d['mpi2']-d['s']-d['t'])
        self.assertEqual(sp.simplify(component-d['component_vector']),0)
        # Unequal light spurions are actually present, rather than silently
        # replaced by an isospin-symmetric pion mass.
        self.assertNotEqual(sp.diff(d['component_vector'],d['chi_entries'][0]),0)
        self.assertNotEqual(sp.diff(d['component_vector'],d['chi_entries'][1]),0)

    def test_photon_response_is_linear_vector_sum(self):
        d=derive()
        first_residue,second_residue,first_mass,second_mass=sp.symbols('R_1 R_2 M_1 M_2',nonzero=True)
        first_generator,second_generator=(sp.MatrixSymbol(name,3,3) for name in ('T_1','T_2'))
        response=first_residue*first_generator/(first_mass**2-d['y'])+second_residue*second_generator/(second_mass**2-d['y'])
        directly_mixed=sp.sqrt(2)*d['mi']*(first_residue*d['vector_wave'].subs(d['vector'],first_generator)/(first_mass**2-d['y'])
                         +second_residue*d['vector_wave'].subs(d['vector'],second_generator)/(second_mass**2-d['y']))
        difference=sp.expand(d['photon_wave'].subs(d['response'],response)-directly_mixed).doit().expand()
        self.assertEqual(sp.simplify(difference),0)

    def test_charged_isospin_limit_and_pion_specialization(self):
        d=derive()
        source={d['axial_entries'][0]:-d['fa']/(2*d['fpi']),
                d['axial_entries'][1]:d['fa']/(2*d['fpi']),d['axial_entries'][2]:0,
                d['pseudoscalar_entries'][0]:-d['fa']/(2*d['fpi']),
                d['pseudoscalar_entries'][1]:d['fa']/(2*d['fpi']),d['pseudoscalar_entries'][2]:0,
                d['chi_entries'][0]:d['mpi2'],d['chi_entries'][1]:d['mpi2'],
                d['ma2']:d['mpi2'],d['vector_entries'][1]:d['vector_entries'][0]}
        g=d['gs']
        expected=-16*d['vector_entries'][0]/(d['mi']*d['fpi']**3*d['mv'])*(
            (g[0]-g[1]-g[2])*(d['y']-3*d['mpi2'])+3*d['mpi2']*(2*g[3]+g[4]))
        self.assertEqual(sp.simplify(d['component_vector'].subs(source)-expected),0)

    def test_angular_projection_and_strange_source(self):
        d=derive()
        self.assertEqual(sp.simplify((d['projected']-d['photon_wave']).subs(d['component_replacements']).doit()),0)
        self.assertEqual(sp.diff(d['component_vector'],d['t']),0)
        for symbol in (d['axial_entries'][2],d['pseudoscalar_entries'][2],d['chi_entries'][2]):
            self.assertEqual(sp.diff(d['component_vector'],symbol),0)
        self.assertEqual(sp.simplify(d['charged_pion_mass_squared'].subs(d['component_replacements']).doit()
                                    -(d['chi_entries'][0]+d['chi_entries'][1])/2),0)

    def test_compact_trace_formula(self):
        d=derive()
        difference=(d['compact_vector']-d['vector_wave']).subs(d['component_replacements']).doit()
        self.assertEqual(sp.simplify(difference),0)
        photon_difference=(d['compact_photon']-d['photon_wave']).subs(d['component_replacements']).doit()
        self.assertEqual(sp.simplify(photon_difference),0)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--test',action='store_true')
    arguments=parser.parse_args()
    if arguments.test:
        unittest.main(argv=['vppp'])
        return
    d=derive()
    print('Tree level; diagonal neutral T_i and no-mixing ALP sources; arbitrary diagonal chi; natural units.\n')
    for label,value in [('Derivative sources',sp.Tuple(*d['derivative'])),
                        ('Mass sources with chi_minus=-i Y',sp.Tuple(*d['mass'])),
                        ('Compact local vector scalar amplitude (divide i epsilon)',d['compact_vector']),
                        ('Compact on-shell external-vector scalar (y=M_i squared)',d['compact_vector'].subs(d['y'],d['mi']**2)),
                        ('Compact photon-mixed charge-stripped current and P wave',d['compact_photon']),
                        ('LO charged-pion mass squared, equal to (chi_u+chi_d)/2',d['charged_pion_mass_squared']),
                        ]:
        print(label+':\n'+sp.latex(value)+'\n')


if __name__=='__main__':
    main()
