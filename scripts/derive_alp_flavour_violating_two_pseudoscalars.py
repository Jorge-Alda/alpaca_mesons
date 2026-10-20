"""Derive a -> PP flavour-changing amplitudes from the LO Alda:2026nrd action.

Run: .venv/bin/python scripts/derive_alp_flavour_violating_two_pseudoscalars.py
Tests: .venv/bin/python scripts/derive_alp_flavour_violating_two_pseudoscalars.py --test
Natural units; tree level, O(a Phi^2/f_a), no resonance or WZW operators.
Fourier convention exp(-i k x); S=1+iT. Charged outgoing states contract
with conjugate fields, hence their generators are Hermitian conjugated.
Source: Alda:2026nrd arXiv v1 Appendix A, especially C_phi[A]=[phi,A].
The companion derive_alp_phi2_lo_chipt.py uses the same source orientation.
"""
import sys
import unittest
from itertools import combinations_with_replacement
import sympy as sp

try:
    from scripts.symbolic_utils import scalar_coefficient, matrix_unit
except ModuleNotFoundError:
    from symbolic_utils import scalar_coefficient, matrix_unit


def comm(first, second):
    return first * second - second * first


def anti(first, second):
    return first * second + second * first


def show(label, expression):
    print(f"{label}:\n{sp.latex(expression)}\n")


f_pi, f_a, B_0 = sp.symbols('f_pi f_a B_0', positive=True)
a, da, x, y, dx, dy = sp.symbols('a partial_mu_a phi_i phi_j partial_mu_phi_i partial_mu_phi_j')
rho, epsilon = sp.symbols('rho epsilon')
Phi = sp.MatrixSymbol('Phi', 3, 3)
dPhi = sp.MatrixSymbol('partial^mu_Phi', 3, 3)
cV = sp.MatrixSymbol('c_q^V', 3, 3)
cS = sp.MatrixSymbol('c_q^S', 3, 3)
mq = sp.MatrixSymbol('m_q', 3, 3)
Ti = sp.MatrixSymbol('lambda_i', 3, 3)
Tj = sp.MatrixSymbol('lambda_j', 3, 3)
phi = rho * Phi / (2 * f_pi)
u_plain = -rho * dPhi / f_pi
# Source definition is C_phi[A]=[phi,A], NOT [A,phi].
u_vector = -2 * sp.I * epsilon * da / f_a * comm(phi, cV)
kinetic_raw = sp.Trace(scalar_coefficient(scalar_coefficient(sp.expand(
    f_pi**2 / 4 * (u_plain*u_vector + u_vector*u_plain)), epsilon, 1), rho, 2))
chi_scalar = 2*sp.I*B_0*epsilon*a/f_a*comm(cS, mq)
scalar_raw = sp.Trace(scalar_coefficient(scalar_coefficient(sp.expand(
    f_pi**2/4*(2*chi_scalar - anti(phi, anti(phi, chi_scalar)))), epsilon, 1), rho, 2))
# Cyclicity reductions of the raw traces. Equality is tested on full generic
# component matrices, rather than taken as a separate physical input.
kinetic = -sp.I*da/(2*f_a)*sp.Trace(dPhi*comm(cV, Phi))
scalar = -sp.I*B_0*a/(2*f_a)*sp.Trace(Phi**2*comm(cS, mq))
component_kinetic = sp.expand(-sp.I*da/(2*f_a)*(dx*Ti+dy*Tj)*comm(cV,x*Ti+y*Tj))
component_scalar = sp.expand(-sp.I*B_0*a/(2*f_a)*(x**2*Ti**2+x*y*(Ti*Tj+Tj*Ti)+y**2*Tj**2)*comm(cS,mq))
# The quadratic product is distributed explicitly because SymPy MatPow does
# not expand a power of MatAdd under expand(); matrix algebra remains SymPy.
# Explicit coefficient differentiation performs field extraction before Fourier
# replacement. There are two derivative placements for distinct fields.
coeff_x_dy = sp.Trace(scalar_coefficient(scalar_coefficient(scalar_coefficient(component_kinetic, da, 1), x, 1), dy, 1))
coeff_y_dx = sp.Trace(scalar_coefficient(scalar_coefficient(scalar_coefficient(component_kinetic, da, 1), y, 1), dx, 1))
coeff_a_x_y = sp.Trace(scalar_coefficient(scalar_coefficient(scalar_coefficient(component_scalar, a, 1), x, 1), y, 1))
ka_ki, ka_kj = sp.symbols('k_a_dot_k_i k_a_dot_k_j', real=True)
vertex_raw = sp.expand(sp.I*((-sp.I)**2*(ka_kj*coeff_x_dy+ka_ki*coeff_y_dx)+coeff_a_x_y))
# Trace cyclicity gives these generator expressions from extracted coefficients.
vector_trace = sp.Trace(cV*comm(Ti,Tj))
scalar_trace = sp.Trace(anti(Ti,Tj)*comm(cS,mq))
vertex = (-vector_trace*(ka_kj-ka_ki)+B_0*scalar_trace)/(2*f_a)
m_a, m_i, m_j = sp.symbols('m_a m_i m_j', positive=True)
# k_a=p_a, k_i=-p_i, k_j=-p_j, all incoming at the rule.
kinematics = {ka_ki:-(m_a**2+m_i**2-m_j**2)/2,
              ka_kj:-(m_a**2+m_j**2-m_i**2)/2}
i_amplitude = sp.expand(vertex.subs(kinematics))


def E(row, column):
    return matrix_unit(row, column)


theta_P = sp.Symbol('theta_P', real=True)
eta8 = sp.diag(1,1,-2)/sp.sqrt(3)
eta0 = sp.sqrt(sp.Rational(2,3))*sp.eye(3)
generators = {
    'pi+':sp.sqrt(2)*E(0,1), 'pi-':sp.sqrt(2)*E(1,0),
    'pi0':sp.diag(1,-1,0), 'K+':sp.sqrt(2)*E(0,2),
    'K-':sp.sqrt(2)*E(2,0), 'K0':sp.sqrt(2)*E(1,2),
    'K0bar':sp.sqrt(2)*E(2,1),
    'eta':sp.cos(theta_P)*eta8-sp.sin(theta_P)*eta0,
    'eta_prime':sp.sin(theta_P)*eta8+sp.cos(theta_P)*eta0,
}
charges = dict(zip(generators, (1,-1,0,1,-1,0,0,0,0)))
conjugate_names = {'pi+':'pi-', 'pi-':'pi+', 'K+':'K-', 'K-':'K+',
                   'K0':'K0bar', 'K0bar':'K0', 'pi0':'pi0',
                   'eta':'eta', 'eta_prime':'eta_prime'}
m_u,m_d,m_s = sp.symbols('m_u m_d m_s', real=True)
mq3 = sp.diag(m_u,m_d,m_s)
# c_ds is row d, column s; Hermiticity fixes c_sd=c_ds^*.
cV_ds,cS_ds = sp.symbols('c_ds^V c_ds^S')
cVS_ds = sp.Symbol('c_ds^{VS}')
cV3 = cV_ds*E(1,2)+sp.conjugate(cV_ds)*E(2,1)
cS3 = cS_ds*E(1,2)+sp.conjugate(cS_ds)*E(2,1)
masses = {name:sp.Symbol('m_'+name.replace('+','_plus').replace('-','_minus'), positive=True)
          for name in generators}
for name in ('pi-','K-','K0bar'):
    masses[name] = masses[conjugate_names[name]]


def extracted_channel(first, second, source_v=cV3, source_s=cS3):
    """Extract the action with outgoing LSZ matrices, then construct i M."""
    first_matrix = generators[first].H
    second_matrix = generators[second].H
    field = x*first_matrix+y*second_matrix
    derivative = dx*first_matrix+dy*second_matrix
    action_v = sp.expand(-sp.I*da/(2*f_a)*sp.trace(derivative*comm(source_v,field)))
    action_s = sp.expand(-sp.I*B_0*a/(2*f_a)*sp.trace(field**2*comm(source_s,mq3)))
    result = sp.I*((-sp.I)**2*(ka_kj*sp.diff(action_v,da,x,dy)
                  +ka_ki*sp.diff(action_v,da,y,dx))+sp.diff(action_s,a,x,y))
    return sp.simplify(sp.expand(result.subs(kinematics).subs(
        {m_i:masses[first],m_j:masses[second]})))


channels = [('K+','pi-'),('K0','pi0'),('K0','eta'),('K0','eta_prime')]
# LO neutral mass map, valid also with isospin violation. For displaying a
# single eta angle below set m_u=m_d; arbitrary theta and masses cannot be
# imposed independently while claiming the LO EOM identity.
m_0_sq = sp.Symbol('m_0^2', real=True)


def mass_map(matrix):
    return B_0*anti(mq3,matrix)+m_0_sq*sp.trace(matrix)*sp.eye(3)/3


neutral_mass_matrix = sp.simplify(sp.Matrix([
    [sp.trace(left*mass_map(right))/2 for right in (eta8,eta0)]
    for left in (eta8,eta0)]).subs(m_u,m_d))
rotation = sp.Matrix([[sp.cos(theta_P),sp.sin(theta_P)],[-sp.sin(theta_P),sp.cos(theta_P)]])
rotated_neutral_mass = sp.simplify(rotation.T*neutral_mass_matrix*rotation)
# The following exact matrix identity is the EOM reason cV-cS enters. The
# singlet anomaly drops out of commutators, without neglecting m_0^2.
generic_first = sp.Matrix(3,3,sp.symbols('t_0:9'))
generic_second = sp.Matrix(3,3,sp.symbols('r_0:9'))
eom_identity_left = comm(mass_map(generic_first),generic_second)-comm(generic_first,mass_map(generic_second))
eom_identity_right = B_0*comm(mq3,anti(generic_first,generic_second))


def run_derivation():
    print('LO even-intrinsic-parity action; natural units; O(a Phi^2/f_a); Fourier exp(-ikx).')
    show('Alda:2026nrd commutator map C_phi[c_q^V]=[phi,c_q^V]',comm(phi,cV))
    show('Raw kinetic extraction',kinetic_raw)
    show('Kinetic interaction after trace cyclicity',kinetic)
    show('Raw scalar extraction',scalar_raw)
    show('Scalar interaction after trace cyclicity',scalar)
    show('Field-extracted vector coefficient phi_i partial phi_j',coeff_x_dy)
    show('Field-extracted scalar coefficient a phi_i phi_j',coeff_a_x_y)
    show('All-incoming vertex i M before trace cyclicity',vertex_raw)
    show('All-incoming vertex i M after trace cyclicity',vertex)
    show('Decay momentum replacement k_a dot k_i',kinematics[ka_ki])
    show('Decay momentum replacement k_a dot k_j',kinematics[ka_kj])
    show('Physical i M; replace lambda_i,j by outgoing conjugate generators',i_amplitude)
    show('Hermitian d-s vector source in u,d,s order',cV3)
    show('Hermitian d-s scalar source in u,d,s order',cS3)
    neutral_diagonal = sp.diag(*sp.symbols(r'(lambda_P)_{uu} (lambda_P)_{dd} (lambda_P)_{ss}', real=True))
    generators['P'] = neutral_diagonal
    masses['P'] = sp.Symbol('m_P',positive=True)
    show('General neutral outgoing generator, including isospin violation', neutral_diagonal)
    show('i M(a -> K0 P), general neutral generator',extracted_channel('K0','P'))
    show('i M(a -> K0 P), scalar-free source cVS',extracted_channel('K0','P',cVS_ds*E(1,2),sp.zeros(3)))
    del generators['P']
    del masses['P']
    for first,second in channels:
        show(f'Outgoing generators {first}, {second}',sp.Tuple(generators[first].H,generators[second].H))
        show(f'i M(a -> {first} {second}), general scalar and vector sources',extracted_channel(first,second))
        show(f'i M(a -> {first} {second}), cS=0 and cV=cVS',extracted_channel(first,second, cVS_ds*E(1,2),sp.zeros(3)))
        conjugate_pair=tuple(conjugate_names[name] for name in (first,second))
        show(f'i M(a -> {conjugate_pair[0]} {conjugate_pair[1]})',extracted_channel(*conjugate_pair))
    show('LO eta8 eta0 squared-mass matrix, m_u=m_d',neutral_mass_matrix)
    show('Rotated neutral squared-mass matrix; off-diagonal entry must vanish',rotated_neutral_mass)
    show('Exact mass-map commutator identity, left minus right',sp.simplify(eom_identity_left-eom_identity_right))
    show('Mass-eigenstate trace identity', sp.Eq((m_i**2-m_j**2)*sp.Trace(cS*comm(sp.Adjoint(Ti),sp.Adjoint(Tj))), B_0*sp.Trace(anti(sp.Adjoint(Ti),sp.Adjoint(Tj))*comm(cS,mq))))
    show('Physical i M after LO EOM, c_qVS=c_qV-c_qS', -(m_i**2-m_j**2)/(2*f_a)*sp.Trace((cV-cS)*comm(sp.Adjoint(Ti),sp.Adjoint(Tj))))
    amplitude_squared = sp.Abs(sp.Symbol(r'\mathcal{M}'))**2
    phase_space = sp.sqrt((m_a**2-(m_i+m_j)**2)*(m_a**2-(m_i-m_j)**2))
    show('Distinct-state two-body width, m_a>m_i+m_j', phase_space*amplitude_squared/(16*sp.pi*m_a**3))


class DerivationTests(unittest.TestCase):
    def test_outgoing_conjugation(self):
        for name,matrix in generators.items():
            self.assertEqual(matrix.H,generators[conjugate_names[name]])
        self.assertTrue(extracted_channel('K+','pi-').has(cV_ds))
        self.assertFalse(extracted_channel('K+','pi-').has(sp.conjugate(cV_ds)))

    def test_raw_trace_reductions(self):
        full = sp.Matrix(3,3,sp.symbols('p_0:9'))
        deriv = sp.Matrix(3,3,sp.symbols('d_0:9'))
        source = sp.Matrix(3,3,sp.symbols('c_0:9'))
        replacements={Phi:full,dPhi:deriv,cV:source,cS:source,mq:mq3}
        for raw,reduced in ((kinetic_raw,kinetic),(scalar_raw,scalar),(vertex_raw,vertex)):
            substitutions=dict(replacements, **{})
            substitutions.update({Ti:generators['K-'],Tj:generators['pi+']})
            self.assertEqual(sp.expand((raw-reduced).doit().subs(substitutions).doit()),0)

    def test_all_flavour_changing_channels(self):
        found=set()
        for first,second in combinations_with_replacement(generators,2):
            if charges[first]+charges[second]==0 and extracted_channel(first,second)!=0:
                found.add(frozenset((first,second)))
        expected={frozenset(pair) for pair in channels}
        expected|={frozenset(conjugate_names[name] for name in pair) for pair in channels}
        self.assertEqual(found,expected)

    def test_hermitian_conjugate_amplitudes(self):
        for pair in channels:
            conjugate_pair=tuple(conjugate_names[name] for name in pair)
            # i M of conjugate channel is minus the complex conjugate of i M;
            # equivalently M_conjugate=M^* in this real-field phase convention.
            self.assertEqual(sp.simplify(extracted_channel(*conjugate_pair)+sp.conjugate(extracted_channel(*pair))),0)

    def test_exact_mass_map_identity(self):
        self.assertEqual(sp.expand(eom_identity_left-eom_identity_right),sp.zeros(3))

    def test_lo_basis_invariance_charged(self):
        raw=extracted_channel('K+','pi-').subs({masses['K+']**2:B_0*(m_u+m_s),masses['pi-']**2:B_0*(m_u+m_d)})
        invariant=extracted_channel('K+','pi-',(cV_ds-cS_ds)*E(1,2),sp.zeros(3)).subs({masses['K+']**2:B_0*(m_u+m_s),masses['pi-']**2:B_0*(m_u+m_d)})
        self.assertEqual(sp.simplify(raw-invariant),0)

    def test_neutral_eigenstate_basis_invariance(self):
        # Exact nondegenerate LO eigenvectors with nonzero anomaly, not a
        # separately chosen phenomenological mixing angle and mass pair.
        source_s = cS_ds*E(1,2)
        source_v = cV_ds*E(1,2)
        parameters={B_0:1,m_u:1,m_d:1,m_s:4,m_0_sq:6}
        matrix=neutral_mass_matrix.subs(parameters)
        for eigenvalue,multiplicity,eigenvectors in matrix.eigenvects():
            coefficients=eigenvectors[0]
            neutral=coefficients[0]*eta8+coefficients[1]*eta0
            first=generators['K0'].H
            delta=5-eigenvalue
            raw=(-delta*sp.trace(source_v*comm(first,neutral))
                 +sp.trace(anti(first,neutral)*comm(source_s,mq3.subs(parameters))))/(2*f_a)
            invariant=-delta*sp.trace((source_v-source_s)*comm(first,neutral))/(2*f_a)
            self.assertEqual(sp.simplify(raw-invariant),0)

    def test_diagonal_equal_mass_zero(self):
        diagonal=sp.diag(*sp.symbols('c_u c_d c_s',real=True))
        for pair in (('pi+','pi-'),('K+','K-'),('K0','K0bar')):
            self.assertEqual(extracted_channel(*pair,diagonal,diagonal),0)


if __name__=='__main__':
    if '--test' in sys.argv:
        unittest.main(argv=[sys.argv[0]])
    else:
        run_derivation()
