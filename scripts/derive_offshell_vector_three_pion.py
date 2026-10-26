"""Derive the transverse sequential-rho contribution to V_i^* -> 3pi.

Natural units, (+---), epsilon^{0123}=+1. This is NOT the complete
three-pion amplitude. The tensor source is explicitly projected as
V^{mu nu}=(q^mu epsilon^nu-q^nu epsilon^mu)/M_i even for q² != M_i².
The fixed mass normalization defines an off-shell convention, not an LSZ
state. Only the single-trace d_i and G_V terms, equal pion masses, and
isospin-symmetric diagonal generators are retained. Local V3pi operators
and pion rescattering remain additional inputs.
The internal rho denominator uses the retained octet mass M_V of
Ruiz-Femenia:2003jdx, not a separately substituted physical rho mass.

Run: .venv/bin/python scripts/derive_offshell_vector_three_pion.py
"""
from itertools import permutations
import unittest
import sympy as sp

vector_energy, pair_energy, pair_momentum, pion_transverse_momentum = sp.symbols(
    'E_V E_pair p_pair p_perp', nonzero=True, real=True)
vector_mass, vector_virtuality, pair_virtuality, pion_mass_squared = sp.symbols(
    'M_i s s_pair m_pi_squared', positive=True)
d1, d2, d3, d4, gv, decay_constant, octet_vector_mass = sp.symbols(
    'd_1 d_2 d_3 d_4 G_V F M_V', real=True)
metric = (1, -1, -1, -1)
index_permutations = tuple(permutations(range(4)))

def wedge(first, second):
    """Upper-index antisymmetric tensor from two upper-index vectors."""
    return first * second.T - second * first.T

def show(label, expression):
    print(f'{label}:\n{sp.latex(expression)}\n')

def derive_contractions(polarization=None):
    """Contract both vector placements in each VVP operator explicitly.

    The equal-mass VPP current p1 wedge p2 belongs entirely to the
    propagating tensor subspace: r wedge (r contraction current)=r² current.
    Thus the full antisymmetric tensor propagator gives its pole part;
    its nonpropagating projector annihilates this current. Overall Fourier
    factors are restored after these real Lorentz contractions.
    """
    q = sp.Matrix([vector_energy, 0, 0, 0])
    r = sp.Matrix([pair_energy, 0, 0, pair_momentum])
    spectator = q-r
    first_pion = sp.Matrix([pair_energy/2, pion_transverse_momentum, 0, pair_momentum/2])
    second_pion = r-first_pion
    if polarization is None:
        polarization = sp.Matrix([0, 0, 1, 0])
    external_tensor = wedge(q, polarization)/vector_mass
    pair_tensor = wedge(first_pion, second_pion)
    reference = -sp.det(sp.Matrix.hstack(first_pion, second_pion, spectator, polarization))

    def contraction(first, second, operator, differentiated_momentum=None):
        terms = []
        for mu, nu, rho, sigma in index_permutations:
            epsilon = -sp.LeviCivita(mu, nu, rho, sigma)
            if operator == 2:
                terms.append(epsilon*first[mu,nu]*second[rho,sigma])
                continue
            for alpha in range(4):
                derivative = {1: spectator[alpha]*spectator[sigma],
                              3: (differentiated_momentum[alpha]*spectator[sigma]
                                  if differentiated_momentum is not None else 0),
                              4: (spectator[alpha]*differentiated_momentum[sigma]
                                  if differentiated_momentum is not None else 0)}[operator]
                terms.append(epsilon*first[mu,nu]*second[rho,alpha]*metric[alpha]*derivative)
        return sp.factor(sum(terms))

    contractions = {}
    for operator in (1, 2):
        contractions[operator] = sp.factor(contraction(external_tensor,pair_tensor,operator)
                                          +contraction(pair_tensor,external_tensor,operator))
    for operator in (3, 4):
        contractions[operator] = sp.factor(contraction(external_tensor,pair_tensor,operator,q)
                                          -contraction(pair_tensor,external_tensor,operator,r))
    return contractions, reference

def derive():
    contractions, reference = derive_contractions()
    substitutions = {pair_momentum**2: pair_energy**2-pair_virtuality}
    energy_substitution = {pair_energy: (vector_energy**2+pair_virtuality-pion_mass_squared)/(2*vector_energy)}
    coefficients = {}
    for operator, expression in contractions.items():
        ratio = sp.factor(expression/reference)
        ratio = sp.factor(sp.expand(ratio).subs(substitutions).subs(energy_substitution))
        coefficients[operator] = sp.factor(ratio.subs(vector_energy**2,vector_virtuality))
    return coefficients

def derive_scalar_amplitude():
    """Restore Fourier factors and the normalized even-parity vertex.

    O1 has two outgoing pion derivatives; O3 has one incoming vector
    derivative and one outgoing pion derivative, giving opposite signs.
    i chi_- gives twice m_pi² relative to the O1 linear vielbein factor.
    The common 4 G_V/F³ combines the linear vielbein, VPP coefficient,
    and the normalized omega flavor trace in Ruiz-Femenia:2003jdx.
    """
    coefficients = derive()
    numerator = sp.factor(vector_mass/2*(d1*coefficients[1]
                      +2*pion_mass_squared*d2*coefficients[2]-d3*coefficients[3]))
    pair_invariants = sp.symbols('s_12 s_13 s_23', real=True)
    scalar = 8*gv/(vector_mass*decay_constant**3)*sum(
        numerator.subs(pair_virtuality,invariant)/(octet_vector_mass**2-invariant)
        for invariant in pair_invariants)
    return numerator, scalar

def derive_flavor_factor():
    """Normalize the isoscalar single-trace flavor factor to ideal omega.

    Explicit components are appropriate here because the supplied pion
    charge states and neutral-generator convention fix a component basis.
    An isoscalar generator is diag(t_n,t_n,t_s); no strange pion entry
    appears in this retained single-trace tree contribution.
    """
    nonstrange_entry, strange_entry = sp.symbols('t_n t_s', real=True)
    generator = sp.diag(nonstrange_entry,nonstrange_entry,strange_entry)
    neutral_rho = sp.diag(1,-1,0)/sp.sqrt(2)
    neutral_pion = neutral_rho
    charged_pion_plus = sp.zeros(3); charged_pion_plus[0,1]=1
    charged_pion_minus = charged_pion_plus.T
    trace = sp.trace((generator*neutral_rho+neutral_rho*generator)*neutral_pion)
    trace *= sp.trace(neutral_rho*(charged_pion_plus*charged_pion_minus
                                 -charged_pion_minus*charged_pion_plus))
    omega_trace = trace.subs(nonstrange_entry,1/sp.sqrt(2))
    return sp.simplify(trace/omega_trace)

class TestOffshellVectorThreePion(unittest.TestCase):
    """Independent source, Ward, projector, and flavor checks."""
    def test_omega_on_shell_source_formula(self):
        numerator,_ = derive_scalar_amplitude()
        expected = pion_mass_squared*(d1+8*d2-d3)+(vector_mass**2+pair_virtuality)*d3
        self.assertEqual(sp.expand(numerator.subs(vector_virtuality,vector_mass**2)-expected),0)

    def test_transversality(self):
        q = sp.Matrix([vector_energy,0,0,0])
        contractions,_ = derive_contractions(q)
        self.assertTrue(all(expression == 0 for expression in contractions.values()))

    def test_pair_current_tensor_projector(self):
        r = sp.Matrix([pair_energy,0,0,pair_momentum])
        first = sp.Matrix([pair_energy/2,pion_transverse_momentum,0,pair_momentum/2])
        current = wedge(first,r-first)
        r_lower = sp.diag(*metric)*r
        contracted = current.T*r_lower
        longitudinal_numerator = wedge(r,contracted)
        self.assertEqual(sp.simplify(longitudinal_numerator-(pair_energy**2-pair_momentum**2)*current),sp.zeros(4))

    def test_ideal_flavor_limits(self):
        factor = derive_flavor_factor()
        entry = next(symbol for symbol in factor.free_symbols if symbol.name=='t_n')
        self.assertEqual(factor.subs(entry,1/sp.sqrt(2)),1)
        self.assertEqual(factor.subs(entry,0),0)

if __name__ == '__main__':
    for operator, expression in derive().items():
        show(f'Lorentz contraction O{operator} / epsilon(p1,p2,p3,epsilon)', expression)
    numerator, scalar = derive_scalar_amplitude()
    show('Derived numerator for each rho pair',numerator)
    show('Sequential-rho scalar with fixed M_i tensor-source normalization',scalar)
    show('Isoscalar single-trace flavor multiplier relative to ideal omega',derive_flavor_factor())
    print('Caveats: equal pion masses; no direct V3pi contacts; no rescattering; fixed-M_i source projection.')
    unittest.main(argv=['derive_offshell_vector_three_pion'],exit=True)
