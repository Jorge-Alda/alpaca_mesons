"""Derive the five single-trace VPPP charged isoscalar contact coefficients.

Inputs: memory/rcht/odd-parity/vector-three-pion-contact.md and the five
GomezDumm:2010ra Eq. (13) orderings already used by
scripts/derive_wzw_gamma_three_pion_and_vector_contact.py.
Natural units, metric +---, epsilon^0123=+1, U=exp(i Phi/f_pi).
Outgoing pion derivatives give +i p; the incoming vector has momentum
q=p_++p_-+p_0 and tensor (q^mu eps^nu-q^nu eps^mu)/M_i.
Output coefficients multiply epsilon_lower(p_+,p_-,p_0,eps) in L,
before the common Feynman-vertex factor i. L=sum g_j O_j/M_V.
Equal on-shell pion masses and a single-trace isoscalar source are assumed.
Run: .venv/bin/python scripts/derive_vector_three_pion_contacts.py
"""
import itertools
import sympy as sp
from sympy.combinatorics import Permutation

pion_decay_constant, vector_mass, resonance_mass = sp.symbols('f_pi M_i M_V', positive=True)
non_strange_entry, strange_entry = sp.symbols('t_n t_s', real=True)
pion_mass_squared, total_invariant = sp.symbols('m_pi_squared s', real=True)
pair_12, pair_13, pair_23 = sp.symbols('s_12 s_13 s_23', real=True)
couplings = sp.symbols('g_1:6', real=True)
plus = sp.sqrt(2)*sp.Matrix([[0,1,0],[0,0,0],[0,0,0]])
pion_generators = (plus, plus.T, sp.diag(1,-1,0))
isoscalar_generator = sp.diag(non_strange_entry,non_strange_entry,strange_entry)
pair_invariants = {(0,1):pair_12,(0,2):pair_13,(1,2):pair_23}


def pion_dot(first, second):
    """On-shell Minkowski scalar products from pair invariants."""
    if first == second:
        return pion_mass_squared
    return (pair_invariants[tuple(sorted((first,second)))]-2*pion_mass_squared)/2


def epsilon_projection(first, second):
    """Reduce epsilon(q,eps,p_first,p_second) to the charged basis.

    SymPy has no matrix-valued exterior coefficient extractor; this finite
    permutation sign is the multilinear epsilon expansion with q=sum p.
    """
    missing = next(index for index in range(3) if index not in (first,second))
    return Permutation([missing,3,first,second]).signature()


def derive_contact_coefficients(source):
    """Differentiate the cubic pion monomials by all Wick assignments."""
    coefficients = [sp.S.Zero for unused in range(5)]
    for first,second,third in itertools.permutations(range(3)):
        first_matrix,second_matrix,third_matrix = (pion_generators[index] for index in (first,second,third))
        # h(first) u(second) u(third): h=-2 ddPhi/f, u=-dPhi/f.
        coefficients[0] += -4*sp.I*pion_dot(first,second)*epsilon_projection(first,third)*sp.trace(source*(first_matrix*second_matrix*third_matrix-third_matrix*second_matrix*first_matrix))
        coefficients[1] += -4*sp.I*pion_dot(first,third)*epsilon_projection(first,second)*sp.trace(source*(first_matrix*second_matrix*third_matrix-third_matrix*second_matrix*first_matrix))
        # u(first) h(second) u(third).
        coefficients[2] += -4*sp.I*pion_dot(first,second)*epsilon_projection(second,third)*sp.trace(source*(first_matrix*second_matrix*third_matrix-third_matrix*second_matrix*first_matrix))
        # chi_minus=-2 i m_pi^2 Phi/f in the pion subspace.
        coefficients[3] += 4*sp.I*pion_mass_squared*epsilon_projection(first,second)*sp.trace((source*first_matrix*second_matrix+first_matrix*second_matrix*source)*third_matrix)
        # O5=epsilon tr[V u_beta chi_minus u_alpha].
        coefficients[4] += 4*sp.I*pion_mass_squared*epsilon_projection(third,first)*sp.trace(source*first_matrix*second_matrix*third_matrix)
    mass_shell = {pair_23:total_invariant+3*pion_mass_squared-pair_12-pair_13}
    return tuple(sp.simplify(sp.expand(value).subs(mass_shell)/(vector_mass*pion_decay_constant**3)) for value in coefficients)

operator_coefficients = derive_contact_coefficients(isoscalar_generator)
lagrangian_coefficient = sp.factor(sum(coupling*coefficient for coupling,coefficient in zip(couplings,operator_coefficients))/resonance_mass)
# A scalar independent of the Dalitz variables multiplying the geometric
# current is already pure total J=1; averaging n_i n_j gives delta_ij/3.
orientation_cosine, orientation_azimuth = sp.symbols("z phi", real=True)
unit_normal=sp.Matrix([sp.sqrt(1-orientation_cosine**2)*sp.cos(orientation_azimuth),sp.sqrt(1-orientation_cosine**2)*sp.sin(orientation_azimuth),orientation_cosine])
orientation_second_moment=(unit_normal*unit_normal.T).applyfunc(lambda value: sp.integrate(value,(orientation_azimuth,0,2*sp.pi),(orientation_cosine,-1,1))/(4*sp.pi))


def show(label,expression):
    """Print the deterministic derivation in LaTeX."""
    print(f'{label}:\n{sp.latex(expression)}\n')


if __name__ == '__main__':
    show('On-shell pair identity',sp.Eq(pair_12+pair_13+pair_23,total_invariant+3*pion_mass_squared))
    show('Coefficients of O1 through O5 before the common vertex i',sp.Tuple(*operator_coefficients))
    show('Coefficient in L=sum g_j O_j/M_V',lagrangian_coefficient)
    show('Total-orientation unit-normal second moment',orientation_second_moment)


def derive_tensor_source_coefficients(source):
    """Coefficients of epsilon_{mu nu alpha beta} p_a^alpha p_b^beta.

    The result contracts V^{mu nu} directly, with no vector polarization
    or propagator assumption. Keep all three independent pion two-forms.
    """
    coefficients = [[sp.S.Zero]*3 for unused in range(5)]
    pairs = ((0,1),(0,2),(1,2))
    def add(operator,first,second,value):
        pair = tuple(sorted((first,second)))
        coefficients[operator][pairs.index(pair)] += value*(1 if first<second else -1)
    for first,second,third in itertools.permutations(range(3)):
        first_matrix,second_matrix,third_matrix = (pion_generators[index] for index in (first,second,third))
        flavor = sp.trace(source*(first_matrix*second_matrix*third_matrix-third_matrix*second_matrix*first_matrix))
        add(0,first,third,-2*sp.I*pion_dot(first,second)*flavor)
        add(1,first,second,-2*sp.I*pion_dot(first,third)*flavor)
        add(2,second,third,-2*sp.I*pion_dot(first,second)*flavor)
        add(3,first,second,2*sp.I*pion_mass_squared*sp.trace((source*first_matrix*second_matrix+first_matrix*second_matrix*source)*third_matrix))
        add(4,third,first,2*sp.I*pion_mass_squared*sp.trace(source*first_matrix*second_matrix*third_matrix))
    return sp.Matrix([[sp.factor(value/pion_decay_constant**3) for value in row] for row in coefficients])


tensor_source_coefficients = derive_tensor_source_coefficients(isoscalar_generator)
if __name__ == '__main__':
    show('Unprojected tensor coefficients; rows O1..5, columns dual(p+ p-), dual(p+ p0), dual(p- p0)',tensor_source_coefficients)

# Free-tensor response for a completely general three-pion source in its
# rest frame. This follows the inverse Hessian established in the linked
# pole-factorization record; no restriction to a wedge source is made.
source_01,source_02,source_03,source_12,source_13,source_23=sp.symbols('B_01 B_02 B_03 B_12 B_13 B_23')
tensor_source=sp.Matrix([[0,source_01,source_02,source_03],[-source_01,0,source_12,source_13],[-source_02,-source_12,0,source_23],[-source_03,-source_13,-source_23,0]])
rest_momentum=sp.Matrix([sp.sqrt(total_invariant),0,0,0])
metric=sp.diag(1,-1,-1,-1)
source_divergence=tensor_source.T*metric*rest_momentum
spin_one_source=sp.simplify((rest_momentum*source_divergence.T-source_divergence*rest_momentum.T)/total_invariant)
local_source=sp.simplify(tensor_source-spin_one_source)
tensor_response=2*spin_one_source/(total_invariant-vector_mass**2)-2*local_source/vector_mass**2
if __name__ == '__main__':
    show('Spin-one tensor-source projection in the three-pion rest frame',spin_one_source)
    show('Local tensor-source projection in the three-pion rest frame',local_source)
    show('Full free antisymmetric-tensor response',tensor_response)
