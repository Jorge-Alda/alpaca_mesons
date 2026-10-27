"""Derive the five VPPP charged-pion projections using SymPy.

Invocation: .venv/bin/python scripts/derive_vppp_isoscalar_three_pion.py
Natural units. Exact isospin, on-shell pions, cubic pion order, tree level.
Inputs: memory/extracts/GomezDumm:2010ra/vppp-basis-and-omega-three-pion.md;
memory/extracts/Alda:2026nrd/appendix-a-expansion.md;
memory/rcht/odd-parity/levi-civita-convention.md.
The displayed scalar divides out i epsilon(p+,p-,p0,epsilon) and uses
V^{mu nu}=(q^mu epsilon^nu-q^nu epsilon^mu)/M_i. All pion derivatives
carry +i p. The overall vector LSZ phase is common and is not fitted.
Explicit flavor matrices are required for the requested charged channel.
"""

import itertools
import unittest
import sympy as sp

pion_mass_squared, total_invariant = sp.symbols('m_pi_squared s', real=True)
pair_12, pair_13, pair_23 = sp.symbols('s_12 s_13 s_23', real=True)
non_strange_entry, strange_entry = sp.symbols('t_n t_s', real=True)
pion_decay_constant, vector_mass, basis_mass = sp.symbols('f_pi M_i M_V', positive=True)
couplings = sp.symbols('g_1:6', real=True)
pion_generators = (sp.sqrt(2)*sp.Matrix([[0,1,0],[0,0,0],[0,0,0]]),
                   sp.sqrt(2)*sp.Matrix([[0,0,0],[1,0,0],[0,0,0]]),
                   sp.diag(1,-1,0))
isoscalar_generator = sp.diag(non_strange_entry, non_strange_entry, strange_entry)
momentum_basis = tuple(sp.eye(4)[:, index] for index in range(4))
total_momentum = sum(momentum_basis[:3], sp.zeros(4,1))
polarization = momentum_basis[3]
pion_scalar_products = sp.Matrix([
    [pion_mass_squared,(pair_12-2*pion_mass_squared)/2,(pair_13-2*pion_mass_squared)/2],
    [(pair_12-2*pion_mass_squared)/2,pion_mass_squared,(pair_23-2*pion_mass_squared)/2],
    [(pair_13-2*pion_mass_squared)/2,(pair_23-2*pion_mass_squared)/2,pion_mass_squared]])


def epsilon_ratio(first_index, second_index):
    """Evaluate an epsilon ratio by SymPy determinant, independent of orientation."""
    return sp.det(sp.Matrix.hstack(total_momentum, polarization,
                                  momentum_basis[first_index], momentum_basis[second_index]))


def derive_operator_projections(source=isoscalar_generator):
    """Expand all six pion assignments of the five literal operator words.

    The two h derivatives give -2 i^2/f_pi; each u gives -i/f_pi;
    chi_minus gives -2 i m_pi^2/f_pi. Division by the common i is explicit.
    SymPy has no matrix-valued field derivative with momentum slots, so this
    finite Wick-assignment loop keeps those slots while using native traces.
    """
    result = sp.zeros(1,5)
    for first_index, second_index, third_index in itertools.permutations(range(3)):
        first, second, third = (pion_generators[index] for index in
                                (first_index,second_index,third_index))
        trace_123 = sp.trace(source*first*second*third)
        trace_321 = sp.trace(source*third*second*first)
        # Operators 1 and 2: h is in the first pion slot, beta in the third.
        result[0,0] += (-4/pion_decay_constant**3/vector_mass
                       *epsilon_ratio(first_index,third_index)
                       *pion_scalar_products[first_index,second_index]
                       *(trace_123-trace_321))
        result[0,1] += (-4/pion_decay_constant**3/vector_mass
                       *epsilon_ratio(first_index,second_index)
                       *pion_scalar_products[first_index,third_index]
                       *(trace_123-trace_321))
        # Operator 3: h is in the second pion slot, beta in the third.
        result[0,2] += (-4/pion_decay_constant**3/vector_mass
                       *epsilon_ratio(second_index,third_index)
                       *pion_scalar_products[first_index,second_index]
                       *(trace_123-trace_321))
        result[0,3] += (4*pion_mass_squared/pion_decay_constant**3/vector_mass
                       *epsilon_ratio(first_index,second_index)
                       *sp.trace((source*first*second+first*second*source)*third))
        result[0,4] += (4*pion_mass_squared/pion_decay_constant**3/vector_mass
                       *epsilon_ratio(first_index,second_index)
                       *sp.trace(first*source*second*third))
    return result.applyfunc(sp.simplify)


def derive_results():
    """Reduce momentum conservation and extract independent coefficient rows."""
    raw_projection = derive_operator_projections()
    projection = raw_projection.subs(pair_23,total_invariant+3*pion_mass_squared-pair_12-pair_13).applyfunc(sp.simplify)
    amplitude = sp.expand((projection*sp.Matrix(couplings))[0]/basis_mass)
    coefficient_matrix = sp.Matrix([[sp.expand(projection[0,index]).coeff(total_invariant) for index in range(5)],
                                   [sp.expand(projection[0,index]).coeff(pion_mass_squared) for index in range(5)]])
    return raw_projection, projection, amplitude, coefficient_matrix


class ProjectionTests(unittest.TestCase):
    """Independent physical degeneracy and published on-shell checks."""
    def test_published_omega_polynomial(self):
        _,_,amplitude,_ = derive_results()
        expected = 8*sp.sqrt(2)/(vector_mass*pion_decay_constant**3*basis_mass)*((couplings[0]-couplings[1]-couplings[2])*(total_invariant-3*pion_mass_squared)+3*pion_mass_squared*(2*couplings[3]+couplings[4]))
        self.assertEqual(sp.simplify(amplitude.subs(non_strange_entry,1/sp.sqrt(2))-expected),0)

    def test_strange_and_isovector_decouple(self):
        self.assertEqual(derive_operator_projections(sp.diag(0,0,1)),sp.zeros(1,5))
        self.assertEqual(derive_operator_projections(sp.diag(1,-1,0)),sp.zeros(1,5))

    def test_rank_and_nullspace(self):
        _,projection,_,matrix = derive_results()
        self.assertEqual(matrix.rank(),2)
        self.assertEqual(projection.rank(),1)
        chiral_projection = projection.subs(pion_mass_squared,0)
        self.assertEqual(chiral_projection.rank(),1)
        for direction in (sp.Matrix([1,1,0,0,0]),sp.Matrix([1,0,1,0,0]),sp.Matrix([0,0,0,1,-2])):
            self.assertEqual(matrix*direction,sp.zeros(2,1))

    def test_two_energies_and_p_wave(self):
        _,projection,amplitude,_ = derive_results()
        first_energy, second_energy = sp.symbols('s_a s_b', real=True)
        energy_rows = projection.subs(total_invariant,first_energy).col_join(
            projection.subs(total_invariant,second_energy))
        self.assertEqual(energy_rows.rank(),2)
        self.assertEqual(energy_rows.subs(second_energy,first_energy).rank(),1)
        self.assertEqual(energy_rows.subs(pion_mass_squared,0).rank(),1)
        angle_cosine = sp.Symbol('z', real=True)
        self.assertEqual(sp.simplify(sp.integrate(
            sp.Rational(3,4)*(1-angle_cosine**2)*amplitude,
            (angle_cosine,-1,1))-amplitude),0)


def main():
    """Print each analytical step using the standard SymPy LaTeX renderer."""
    print('Exact isospin; pion external legs on shell; tree level; natural units.\n')
    raw,projection,amplitude,matrix = derive_results()
    angle_cosine = sp.Symbol('z', real=True)
    p_wave = sp.expand(sp.integrate(
        sp.Rational(3,4)*(1-angle_cosine**2)*amplitude,
        (angle_cosine,-1,1)))
    first_energy, second_energy = sp.symbols('s_a s_b', real=True)
    two_energy_projection = projection.subs(total_invariant,first_energy).col_join(
        projection.subs(total_invariant,second_energy))
    for label,expression in [('Pion scalar-product matrix',pion_scalar_products),
                             ('Five projections before total-invariant substitution',raw),
                             ('Five projections after momentum conservation',projection),
                             ('Local scalar including g_i/M_V',amplitude),
                             ('Local P-wave projection (3/4) integral (1-z^2) F dz',p_wave),
                             ('Rows of coefficients of s and m_pi_squared',matrix),
                             ('Coefficient rank',sp.Integer(matrix.rank())),
                             ('Polynomial coefficient nullspace',sp.Matrix.hstack(*matrix.nullspace())),
                             ('Fixed-energy coefficient rank',sp.Integer(projection.rank())),
                             ('Fixed-energy coefficient nullspace',sp.Matrix.hstack(*projection.nullspace())),
                             ('Two-distinct-energy coefficient matrix',two_energy_projection),
                             ('Two-distinct-energy coefficient rank (m_pi_squared nonzero)',sp.Integer(two_energy_projection.rank())),
                             ('Two-distinct-energy coefficient nullspace',sp.Matrix.hstack(*two_energy_projection.nullspace())),
                             ('Two-distinct-energy rank in the chiral limit',sp.Integer(two_energy_projection.subs(pion_mass_squared,0).rank()))]:
        print(f'{label}:\n{sp.latex(expression)}\n')


if __name__ == '__main__':
    main()
