"""Enumerate a -> 4P vector trees and derive their odd vertex polynomial.

Run: .venv/bin/python scripts/derive_alp_four_pseudoscalar.py
Natural units (GeV), metric (+---), epsilon^{0123}=+1. Pseudoscalar
normalization is Phi=sum lambda_i P_i, Tr(lambda_i lambda_j)=2 delta_ij;
vector generators have unit trace normalization. The ALP sources are already
in the no-mixing basis. This calculation retains the original preliminary-fit
Lagrangian at first order in ALP sources and zero electromagnetic coupling.
It does not add independent VPPP operators. The generic tensor contractions
printed below specify topology and combinatorics; they are not numerical
channel amplitudes.
"""
from dataclasses import dataclass
from itertools import combinations, combinations_with_replacement, product

import sympy as sp

from derive_alp_vector_two_pseudoscalar import derive as derive_vvp
from derive_radiative_tff_amplitude import m, s1, s2


@dataclass(frozen=True)
class Vertex:
    """Field multiplicities after selecting exactly one ALP insertion."""

    name: str
    vectors: int
    pseudoscalars: int


# Both odd and even aVV vertices share field multiplicities. Their distinct
# flavour and Lorentz tensors must subsequently be retained separately.
SOURCE_VERTICES = (
    Vertex('aPP', 0, 2), Vertex('aPPP', 0, 3), Vertex('aPPPP', 0, 4),
    Vertex('aVP', 1, 1), Vertex('aVPP', 1, 2), Vertex('aVPPP', 1, 3),
    Vertex('aVV', 2, 0), Vertex('aVVP', 2, 1), Vertex('aVVPP', 2, 2),
)
STRONG_VERTICES = (
    Vertex('PPPP', 0, 4), Vertex('VPP', 1, 2),
    Vertex('VPPPP', 1, 4), Vertex('VVP', 2, 1),
    Vertex('VVPP', 2, 2), Vertex('VVPPP', 2, 3),
)


def enumerate_vector_trees():
    """Obtain every tree with one source, four external P and internal V.

    Every interaction has at least three total fields, including the ALP.
    For a tree with five external fields this bounds the number of vertices
    by three. A source aP two-point insertion is absent in the no-mixing
    basis. Quadratic masses are included in propagators. The small explicit
    graph enumeration is used because SymPy has no typed Feynman-tree
    enumerator; all analytical vertex polynomials below use SymPy.
    """
    results = {}
    for source in SOURCE_VERTICES:
        for strong_count in (1, 2):
            for strong in combinations_with_replacement(STRONG_VERTICES, strong_count):
                vertices = (source,) + strong
                vertex_count = len(vertices)
                possible_edges = tuple(combinations(range(vertex_count), 2))
                for edges in combinations(possible_edges, vertex_count - 1):
                    reached = {0}
                    for _ in range(vertex_count):
                        reached.update(right for left, right in edges if left in reached)
                        reached.update(left for left, right in edges if right in reached)
                    if len(reached) != vertex_count:
                        continue
                    for species in product(('P', 'V'), repeat=len(edges)):
                        if 'V' not in species:
                            continue
                        remaining_v = [vertex.vectors for vertex in vertices]
                        remaining_p = [vertex.pseudoscalars for vertex in vertices]
                        for (left, right), particle in zip(edges, species):
                            remaining = remaining_v if particle == 'V' else remaining_p
                            remaining[left] -= 1
                            remaining[right] -= 1
                        if min(remaining_v + remaining_p) < 0:
                            continue
                        if any(remaining_v) or sum(remaining_p) != 4:
                            continue
                        key = tuple(vertex.name for vertex in vertices)
                        # Equivalent internal labels are counted only once here.
                        # External assignments are summed separately, without a
                        # symmetry factor in the amplitude.
                        results.setdefault(key, (vertices, edges, species, tuple(remaining_p)))
    return tuple(results[key] for key in sorted(results))


def derive_odd_numerator():
    """Derive the two off-shell vector numerator from the operator contraction."""
    calculation = derive_vvp()
    source_mass_squared, first_virtuality, second_virtuality = sp.symbols(
        'm_X_squared s_R s_S', real=True)
    numerator = sp.expand(calculation['strong'].subs({
        m**2: source_mass_squared, s1: first_virtuality, s2: second_virtuality}))
    coefficient_row = sp.Matrix([[sp.diff(numerator, coefficient)
                                 for coefficient in sp.symbols('c_1:8 d_1:5', real=True)]])
    return numerator, sp.simplify(coefficient_row * calculation['nulls'])


def derive_tensor_response():
    """Derive the propagator on a PP current using the existing tensor Hessian."""
    from derive_alp_vector_two_pseudoscalar import derive_tensor_exchange
    hessian, response, expected = derive_tensor_exchange()
    return hessian, response, sp.simplify(response - expected)


def derive_pairing_signs():
    """Reduce each two-pair epsilon contraction to epsilon(p1,p2,p3,p4)."""
    # A symbolic 4x4 matrix supplies four independent external momenta;
    # determinant multilinearity gives the coefficient without imposing
    # masses or kinematic restrictions.
    momenta = sp.eye(4)
    results = []
    for second in (1, 2, 3):
        first_pair = (0, second)
        second_pair = tuple(index for index in range(4) if index not in first_pair)
        ordered = first_pair + second_pair
        coefficient = sp.det(sp.Matrix.vstack(*(momenta[index, :] for index in ordered)))
        results.append((tuple(index + 1 for index in ordered), coefficient))
    return results


def derive_scalar_mass_vertices():
    """Extract the fitted vector mass terms after removing the vector source.

    The change of flavour variables moves a pure-gradient c_q^VS source
    into delta chi_plus=-2 i a [c_q^VS,chi]/f_a. This is a representation
    change at first order in a, not an additional phenomenological operator.
    The invariant M_V and Delta_V masses have zero variation. Both retained
    quark-mass terms, including kappa_V, are varied explicitly.
    """
    source, mass, first, second = [sp.MatrixSymbol(name, 3, 3) for name in
                                 ('c_q^VS', 'chi', 'T_R', 'T_S')]
    alp_field, first_field, second_field = sp.symbols('a V_R V_S')
    alp_scale, mass_coefficient, double_trace_coefficient = sp.symbols(
        'f_a e_m^V kappa_V', real=True)
    vector = first_field * first + second_field * second
    delta_chi_plus = -2 * sp.I * alp_field * (source * mass - mass * source) / alp_scale
    mass_lagrangian = (-mass_coefficient * sp.Trace(vector * vector * delta_chi_plus)
                       + double_trace_coefficient * sp.Trace(vector)
                       * sp.Trace(vector * delta_chi_plus) / 4)
    # Polarization extracts a V_R V_S without differentiating matrix traces;
    # SymPy Trace's scalar derivative is unsupported for a product of traces.
    vertex = sum((-1)**(3 - sum(values)) * mass_lagrangian.subs(
        dict(zip((alp_field, first_field, second_field), values))).doit()
        for values in product((0, 1), repeat=3))
    return expand_traces(vertex)


def expand_traces(expression):
    """Apply SymPy's matrix expansion and trace linearity to display results."""
    return sp.expand(expression.replace(
        lambda term: isinstance(term, sp.Trace),
        lambda term: sp.Trace(sp.expand(term.arg)).doit()))


def show(label, expression):
    """Print deterministic symbolic expressions in the standard LaTeX format."""
    print(f'{label}:\n{sp.latex(expression)}\n')


def main():
    """Print the field-counting derivation and independently reduced vertices."""
    print(__doc__)
    for vertices, edges, species, external_counts in enumerate_vector_trees():
        print('Tree:', ' + '.join(vertex.name for vertex in vertices))
        print('Internal lines:', tuple((vertices[left].name, vertices[right].name, particle)
                                       for (left, right), particle in zip(edges, species)))
        print('External P at each vertex:', external_counts)
    numerator, projection = derive_odd_numerator()
    show('VVP numerator: D=Tr({T_R,T_S} U_X)/2; P=Tr({T_R,T_S} Y_X)/2', numerator)
    show('Original preliminary-fit coefficient null projection', projection)
    _, response, residual = derive_tensor_response()
    show('Full tensor propagator response on J times p1 wedge p2', response)
    show('Response minus 2 J (p1 wedge p2)/(q squared - M_R squared)', residual)
    for pairing, sign in derive_pairing_signs():
        show(f'Epsilon ordering {pairing} relative to (1,2,3,4)', sign)
    from derive_alp_vector_two_pseudoscalar import derive_general_vector_source_vertices
    show('Rotated scalar-source aPP Lagrangian coefficient',
         expand_traces(derive_general_vector_source_vertices()))
    show('Rotated scalar-source aVV Lagrangian coefficient, contracted Lorentz indices',
         derive_scalar_mass_vertices())
    print('aVV has separate odd VVP-source and even vector-source tensors. '
          'aVP is the axial G_V source; aVPP and even aVV arise from the '
          'vector source. aPP followed by VPP and VPP has one internal '
          'pseudoscalar and requires the LO chiral source vertex. '
          'The higher field terms in the same retained interactions cannot '
          'contribute additional vector trees at this external multiplicity. '
          'Pure pseudoscalar trees and electromagnetic exchanges are outside '
          'this vector-exchange subset. Each allowed flavour generator and '
          'each inequivalent assignment of four external labelled legs must '
          'be summed. Identical-particle phase-space factors do not multiply '
          'the amplitude.')


if __name__ == '__main__':
    import sys
    import unittest

    class DerivationTests(unittest.TestCase):
        """Focused checks of topology completeness and off-shell contractions."""

        def test_tree_inventory(self):
            expected = {('aPP', 'VPP', 'VPP'), ('aVP', 'VPP', 'VVP'),
                        ('aVPP', 'VPP'), ('aVV', 'VPP', 'VPP')}
            self.assertEqual({tuple(vertex.name for vertex in result[0])
                              for result in enumerate_vector_trees()}, expected)

        def test_pairing_parity(self):
            self.assertEqual([sign for _, sign in derive_pairing_signs()], [1, -1, 1])

        def test_off_shell_null_projection(self):
            self.assertEqual(derive_odd_numerator()[1], sp.zeros(1, 4))

        def test_full_tensor_response(self):
            self.assertEqual(derive_tensor_response()[2], sp.zeros(6, 1))

        def test_conserved_vector_source_mass_vertex(self):
            source = sp.MatrixSymbol('c_q^VS', 3, 3)
            self.assertEqual(sp.simplify(derive_scalar_mass_vertices().subs(
                source, sp.Identity(3)).doit()), 0)

    if '--test' in sys.argv:
        unittest.main(argv=[sys.argv[0]])
    else:
        main()
