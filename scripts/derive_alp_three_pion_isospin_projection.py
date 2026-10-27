"""Derive the I=0 charged component of the LO Born current.

Invocation: .venv/bin/python scripts/derive_alp_three_pion_isospin_projection.py
Inputs: memory/alp/lepton-pair/alp-three-pion-current-born-amplitude.md.
Natural units; metric (+---).
The isospin projection assumes an exact pion triplet, common pion masses,
T_0=lambda_3, and fixed external source coefficients under permutations.
It is an isospin projection, not a J=1 angular projection.
"""
from itertools import permutations
import sympy as sp
from derive_alp_virtual_photon_three_pion_born import derive as derive_born


def derive():
    """Apply the normalized three-pion isospin antisymmetrizer."""
    born = derive_born()
    vertex = born['isospin_source_amplitude']
    axial_source = sp.MatrixSymbol('c_q^A', 3, 3)
    pseudoscalar_mass_source = sp.MatrixSymbol('X_P', 3, 3)
    positive_generator = sp.zeros(3)
    positive_generator[0, 1] = sp.sqrt(2)
    negative_generator = positive_generator.T
    neutral_generator = sp.diag(1, -1, 0)
    generators = (positive_generator, negative_generator, neutral_generator)
    pion_fields = sp.symbols('pi_plus pi_minus pi_zero')
    pion_derivatives = sp.symbols('dpi_plus dpi_minus dpi_zero')
    phi = sum((field * generator for field, generator in zip(pion_fields, generators)), sp.zeros(3))
    derivative_phi = sum((field * generator for field, generator in zip(pion_derivatives, generators)), sp.zeros(3))
    commutator = phi * derivative_phi - derivative_phi * phi
    double_commutator = phi * commutator - commutator * phi
    kinetic_matrix = double_commutator.applyfunc(lambda entry: sp.expand(entry).coeff(pion_derivatives[0] * pion_fields[1] * pion_fields[2]))
    mass_matrix = (phi**3).applyfunc(lambda entry: sp.expand(entry).coeff(pion_fields[0] * pion_fields[1] * pion_fields[2]))
    kinetic_trace = sp.Trace(axial_source * kinetic_matrix)
    mass_trace = sp.Trace(pseudoscalar_mass_source * mass_matrix)
    trace_vertex = (kinetic_trace * (sp.Symbol('m_a2') - 3 * sp.Symbol('p_dot_p0')) - sp.Symbol('B_0', nonzero=True) * mass_trace) / (6 * sp.Symbol('f_a', nonzero=True) * sp.Symbol('f_pi', nonzero=True))
    source_components = {axial_source: sp.diag(*sp.symbols('c_Au c_Ad c_As')), pseudoscalar_mass_source: sp.diag(*sp.symbols('P_u P_d P_s'))}
    trace_vertex_difference = sp.factor(trace_vertex.subs(source_components).doit() - vertex)
    electric_charge = sp.Symbol('e')
    photon_squared = sp.Symbol('q2')
    alp_products = sp.symbols('p_dot_pplus p_dot_pminus p_dot_p0')
    photon_products = sp.symbols('q_dot_pplus q_dot_pminus q_dot_p0')
    denominators = tuple(photon_squared - 2 * product for product in photon_products)
    # Independent formal vector components suffice because all operations
    # are linear. No custom tensor algebra or angular projection is used.
    pion_components = sp.symbols('pplus_mu pminus_mu pzero_mu')
    photon_component = sp.Symbol('q_mu')
    photon_vertices = tuple(2 * component - photon_component for component in pion_components)
    neutral_product = sp.Symbol('p_dot_p0')
    def ordered_current(order):
        positive, negative, neutral = order
        return electric_charge * vertex.subs(neutral_product, alp_products[neutral]) * (
            photon_vertices[positive] / denominators[positive]
            - photon_vertices[negative] / denominators[negative])
    projected = 0
    trace_projected = 0
    for order in permutations(range(3)):
        parity = (-1)**sum(order[first] > order[second] for first in range(3) for second in range(first+1, 3))
        term = ordered_current(order)
        projected += parity * term / 6
        positive, negative, neutral = order
        trace_projected += parity * electric_charge * trace_vertex.subs(neutral_product, alp_products[neutral]) * (photon_vertices[positive] / denominators[positive] - photon_vertices[negative] / denominators[negative]) / 6
    ward = sp.factor(projected.subs(dict(zip(pion_components + (photon_component,), photon_products + (photon_squared,)))))
    coefficients = tuple(sp.factor(sp.diff(projected, component)) for component in pion_components)
    # The constant part cancels before any values of the source matrices
    # are selected. Differentiate with respect to the axial/pseudoscalar inputs.
    mass_source_coefficients = tuple(sp.factor(sp.diff(projected, source)) for source in sp.symbols('P_u P_d'))
    compact_projected = sp.Add(*(sp.factor(sp.diff(projected, component)) * component for component in pion_components + (photon_component,)), evaluate=False)
    trace_projected = sp.collect(sp.expand(trace_projected), (kinetic_trace, mass_trace))
    trace_projection_difference = sp.factor(trace_projected.subs(source_components).doit() - projected)
    return dict(kinetic_trace=kinetic_trace, mass_trace=mass_trace, trace_vertex=trace_vertex, trace_projected=trace_projected, trace_vertex_difference=trace_vertex_difference, trace_projection_difference=trace_projection_difference, vertex=vertex, projected=compact_projected, coefficients=coefficients,
                ward=ward, mass_source_coefficients=mass_source_coefficients,
                photon_component_coefficient=sp.factor(sp.diff(projected, photon_component)))


def main():
    """Print the executable derivation in deterministic LaTeX form."""
    print('Exact-isospin pion triplet; natural units; e retained. No J=1 projection or rescattering.')
    results = derive()
    print('T_0=diag(1,-1,0); X_P={c_q^P,m_q}; source matrices diagonal in the no-mixing basis.')
    for label in ('kinetic_trace', 'mass_trace', 'trace_vertex', 'trace_projected', 'trace_vertex_difference', 'trace_projection_difference', 'vertex', 'projected', 'coefficients', 'mass_source_coefficients', 'photon_component_coefficient', 'ward'):
        print(f'{label}:\n{sp.latex(results[label])}\n')
    print('P0 J is the charged component: (1/6) sum of six signed permutations; even +, odd -.')
    print('This is not the normalized I=0 state amplitude. No J=1 angular projection is performed.')


if __name__ == '__main__':
    main()
