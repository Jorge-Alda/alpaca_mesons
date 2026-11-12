#!/usr/bin/env python3
"""LO pole-projected amplitudes from CCWZ vertices, with exact mixed exchange.

Run: .venv/bin/python scripts/derive_alp_pole_projected_amplitudes.py
Tests: .venv/bin/python -m unittest discover -s scripts -p test_derive_alp_pole_projected_amplitudes.py
All momenta are incoming. Gram entries are p_i dot p_j in a common mass unit;
vertices omit the overall i. Four-leg reduced amplitude is contact minus exchange.
"""
import sympy as sp
from derive_alp_lo_joint_vertices import vertex
from solve_alp_meson_spectrum import inverse_propagator, solve_spectrum_exact


def schur_propagator(kinetic, mass_squared, squared_momentum):
    """Analytic block inverse, valid where the meson block and Schur complement invert."""
    inverse = inverse_propagator(kinetic, mass_squared, squared_momentum)
    meson = inverse[:-1, :-1]
    mixing = inverse[:-1, -1:]
    meson_inverse = meson.inv()
    denominator = (inverse[-1, -1] - (mixing.T * meson_inverse * mixing)[0])
    overlap = -meson_inverse * mixing
    unnormalized = overlap.col_join(sp.ones(1, 1))
    meson_part = meson_inverse.row_join(sp.zeros(meson.rows, 1)).col_join(sp.zeros(1, meson.rows + 1))
    return meson_part + unnormalized * unnormalized.T / denominator


def _scalar_products(expression, momenta, gram):
    """Replace q_i q_j placeholders by momentum scalar products."""
    symbols = sp.symbols(f'_q0:{len(momenta)}', real=True)
    polynomial = sp.Poly(sp.expand(expression), *symbols)
    output = 0
    for powers, coefficient in polynomial.terms():
        indices = [i for i, power in enumerate(powers) for _ in range(power)]
        if not indices:
            output += coefficient
        elif len(indices) == 2:
            output += coefficient * (momenta[indices[0]].T * gram * momenta[indices[1]])[0]
        else:
            raise ValueError('expected zero or two derivatives in a vertex')
    return sp.expand(output)


def ccwz_vertex(meson_basis, sources, directions, momenta, gram):
    """Compute one projected CCWZ vertex in the specified field directions.

    sources is (chi derivatives, right derivatives, left derivatives, theta
    derivatives, f_pi, f_a, anomaly mass squared, bare ALP mass squared).
    Each direction has meson components followed by the ALP component.
    """
    zero = sp.zeros(meson_basis[0].rows)
    generators = [sum((direction[i] * basis for i, basis in enumerate(meson_basis)), zero)
                  for direction in directions]
    momenta_symbols = sp.symbols(f'_q0:{len(directions)}', real=True)
    expression = vertex(generators, [direction[-1] for direction in directions],
                        momenta_symbols, *sources)
    return _scalar_products(expression, momenta, gram)


def projected_amplitude(spectrum, pole_indices, gram, vertex_function, pole_vectors=None):
    """Three- or four-external-leg tree amplitude using exact mixed propagation.

    vertex_function(directions, momenta, gram) returns the reduced amputated
    vertex. Directions are field-space columns, momenta are coefficients in the
    basis of all-incoming external momenta, and gram contains their products.
    Pole indices identify parent then daughters; gram must obey momentum
    conservation and place all external legs on their chosen poles. Optionally
    supply pole_vectors to select normalized linear combinations inside a
    degenerate subspace, including complex charged states. Directions use the
    all-incoming field convention: cross outgoing momenta and use the
    conjugate field direction appropriate to the outgoing charged particle.
    """
    count = len(pole_indices)
    if count not in (3, 4):
        raise ValueError('three or four external legs required')
    directions = ([sp.Matrix(spectrum.pole_vector(i)) for i in pole_indices]
                  if pole_vectors is None else [sp.Matrix(v) for v in pole_vectors])
    if len(directions) != count:
        raise ValueError('one pole vector per external leg required')
    momenta = [sp.eye(count)[:, i] for i in range(count)]
    contact = vertex_function(directions, momenta, gram)
    if count == 3:
        return sp.expand(contact)
    total = contact
    basis = [sp.eye(spectrum.size)[:, i] for i in range(spectrum.size)]
    for first, second in (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))):
        internal = -sum((momenta[i] for i in first), sp.zeros(count, 1))
        squared_momentum = (internal.T * gram * internal)[0]
        propagator = inverse_propagator(sp.Matrix(spectrum.kinetic), sp.Matrix(spectrum.mass_squared), squared_momentum).inv()
        left = sp.Matrix([vertex_function([directions[i] for i in first] + [direction],
                                          [momenta[i] for i in first] + [internal], gram) for direction in basis])
        right = sp.Matrix([vertex_function([directions[i] for i in second] + [direction],
                                           [momenta[i] for i in second] + [-internal], gram) for direction in basis])
        total -= (left.T * propagator * right)[0]
    return sp.expand(total)


def main():
    s, meson_mass_squared, alp_mass_squared, mass_mixing, kinetic_mixing = sp.symbols(
        's m_P_squared m_a_squared M_Pa_squared K_Pa', real=True)
    kinetic = sp.Matrix([[1, kinetic_mixing], [kinetic_mixing, 1 + kinetic_mixing**2]])
    mass_squared = sp.Matrix([[meson_mass_squared, mass_mixing], [mass_mixing, alp_mass_squared]])
    inverse = inverse_propagator(kinetic, mass_squared, s)
    overlap = -inverse[0, 1] / inverse[0, 0]
    denominator = sp.factor(inverse[1, 1] - inverse[1, 0]*inverse[0, 1]/inverse[0, 0])
    print('Two-field analytical overlap and Schur complement:')
    sp.pprint(overlap)
    sp.pprint(denominator)
    print('Schur inverse equals direct inverse:',
          sp.simplify(schur_propagator(kinetic, mass_squared, s) - inverse.inv()) == sp.zeros(2))
    spectrum = solve_spectrum_exact(sp.eye(2), sp.diag(1, 4))
    print('Example pole-sum equality:',
          sp.simplify(schur_propagator(sp.eye(2), sp.diag(1, 4), s) - spectrum.propagator_from_poles(s)) == sp.zeros(2))

if __name__ == '__main__':
    main()
