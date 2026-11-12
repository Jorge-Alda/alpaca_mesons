#!/usr/bin/env python3
"""Pole-aligned symmetric chiral map and its induced LO vertices.

Run: .venv/bin/python scripts/derive_alp_pole_aligned_vertices.py
Natural units; all momenta incoming; CCWZ vertices omit the overall i.
The input meson generators form a complete orthonormal U(N) basis.
"""
from itertools import combinations

import sympy as sp

from derive_alp_lo_joint_vertices import quadratic_blocks
from derive_alp_pole_projected_amplitudes import ccwz_vertex, projected_amplitude
from solve_alp_meson_spectrum import PoleSpectrum, solve_spectrum_exact


def alignment_matrix(pole_vector):
    """Return the original-to-primed linear map for a regular parent pole."""
    pole_vector = sp.Matrix(pole_vector)
    if pole_vector[-1] == 0:
        raise ValueError('the selected pole has no ALP component')
    result = sp.eye(len(pole_vector))
    for index in range(len(pole_vector)-1):
        result[index, -1] = pole_vector[index] / pole_vector[-1]
    return result


def cubic_map(basis, alignment, directions, f_pi):
    """Third derivative of log(exp(A) exp(B) exp(A)) in primed directions.

    A=i a' w_i lambda_i/(2 f_pi), B=i P'_i lambda_i/f_pi.
    The matrix logarithm is expanded around the identity, on the vacuum branch.
    """
    variables = sp.symbols('t0:3', real=True)
    zero = sp.zeros(basis[0].rows)
    meson = sum((variables[k] * sum((directions[k][j] * basis[j]
                 for j in range(len(basis))), zero) for k in range(3)), zero)
    alp = sum(variables[k] * directions[k][-1] for k in range(3))
    shift = sum((alignment[j, -1] * basis[j] for j in range(len(basis))), zero)
    a_matrix = sp.I * alp * shift / (2*f_pi)
    b_matrix = sp.I * meson / f_pi
    identity = sp.eye(zero.rows)
    series_a = [identity, a_matrix, a_matrix**2/2, a_matrix**3/6]
    series_b = [identity, b_matrix, b_matrix**2/2, b_matrix**3/6]
    degree = [sum((series_a[i]*series_b[j]*series_a[n-i-j]
                   for i in range(n+1) for j in range(n-i+1)), zero)
              for n in range(4)]
    logarithm = degree[3] - (degree[1]*degree[2]+degree[2]*degree[1])/2 + degree[1]**3/3
    coordinate = -sp.I*f_pi*logarithm
    result = [sp.expand(sp.diff(sp.trace(generator*coordinate)/2, *variables))
              for generator in basis]
    return sp.Matrix(result + [0])


def transformed_vertex_parts(basis, sources, alignment, directions, momenta, gram):
    """Linear-map vertex and cubic-map contact of the exact axial map.

    This is complete through four external legs about a stationary vacuum.
    The cubic map has no quadratic term and a possible quartic coordinate
    term multiplies a vanishing one-point vertex.
    """
    directions = [sp.Matrix(direction) for direction in directions]
    count = len(directions)
    if count not in (2, 3, 4):
        raise ValueError('two to four legs supported')
    linear = ccwz_vertex(basis, sources, [alignment*d for d in directions], momenta, gram)
    nonlinear = sp.S.Zero
    if count == 4:
        for indices in combinations(range(4), 3):
            remaining = next(i for i in range(4) if i not in indices)
            mapped = cubic_map(basis, alignment, [directions[i] for i in indices], sources[4])
            total_momentum = sum((momenta[i] for i in indices), sp.zeros(len(momenta[0]), 1))
            nonlinear += ccwz_vertex(basis, sources,
                                      [mapped, alignment*directions[remaining]],
                                      [total_momentum, momenta[remaining]], gram)
    return sp.expand(linear), sp.expand(nonlinear)


def transformed_vertex(basis, sources, alignment, directions, momenta, gram):
    """Complete amputated vertex through four legs in aligned coordinates."""
    return sum(transformed_vertex_parts(basis, sources, alignment, directions, momenta, gram))


def nonlinear_contact_parts(basis, sources, alignment, directions, momenta, gram):
    """Return the kinetic and mass contributions of the cubic-map contact."""
    _, total = transformed_vertex_parts(basis, sources, alignment, directions, momenta, gram)
    _, mass = transformed_vertex_parts(basis, sources, alignment, directions, momenta, sp.zeros(*gram.shape))
    return sp.expand(total-mass), mass


def aligned_amplitude(basis, sources, spectrum, pole_indices, gram):
    """Return the on-shell amplitude in the pole-aligned coordinates.

    The first external leg selects the regular parent; every external pole
    vector and the complete internal propagator are transformed consistently.
    """
    change = alignment_matrix(spectrum.pole_vector(pole_indices[0]))
    inverse = change.inv()
    vectors = inverse * sp.Matrix(spectrum.pole_vectors)
    transformed = PoleSpectrum(spectrum.field_labels, change.T*sp.Matrix(spectrum.kinetic)*change,
                               change.T*sp.Matrix(spectrum.mass_squared)*change,
                               spectrum.pole_squared_masses, vectors,
                               [vectors[:, i]*vectors[:, i].T for i in range(spectrum.size)], True)
    def vertex(directions, momenta, products):
        return transformed_vertex(basis, sources, change, directions, momenta, products)
    return projected_amplitude(transformed, pole_indices, gram, vertex)


def main():
    pauli = [sp.Matrix([[0,1],[1,0]]), sp.Matrix([[0,-sp.I],[sp.I,0]]), sp.diag(1,-1), sp.eye(2)]
    zero = sp.zeros(2)
    sources = ([sp.eye(2), sp.eye(2), pauli[2], pauli[0], zero],
               [pauli[2]/2, pauli[0]/3, zero, zero, zero],
               [-pauli[2]/2, pauli[0]/3, zero, zero, zero], [sp.S.Zero]*5,
               sp.S.One, sp.Integer(2), sp.S.Zero, sp.Rational(236, 15))
    kinetic, mass = quadratic_blocks(pauli, *sources)
    spectrum = solve_spectrum_exact(kinetic, mass)
    indices = [4, 1, 1, 0]
    vectors = [sp.Matrix(spectrum.pole_vector(i)) for i in indices]
    alignment = alignment_matrix(vectors[0])
    directions = [alignment.inv()*v for v in vectors]
    momenta = [sp.eye(4)[:,i] for i in range(4)]
    gram = sp.diag(16, 1, 1, sp.Rational(59,60))
    gram[0,1] = gram[1,0] = -2
    gram[0,2] = gram[2,0] = -3
    # The final row and column are fixed by momentum conservation.
    gram[1,2] = gram[2,1] = (sp.Rational(59,60)-18)/2 + 5
    for i in range(3):
        gram[i,3] = gram[3,i] = -sum(gram[i,j] for j in range(3))
    contact = ccwz_vertex(pauli, sources, vectors, momenta, gram)
    amplitude = aligned_amplitude(pauli, sources, spectrum, indices, gram)
    print('Illustrative exact on-shell nonlinear kinetic and mass contacts:',
          nonlinear_contact_parts(pauli, sources, alignment, directions, momenta, gram))
    print('Linear-map contact:', contact, 'exchange sum:', sp.simplify(amplitude-contact))
    print('Full aligned amplitude:', amplitude)


if __name__ == '__main__':
    main()
