#!/usr/bin/env python3
"""Evaluate the LO CCWZ cubic amplitude for three arbitrary pole directions.

Run: .venv/bin/python scripts/derive_alp_explicit_neutral_amplitude.py
Momenta are incoming, q_i q_j represents p_i.p_j and no Feynman i is included.
Natural units; chi derivatives are derivatives in x=a/f_a, not coefficients.
"""
from itertools import permutations
import sympy as sp


def cubic_amplitude(meson_matrices, alp_components, momenta, chi_derivatives,
                    right_derivatives, left_derivatives, theta_derivatives,
                    f_pi, f_a, anomaly_mass_squared):
    """Evaluate all cubic CCWZ terms for the supplied normalized pole directions."""
    chi1, chi2, chi3 = chi_derivatives[1:4]
    right0, right1 = right_derivatives[:2]
    left0, left1 = left_derivatives[:2]
    dimension = sp.Integer(meson_matrices[0].rows)
    result = 0
    for i, j, k in permutations(range(3)):
        # Derivatives act on incoming waves, hence -i p.
        first = -meson_matrices[i]/f_pi + alp_components[i]*(right0-left0)/f_a
        commutator = (right0+left0)*meson_matrices[k]-meson_matrices[k]*(right0+left0)
        second = (sp.I*alp_components[j]*commutator/(2*f_pi*f_a)
                  + alp_components[j]*alp_components[k]*(right1-left1)/f_a**2)
        result -= f_pi**2*momentum_dot(momenta, i, j)*sp.trace(first*second)/2
        result -= alp_components[i]*sp.trace((chi1+chi1.conjugate().T)*meson_matrices[j]*meson_matrices[k])/(8*f_a)
        result += sp.I*f_pi*alp_components[i]*alp_components[j]*sp.trace(meson_matrices[k]*(chi2.conjugate().T-chi2))/(8*f_a**2)
        result += f_pi**2*alp_components[i]*alp_components[j]*alp_components[k]*sp.trace(chi3+chi3.conjugate().T)/(24*f_a**3)
        result -= f_pi*anomaly_mass_squared*theta_derivatives[2]*alp_components[i]*alp_components[j]*sp.trace(meson_matrices[k])/(2*dimension*f_a**2)
        result -= f_pi**2*anomaly_mass_squared*theta_derivatives[1]*theta_derivatives[2]*alp_components[i]*alp_components[j]*alp_components[k]/(dimension*f_a**3)
    return sp.expand(result)


def momentum_dot(momenta, i, j):
    """Use symbolic products or an explicitly supplied symmetric dot-product matrix."""
    return momenta[i, j] if isinstance(momenta, sp.MatrixBase) else momenta[i]*momenta[j]


def pole_overlap(kinetic, mass_squared, squared_mass):
    """Return w(s) and its positive norm when the meson block is invertible.

    K and M² are the source-derived full quadratic matrices, with the ALP last.
    The supplied squared mass is a simple root of det(s K-M²).
    """
    inverse_propagator = squared_mass*kinetic - mass_squared
    overlap = -inverse_propagator[:-1, :-1].inv()*inverse_propagator[:-1, -1]
    direction = overlap.col_join(sp.ones(1, 1))
    norm = (direction.T*kinetic*direction)[0]
    return overlap, sp.sqrt(norm)


def overlap_flavor_matrix(squared_mass, quark_masses, chi_prime, right_0,
                          left_0, theta_prime, f_pi, f_a,
                          anomaly_mass_squared):
    """Apply the explicit meson-operator inverse to the general mixing source.

    quark_masses are 2 B_0 m_q, not mixed physical pole masses.
    """
    dimension = len(quark_masses)
    identity = sp.eye(dimension)
    source = (f_pi/f_a)*(
        -squared_mass*(right_0-left_0)
        + sp.I*(chi_prime.conjugate().T-chi_prime)/2
        - 2*anomaly_mass_squared*theta_prime*identity/dimension
    )
    denominators = [mass-squared_mass for mass in quark_masses]
    trace = sum(source[i, i]/denominators[i] for i in range(dimension))
    trace /= 1 + anomaly_mass_squared*sum(1/value for value in denominators)/dimension
    return sp.Matrix(dimension, dimension, lambda i, j:
        ((source[i, i]-anomaly_mass_squared*trace/dimension)/denominators[i]
         if i == j else source[i, j]/((quark_masses[i]+quark_masses[j])/2-squared_mass)))


def aligned_sources(parent_matrix, chi_derivatives, right_derivatives,
                    left_derivatives, theta_derivatives, f_pi, f_a):
    """Taylor derivatives after the exact symmetric parent-pole chiral map."""
    x = sp.Symbol('x', real=True)
    order = 3
    identity = sp.eye(parent_matrix.rows)
    generator = sp.I*f_a*parent_matrix/(2*f_pi)
    positive = sum(((x*generator)**n/sp.factorial(n) for n in range(order+1)),
                   sp.zeros(parent_matrix.rows))
    negative = sum(((-x*generator)**n/sp.factorial(n) for n in range(order+1)),
                   sp.zeros(parent_matrix.rows))
    def source_series(coefficients):
        return sum((coefficients[n]*x**n/sp.factorial(n) for n in range(len(coefficients))),
                   sp.zeros(parent_matrix.rows))
    mass = negative*source_series(chi_derivatives)*negative
    right = negative*source_series(right_derivatives)*positive-f_a*parent_matrix/(2*f_pi)
    left = positive*source_series(left_derivatives)*negative+f_a*parent_matrix/(2*f_pi)
    def derivatives(matrix, count):
        return [matrix.diff(x, n).subs(x, 0).applyfunc(sp.expand) for n in range(count)]
    theta = list(theta_derivatives)
    theta[1] += f_a*sp.trace(parent_matrix)/(2*f_pi)
    return derivatives(mass, 4), derivatives(right, 2), derivatives(left, 2), theta


def aligned_vector_scalar_sources(parent_matrix, vector, axial, scalar,
                                  pseudoscalar, f_pi, f_a):
    """Transform source derivatives directly in V/A and S/P variables."""
    dimension = parent_matrix.rows
    def source_derivative(n, first, second, current):
        out = sp.zeros(dimension)
        for k in range(n+1):
            for l in range(n-k+1):
                m = n-k-l
                coefficient = (sp.factorial(n)/(sp.factorial(k)*sp.factorial(l)*sp.factorial(m))
                               *(f_a/(2*f_pi))**(k+m))
                if current:
                    plus = ((-sp.I)**k*sp.I**m+sp.I**k*(-sp.I)**m)/2
                    minus = ((-sp.I)**k*sp.I**m-sp.I**k*(-sp.I)**m)/2
                else:
                    plus = ((-sp.I)**(k+m)+sp.I**(k+m))/2
                    minus = sp.I*((-sp.I)**(k+m)-sp.I**(k+m))/2
                out += coefficient*parent_matrix**k*(plus*first[l]+minus*second[l])*parent_matrix**m
        return out.applyfunc(sp.expand)
    vector_result = [source_derivative(n, vector, axial, True) for n in range(2)]
    axial_result = [source_derivative(n, axial, vector, True) for n in range(2)]
    axial_result[0] -= f_a*parent_matrix/(2*f_pi)
    scalar_result = [source_derivative(n, scalar, pseudoscalar, False) for n in range(4)]
    pseudoscalar_result = [source_derivative(n, pseudoscalar, [-p for p in scalar], False)
                           for n in range(4)]
    return vector_result, axial_result, scalar_result, pseudoscalar_result


def aligned_cubic_terms(first_meson, second_meson, first_alp, second_alp,
                        momenta, chi, right, left, theta, f_pi, f_a,
                        anomaly_mass_squared):
    """Explicit aPP, aaP and aaa terms after parent alignment."""
    p01, p02, p12 = momenta[0, 1], momenta[0, 2], momenta[1, 2]
    dimension = sp.Integer(first_meson.rows)
    commutator = lambda a, b: a*b-b*a
    anticommutator = lambda a, b: a*b+b*a
    two_meson = (
        sp.I*(p01-p02)*sp.trace(first_meson*commutator(right[0]+left[0], second_meson))/(4*f_a)
        - sp.trace((chi[1]+chi[1].conjugate().T)*anticommutator(first_meson, second_meson))/(8*f_a)
    )
    one_meson = (
        - sp.I*f_pi*(second_alp*p02*sp.trace((right[0]-left[0])*commutator(right[0]+left[0], first_meson))
                     + first_alp*p01*sp.trace((right[0]-left[0])*commutator(right[0]+left[0], second_meson)))/(2*f_a**2)
        + f_pi*(second_alp*(p01+p12)*sp.trace(first_meson*(right[1]-left[1]))
                + first_alp*(p02+p12)*sp.trace(second_meson*(right[1]-left[1])))/(2*f_a**2)
        + sp.I*f_pi*(second_alp*sp.trace(first_meson*(chi[2].conjugate().T-chi[2]))
                     + first_alp*sp.trace(second_meson*(chi[2].conjugate().T-chi[2])))/(4*f_a**2)
        - f_pi*anomaly_mass_squared*theta[2]*(second_alp*sp.trace(first_meson)
                                               + first_alp*sp.trace(second_meson))/(dimension*f_a**2)
    )
    three_alp = first_alp*second_alp*(
        - f_pi**2*(p01+p02+p12)*sp.trace((right[0]-left[0])*(right[1]-left[1]))/f_a**3
        + f_pi**2*sp.trace(chi[3]+chi[3].conjugate().T)/(4*f_a**3)
        - 6*f_pi**2*anomaly_mass_squared*theta[1]*theta[2]/(dimension*f_a**3)
    )
    return tuple(map(sp.expand, (two_meson, one_meson, three_alp)))


def aligned_component_terms(first_meson, second_meson, momenta, vector, axial,
                            scalar, pseudoscalar, theta, B_0, f_pi, f_a,
                            anomaly_mass_squared):
    """Explicit flavour sums after pole alignment; both daughters have ALP overlap one.

    The supplied sources are the derivatives *after* the chiral map, with
    chi_n=2 B_0(S_n+i P_n), C_R,n=C_V,n+C_A,n. Scalar and pseudoscalar
    matrices are Hermitian. The dimensions and source normalizations agree
    with aligned_cubic_terms; no on-shell mass identities are substituted.
    """
    p01, p02, p12 = momenta[0, 1], momenta[0, 2], momenta[1, 2]
    dimension = first_meson.rows
    indices = range(dimension)
    vector_term = sum(first_meson[i, j]*(vector[0][j, k]*second_meson[k, i]
                       - second_meson[j, k]*vector[0][k, i])
                      for i in indices for j in indices for k in indices)
    scalar_term = sum(scalar[1][i, j]*(first_meson[j, k]*second_meson[k, i]
                      + second_meson[j, k]*first_meson[k, i])
                      for i in indices for j in indices for k in indices)
    first_commutator = sum(axial[0][i, j]*(vector[0][j, k]*first_meson[k, i]
                           - first_meson[j, k]*vector[0][k, i])
                          for i in indices for j in indices for k in indices)
    second_commutator = sum(axial[0][i, j]*(vector[0][j, k]*second_meson[k, i]
                            - second_meson[j, k]*vector[0][k, i])
                           for i in indices for j in indices for k in indices)
    first_axial = sum(first_meson[i, j]*axial[1][j, i] for i in indices for j in indices)
    second_axial = sum(second_meson[i, j]*axial[1][j, i] for i in indices for j in indices)
    pseudoscalar_term = sum((first_meson[i, j]+second_meson[i, j])*pseudoscalar[2][j, i]
                            for i in indices for j in indices)
    axial_product = sum(axial[0][i, j]*axial[1][j, i] for i in indices for j in indices)
    two = sp.I*(p01-p02)*vector_term/(2*f_a)-B_0*scalar_term/(2*f_a)
    one = (-2*sp.I*f_pi*(p02*first_commutator+p01*second_commutator)/f_a**2
           +f_pi*((p01+p12)*first_axial+(p02+p12)*second_axial)/f_a**2
           +f_pi*B_0*pseudoscalar_term/f_a**2
           -f_pi*anomaly_mass_squared*theta[2]*sum(first_meson[i, i]+second_meson[i, i]
             for i in indices)/(dimension*f_a**2))
    three = (-4*f_pi**2*(p01+p02+p12)*axial_product/f_a**3
             +f_pi**2*B_0*sum(scalar[3][i, i] for i in indices)/f_a**3
             -6*f_pi**2*anomaly_mass_squared*theta[1]*theta[2]/(dimension*f_a**3))
    return tuple(map(sp.expand, (two, one, three)))


def decoupled_kaon_component_amplitude(pion_direction, kaon_mixing_6,
                                       kaon_mixing_7, momenta, chi, right, left,
                                       f_pi, f_a):
    """Explicit ds/sd components for the zero-ALP kaon and a mixed daughter."""
    p01, p02, p12 = momenta[0, 1], momenta[0, 2], momenta[1, 2]
    vector = right[0]+left[0]
    axial = right[0]-left[0]
    mass_source = chi[1]+chi[1].conjugate().T
    second_mass_source = chi[2].conjugate().T-chi[2]
    norm = sp.sqrt(kaon_mixing_6**2+kaon_mixing_7**2)
    result = 0
    for row, column, factor in (
        (2, 1, (-kaon_mixing_7-sp.I*kaon_mixing_6)/norm),
        (1, 2, (-kaon_mixing_7+sp.I*kaon_mixing_6)/norm),
    ):
        result += factor*(
            sp.I*(p01-p02)*(vector*pion_direction-pion_direction*vector)[row, column]/(4*f_a)
            - (pion_direction*mass_source+mass_source*pion_direction)[row, column]/(8*f_a)
            - sp.I*f_pi*p02*(axial*vector-vector*axial)[row, column]/(2*f_a**2)
            + f_pi*(p01+p12)*(right[1]-left[1])[row, column]/(2*f_a**2)
            + sp.I*f_pi*second_mass_source[row, column]/(4*f_a**2)
        )
    return sp.expand(result)


def aligned_three_pole_amplitude(meson_basis, kinetic, mass_squared, squared_masses,
                                 momenta, chi_derivatives, right_derivatives,
                                 left_derivatives, theta_derivatives, f_pi, f_a,
                                 anomaly_mass_squared, zero_alp_directions=None):
    """Evaluate the CCWZ amplitude with exact pole alignment and no vertex tensor.

    External momenta are incoming. The parent must have nonzero ALP component
    and invertible meson block. For a daughter with zero ALP component, pass
    its K-normalized meson-component column in zero_alp_directions keyed by
    its leg number (1 or 2). Other poles use w(s) and the positive square-root
    normalization; their overall signs are a choice of state phase.
    """
    zero_alp_directions = zero_alp_directions or {}
    def source_overlap(squared_mass):
        overlap_matrix = overlap_flavor_matrix(
            squared_mass, list(chi_derivatives[0].diagonal()), chi_derivatives[1],
            right_derivatives[0], left_derivatives[0], theta_derivatives[1],
            f_pi, f_a, anomaly_mass_squared)
        overlap = sp.Matrix([sp.trace(generator*overlap_matrix)/2
                             for generator in meson_basis])
        direction = overlap.col_join(sp.ones(1, 1))
        return overlap, sp.sqrt((direction.T*kinetic*direction)[0])

    parent_overlap, parent_norm = source_overlap(squared_masses[0])
    zero_matrix = sp.zeros(meson_basis[0].rows)
    parent_matrix = sum((parent_overlap[i]*meson_basis[i]
                         for i in range(len(meson_basis))), zero_matrix)
    aligned_matrices = [zero_matrix]
    alp_components = [sp.S.One]
    norms = [parent_norm]
    for leg in (1, 2):
        if leg in zero_alp_directions:
            aligned_overlap = zero_alp_directions[leg]
            alp_components.append(sp.S.Zero)
            norms.append(sp.S.One)
        else:
            overlap, norm = source_overlap(squared_masses[leg])
            aligned_overlap = overlap-parent_overlap
            alp_components.append(sp.S.One)
            norms.append(norm)
        aligned_matrices.append(sum((aligned_overlap[i]*meson_basis[i]
                                     for i in range(len(meson_basis))), zero_matrix))
    mass, right, left, theta = aligned_sources(
        parent_matrix, chi_derivatives, right_derivatives, left_derivatives,
        theta_derivatives, f_pi, f_a)
    return sum(aligned_cubic_terms(aligned_matrices[1], aligned_matrices[2],
                                   alp_components[1], alp_components[2],
                                   momenta, mass, right, left, theta,
                                   f_pi, f_a, anomaly_mass_squared))/sp.prod(norms)


def main():
    f_pi, f_a, B_0, m_u, m_d, m_s = sp.symbols(
        'f_pi f_a B_0 m_u m_d m_s', positive=True)
    coupling = sp.Symbol('c_ds', real=True)
    zero = sp.zeros(3)
    vector = sp.Matrix([[0, 0, 0], [0, 0, sp.I*coupling], [0, -sp.I*coupling, 0]])
    chi = [2*B_0*sp.diag(m_u, m_d, m_s), zero, zero, zero]
    first = sp.diag(1, -1, 0)
    second = sp.Matrix([[0, 0, 0], [0, 0, 1], [0, 1, 0]])
    momentum = sp.symbols('p_a p_1 p_2', real=True)
    gram = sp.Matrix([[momentum[i]*momentum[j] for j in range(3)] for i in range(3)])
    terms = aligned_cubic_terms(first, second, 0, 0, gram, chi,
                                [vector, zero], [vector, zero], [0, 0, 0],
                                f_pi, f_a, sp.S.Zero)
    print('Illustrative restricted d-s vector source; all momenta incoming.')
    print('aPP, aaP, aaa terms, with no overall Feynman-rule factor i:')
    for term in terms:
        print(sp.latex(term))



if __name__ == '__main__':
    main()
