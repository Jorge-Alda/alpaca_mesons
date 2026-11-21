"""Check the finite-mixing pole-vector normalization of latex/finite_mixing.tex.

Invocation:
    .venv/bin/python scripts/check_finite_pole_vector_normalization.py

Checks eq:finite-pole-vector-normalization against eq:finite-kinetic-blocks,
eq:finite-pole-norm, eq:finite-invariant-daughter-pole-vector and
eq:finite-field-expansion, symbolically and for generic Hermitian 3x3 matrices.

Conventions (Alda:2026nrd): <A> = tr A; U(3) basis lambda_0 = sqrt(2/3) 1 and
the Gell-Mann matrices lambda_1..8, with <lambda_i lambda_j> = 2 delta_ij.
The axial source c = c_q^A is dimensionless and Hermitian, including its trace
part. The matrix X stands for M_a^{-1}(C) evaluated at the pole s = m_r^2; it is
dimensionless (C and M_a both have mass dimension two) and is treated as a
generic Hermitian matrix, since the identity only uses Hermiticity of X.
f_pi and f_a have mass dimension one; only their ratio enters.
The pole vector is v = (v_i, v_hat) with
    v_i = f_pi <lambda_i (c + X)> / (f_a sqrt N),  v_hat = 1/sqrt N,
and lambda_r = sum_i v_i lambda_i. No expansion in f_pi/f_a is made.
"""
import sympy as sp

PION_DECAY_CONSTANT, ALP_DECAY_CONSTANT = sp.symbols('f_pi f_a', positive=True)
POLE_NORM = sp.Symbol('N', positive=True)

ALL_GENERATORS = tuple(range(9))
# Generators commuting with the quark charge diag(2/3,-1/3,-1/3): pi^0, eta_8, eta_0, K^0 and Kbar^0 directions.
NEUTRAL_GENERATORS = (0, 3, 6, 7, 8)


def u3_generators():
    """Return lambda_0..lambda_8 with <lambda_i lambda_j> = 2 delta_ij."""
    gell_mann = [sp.sqrt(sp.Rational(2, 3)) * sp.eye(3)]
    entries = [
        {(0, 1): 1, (1, 0): 1}, {(0, 1): -sp.I, (1, 0): sp.I}, {(0, 0): 1, (1, 1): -1},
        {(0, 2): 1, (2, 0): 1}, {(0, 2): -sp.I, (2, 0): sp.I},
        {(1, 2): 1, (2, 1): 1}, {(1, 2): -sp.I, (2, 1): sp.I},
    ]
    for entry in entries:
        matrix = sp.zeros(3)
        for (row, column), value in entry.items():
            matrix[row, column] = value
        gell_mann.append(matrix)
    gell_mann.append(sp.diag(1, 1, -2) / sp.sqrt(3))
    return gell_mann


def hermitian_matrix(name, generator_indices):
    """Generic Hermitian matrix sum_i (a_i/2) lambda_i over the given generators, real a_i."""
    generators = u3_generators()
    components = {index: sp.Symbol(f'{name}_{index}', real=True) for index in generator_indices}
    return sum((components[index] * generators[index] / 2 for index in generator_indices), sp.zeros(3)), components


def trace(matrix):
    return sp.expand(matrix.trace())


def kinetic_matrix(axial_source, generator_indices):
    """K of eq:finite-kinetic-blocks, meson fields first, ALP last."""
    generators = u3_generators()
    ratio = PION_DECAY_CONSTANT / ALP_DECAY_CONSTANT
    size = len(generator_indices)
    kinetic = sp.eye(size + 1)
    overlaps = [trace(axial_source * generators[index]) for index in generator_indices]
    for position, overlap in enumerate(overlaps):
        kinetic[position, size] = kinetic[size, position] = -ratio * overlap
    kinetic[size, size] = 1 + ratio**2 * sum(overlap**2 for overlap in overlaps)
    return kinetic


def normalization_left_side(meson_matrix, axial_source, alp_component):
    """Left-hand side of eq:finite-pole-vector-normalization, with 1/sqrt(N) -> alp_component."""
    ratio = PION_DECAY_CONSTANT / ALP_DECAY_CONSTANT
    return sp.expand(trace(meson_matrix**2) / 2 - 2 * ratio * alp_component * trace(axial_source * meson_matrix)
                     + alp_component**2 * (1 + 2 * ratio**2 * trace(axial_source**2)))


def check_quadratic_form(generator_indices):
    """Step 1: v^T K v equals the normalization left side for arbitrary v_i, v_hat. Returns residual."""
    generators = u3_generators()
    axial_source, _ = hermitian_matrix('c', generator_indices)
    meson_components = sp.symbols(f'v_0:{len(generator_indices)}', real=True)
    alp_component = sp.Symbol('v_a', real=True)
    vector = sp.Matrix(list(meson_components) + [alp_component])
    quadratic_form = sp.expand((vector.T * kinetic_matrix(axial_source, generator_indices) * vector)[0, 0])
    meson_matrix = sum((component * generators[index] for component, index in zip(meson_components, generator_indices)), sp.zeros(3))
    return sp.simplify(quadratic_form - normalization_left_side(meson_matrix, axial_source, alp_component))


def charged_source_residual():
    """Residual of Step 1 if c has charged components but only neutral meson fields are kept.

    K uses sum over neutral i of <c lambda_i>^2, the normalization uses the full <c^2>.
    Returns (residual, expected 2 (f_pi/f_a)^2 v_hat^2 (c_1^2+c_2^2+c_4^2+c_5^2)/2)."""
    generators = u3_generators()
    axial_source, components = hermitian_matrix('c', ALL_GENERATORS)
    meson_components = sp.symbols(f'v_0:{len(NEUTRAL_GENERATORS)}', real=True)
    alp_component = sp.Symbol('v_a', real=True)
    vector = sp.Matrix(list(meson_components) + [alp_component])
    quadratic_form = sp.expand((vector.T * kinetic_matrix(axial_source, NEUTRAL_GENERATORS) * vector)[0, 0])
    meson_matrix = sum((component * generators[index] for component, index in zip(meson_components, NEUTRAL_GENERATORS)), sp.zeros(3))
    residual = sp.simplify(quadratic_form - normalization_left_side(meson_matrix, axial_source, alp_component))
    ratio = PION_DECAY_CONSTANT / ALP_DECAY_CONSTANT
    expected = -ratio**2 * alp_component**2 * sum(components[index]**2 for index in (1, 2, 4, 5))
    return residual, expected


def pole_vector(axial_source, overlap_inverse, generator_indices, pole_norm):
    """v_i and v_hat of eq:finite-invariant-daughter-pole-vector for a given N."""
    generators = u3_generators()
    ratio = PION_DECAY_CONSTANT / ALP_DECAY_CONSTANT
    meson = [ratio / sp.sqrt(pole_norm) * trace(generators[index] * (axial_source + overlap_inverse)) for index in generator_indices]
    return meson, 1 / sp.sqrt(pole_norm)


def check_field_expansion(generator_indices):
    """Step 2: sum_i v_i lambda_i = 2 f_pi (c+X)/(f_a sqrt N). Returns the residual matrix."""
    generators = u3_generators()
    axial_source, _ = hermitian_matrix('c', generator_indices)
    overlap_inverse, _ = hermitian_matrix('X', generator_indices)
    meson, _ = pole_vector(axial_source, overlap_inverse, generator_indices, POLE_NORM)
    meson_matrix = sum((component * generators[index] for component, index in zip(meson, generator_indices)), sp.zeros(3))
    expected = 2 * PION_DECAY_CONSTANT / (ALP_DECAY_CONSTANT * sp.sqrt(POLE_NORM)) * (axial_source + overlap_inverse)
    return sp.simplify(meson_matrix - expected)


def normalization_with_free_norm(generator_indices, axial_is_zero=False):
    """Step 3: left side of the normalization with v substituted and N kept free.

    Returns (left side, N * left side, 1 + 2 f_pi^2 <X^2>/f_a^2)."""
    generators = u3_generators()
    axial_source = sp.zeros(3) if axial_is_zero else hermitian_matrix('c', generator_indices)[0]
    overlap_inverse, _ = hermitian_matrix('X', generator_indices)
    meson, alp_component = pole_vector(axial_source, overlap_inverse, generator_indices, POLE_NORM)
    meson_matrix = sum((component * generators[index] for component, index in zip(meson, generator_indices)), sp.zeros(3))
    left_side = normalization_left_side(meson_matrix, axial_source, alp_component)
    ratio = PION_DECAY_CONSTANT / ALP_DECAY_CONSTANT
    pole_norm_formula = 1 + 2 * ratio**2 * trace(overlap_inverse**2)
    return left_side, sp.expand(sp.simplify(POLE_NORM * left_side)), pole_norm_formula


def check_normalization(generator_indices, axial_is_zero=False):
    """Step 4: residual of the normalization after inserting N of eq:finite-pole-norm."""
    left_side, _, pole_norm_formula = normalization_with_free_norm(generator_indices, axial_is_zero)
    return sp.simplify(left_side.subs(POLE_NORM, pole_norm_formula) - 1)


def check_direct_quadratic_form(generator_indices):
    """Step 5: v^T K v - 1 directly from K and v, with N from eq:finite-pole-norm."""
    axial_source, _ = hermitian_matrix('c', generator_indices)
    overlap_inverse, _ = hermitian_matrix('X', generator_indices)
    ratio = PION_DECAY_CONSTANT / ALP_DECAY_CONSTANT
    pole_norm_formula = 1 + 2 * ratio**2 * trace(overlap_inverse**2)
    meson, alp_component = pole_vector(axial_source, overlap_inverse, generator_indices, pole_norm_formula)
    vector = sp.Matrix(meson + [alp_component])
    return sp.simplify((vector.T * kinetic_matrix(axial_source, generator_indices) * vector)[0, 0] - 1)


def main():
    for label, generator_indices in (('full U(3) basis, generic Hermitian c and X', ALL_GENERATORS),
                                     ('neutral generators 0,3,6,7,8, charge-commuting c and X', NEUTRAL_GENERATORS)):
        print(f'== {label}')
        residual = check_quadratic_form(generator_indices)
        print('Step 1  v^T K v - [normalization left side], arbitrary v:', residual)
        assert residual == 0
        residual = check_field_expansion(generator_indices)
        print('Step 2  sum_i v_i lambda_i - 2 f_pi (c+X)/(f_a sqrt N):', list(residual))
        assert residual == sp.zeros(3)
        _, scaled_left_side, pole_norm_formula = normalization_with_free_norm(generator_indices)
        print('Step 3  N * left side, N free:', sp.factor(scaled_left_side))
        print('        minus (1 + 2 f_pi^2 <X^2>/f_a^2):', sp.simplify(scaled_left_side - pole_norm_formula))
        assert sp.simplify(scaled_left_side - pole_norm_formula) == 0
        residual = check_normalization(generator_indices)
        print('Step 4  left side - 1 with N of eq:finite-pole-norm:', residual)
        assert residual == 0
        residual = check_normalization(generator_indices, axial_is_zero=True)
        print('Step 4b c = 0: (1/2)<lambda_r^2> + 1/N - 1:', residual)
        assert residual == 0
        residual = check_direct_quadratic_form(generator_indices)
        print('Step 5  v^T K v - 1 directly:', residual)
        assert residual == 0
    residual, expected = charged_source_residual()
    print('Caveat  c with charged components, only neutral fields kept, Step 1 residual:', residual)
    assert sp.simplify(residual - expected) == 0
    print('Identity holds.')


if __name__ == '__main__':
    main()
