#!/usr/bin/env python
"""Flavor conventions of the even-parity axial-vector nonet A_{mu nu}.

Run: .venv/bin/python scripts/derive_axial_nonet_conventions.py

Conventions (mirroring the vector nonet of memory/rcht/foundations):
- generators T_a = lambda_a/sqrt(2) (a=1..8), T_0 = 1/sqrt(3), so <T_i T_j> = delta_ij;
- A_{mu nu} = sum_i T_i A^i_{mu nu}; the Lorentz indices are carried as a label only;
- ideal isoscalars f1_N = (u ubar + d dbar)/sqrt(2) with T_N = diag(1,1,0)/sqrt(2)
  and f1_S = s sbar with T_S = diag(0,0,1) (positive strange phase);
- nonideal mixing A_I = R A_phys, i.e. A_phys = R^T A_I, with
  R = ((cos theta, sin theta), (-sin theta, cos theta)) acting on (f1_N, f1_S);
  the physical generators are T_phys = R^T applied to (T_N, T_S) as a column.
No approximation is made: all statements are exact flavor-matrix identities.
Nothing is written to disk; derived matrices are printed with sympy.latex.
"""
import sympy as sp


def gell_mann_matrices():
    """Return the eight Gell-Mann matrices lambda_1..lambda_8."""
    i = sp.I
    return [
        sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]]),
        sp.Matrix([[0, -i, 0], [i, 0, 0], [0, 0, 0]]),
        sp.Matrix([[1, 0, 0], [0, -1, 0], [0, 0, 0]]),
        sp.Matrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]]),
        sp.Matrix([[0, 0, -i], [0, 0, 0], [i, 0, 0]]),
        sp.Matrix([[0, 0, 0], [0, 0, 1], [0, 1, 0]]),
        sp.Matrix([[0, 0, 0], [0, 0, -i], [0, i, 0]]),
        sp.diag(1, 1, -2) / sp.sqrt(3),
    ]


def nonet_generators():
    """Return [T_0, T_1, ..., T_8] with T_0 = 1/sqrt(3) and T_a = lambda_a/sqrt(2)."""
    return [sp.eye(3) / sp.sqrt(3)] + [lam / sp.sqrt(2) for lam in gell_mann_matrices()]


def trace_gram_matrix(generators):
    """Return the matrix of flavor traces <T_i T_j>."""
    size = len(generators)
    return sp.Matrix(size, size, lambda row, col: sp.trace(generators[row] * generators[col]))


def axial_field_symbols():
    """Return the physical-basis (ideal isoscalar) axial-nonet field symbols."""
    names = 'a1_0 a1_p a1_m f1_N f1_S K1_p K1_0 K1_m K1bar_0'
    return dict(zip(names.split(), sp.symbols(names)))


def axial_matrix_explicit(fields):
    """Axial nonet in the ideal isospin form, same layout as the vector matrix."""
    f = fields
    return sp.Matrix([
        [(f['a1_0'] + f['f1_N']) / sp.sqrt(2), f['a1_p'], f['K1_p']],
        [f['a1_m'], (-f['a1_0'] + f['f1_N']) / sp.sqrt(2), f['K1_0']],
        [f['K1_m'], f['K1bar_0'], f['f1_S']],
    ])


def octet_singlet_components(fields):
    """Components A^0..A^8 written through the charge/ideal fields.

    The off-diagonal fields obey A^1 -+ i A^2 = sqrt(2) a1^{+/-}, and similarly
    for (4,5) -> K1^{+/-} and (6,7) -> K1^0, K1bar^0. The isoscalars follow from
    A^8 = (f1_N - sqrt(2) f1_S)/sqrt(3), A^0 = (sqrt(2) f1_N + f1_S)/sqrt(3).
    """
    f = fields
    s2, s3 = sp.sqrt(2), sp.sqrt(3)
    return [
        (s2 * f['f1_N'] + f['f1_S']) / s3,
        (f['a1_p'] + f['a1_m']) / s2,
        sp.I * (f['a1_p'] - f['a1_m']) / s2,
        f['a1_0'],
        (f['K1_p'] + f['K1_m']) / s2,
        sp.I * (f['K1_p'] - f['K1_m']) / s2,
        (f['K1_0'] + f['K1bar_0']) / s2,
        sp.I * (f['K1_0'] - f['K1bar_0']) / s2,
        (f['f1_N'] - s2 * f['f1_S']) / s3,
    ]


def ideal_generators():
    """Return T_N = diag(1,1,0)/sqrt(2) and T_S = diag(0,0,1)."""
    return sp.diag(1, 1, 0) / sp.sqrt(2), sp.diag(0, 0, 1)


def mixing_rotation(angle):
    """Return R with A_I = R A_phys acting on the (f1_N, f1_S) column."""
    return sp.Matrix([[sp.cos(angle), sp.sin(angle)], [-sp.sin(angle), sp.cos(angle)]])


def alda_exponent_comparison():
    """Compare the exponent of u in Alda:2026nrd and in the RChT records.

    Alda:2026nrd: Phi_Alda = sum lambda_P P, u = exp(i Phi_Alda/(2 f_pi)).
    RChT: Phi_R = sum T_P P with T_P = lambda_P/sqrt(2), U = u^2 = exp(i sqrt(2) Phi_R/F),
    hence u = exp(i Phi_R/(sqrt(2) F)). Returns the exponent differences for
    F = f_pi, for u and for U, using generic coefficients for all nine generators
    (lambda_0 = sqrt(2/3) 1 so that <lambda_i lambda_j> = 2 delta_ij).
    """
    f_pi = sp.Symbol('f_pi', positive=True)
    pseudoscalars = sp.symbols('P0:9', real=True)
    lambdas = [sp.sqrt(sp.Rational(2, 3)) * sp.eye(3)] + gell_mann_matrices()
    phi_alda = sum((lam * field for lam, field in zip(lambdas, pseudoscalars)), sp.zeros(3))
    phi_rcht = sum((T * field for T, field in zip(nonet_generators(), pseudoscalars)), sp.zeros(3))
    u_exponent_alda = sp.I * phi_alda / (2 * f_pi)
    u_exponent_rcht = sp.I * phi_rcht / (sp.sqrt(2) * f_pi)
    big_u_exponent_rcht = sp.I * sp.sqrt(2) * phi_rcht / f_pi
    return dict(
        lambda_gram=trace_gram_matrix(lambdas),
        u_exponent_difference=(u_exponent_alda - u_exponent_rcht).applyfunc(sp.simplify),
        big_u_exponent_difference=(2 * u_exponent_alda - big_u_exponent_rcht).applyfunc(sp.simplify),
    )


def derive():
    """Build all convention matrices and the quantities checked by the tests."""
    generators = nonet_generators()
    fields = axial_field_symbols()
    explicit = axial_matrix_explicit(fields)
    components = octet_singlet_components(fields)
    from_generators = sum((T * c for T, c in zip(generators, components)), sp.zeros(3))
    T_N, T_S = ideal_generators()
    isovector_and_kaon = explicit.subs({fields['f1_N']: 0, fields['f1_S']: 0})
    ideal_projection = isovector_and_kaon + T_N * fields['f1_N'] + T_S * fields['f1_S']
    theta = sp.Symbol('theta', real=True)
    R = mixing_rotation(theta)
    f1_1285, f1_1420 = sp.symbols('f1_1285 f1_1420')
    physical_fields = sp.Matrix([f1_1285, f1_1420])
    ideal_fields_from_physical = R * physical_fields
    T_light = sp.simplify(T_N * R[0, 0] + T_S * R[1, 0])
    T_heavy = sp.simplify(T_N * R[0, 1] + T_S * R[1, 1])
    isoscalar_ideal = T_N * ideal_fields_from_physical[0] + T_S * ideal_fields_from_physical[1]
    isoscalar_physical = T_light * f1_1285 + T_heavy * f1_1420
    return dict(
        generators=generators,
        generator_gram=trace_gram_matrix(generators),
        explicit_matrix=explicit,
        matrix_from_generators_minus_explicit=(from_generators - explicit).applyfunc(sp.simplify),
        ideal_gram=trace_gram_matrix([T_N, T_S]),
        ideal_projection_minus_explicit=(ideal_projection - explicit).applyfunc(sp.simplify),
        ideal_singlet_overlaps=sp.Matrix([sp.trace(T_N * generators[0]), sp.trace(T_S * generators[0])]),
        rotation=R,
        rotation_orthogonality=(R.T * R - sp.eye(2)).applyfunc(sp.simplify),
        rotation_determinant=sp.simplify(R.det()),
        physical_generators=(T_light, T_heavy),
        physical_gram=trace_gram_matrix([T_light, T_heavy]).applyfunc(sp.simplify),
        physical_fields_from_ideal=(R.T * ideal_fields_from_physical).applyfunc(sp.simplify),
        isoscalar_rotation_invariance=(isoscalar_ideal - isoscalar_physical).applyfunc(sp.simplify),
        alda=alda_exponent_comparison(),
    )


def main():
    result = derive()
    for label in ('generator_gram', 'explicit_matrix', 'ideal_gram', 'ideal_singlet_overlaps',
                  'rotation', 'physical_gram', 'physical_fields_from_ideal'):
        print(f'{label}:\n{sp.latex(result[label])}\n')
    T_light, T_heavy = result['physical_generators']
    print(f'T_f1(1285):\n{sp.latex(T_light)}\n')
    print(f'T_f1(1420):\n{sp.latex(T_heavy)}\n')
    zero_checks = ('matrix_from_generators_minus_explicit', 'ideal_projection_minus_explicit',
                   'rotation_orthogonality', 'isoscalar_rotation_invariance')
    for label in zero_checks:
        assert result[label] == sp.zeros(*result[label].shape), label
    assert result['generator_gram'] == sp.eye(9)
    assert result['ideal_gram'] == sp.eye(2)
    assert result['physical_gram'] == sp.eye(2)
    assert result['rotation_determinant'] == 1
    assert result['alda']['lambda_gram'] == 2 * sp.eye(9)
    assert result['alda']['u_exponent_difference'] == sp.zeros(3)
    assert result['alda']['big_u_exponent_difference'] == sp.zeros(3)
    print('Assertions passed: nonet normalization, ideal basis, rotation, Alda exponent dictionary.')


if __name__ == '__main__':
    main()
