"""Cartesian frame: size of the difference between the Cartesian and the polar one-loop vacuum conditions, relative to the tree of sigma -> pi+ pi- (task 4.5).

Run:
    PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_cartesian_loop_limits.py
    PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_cartesian_loop_limits

OpenSpec change estimate-sigma-pipi-one-loop-murayama-cartesian.  Units GeV;
eps = f_pi/f_a, u = 1/eps.

Vacuum conditions (assumptions of the two estimates):
  Cartesian: one-loop s tadpole = 0, i.e. <Re phi> = f_a;
  polar:     one-loop sigma tadpole = 0, i.e. <|phi|> = f_a.
With phi = f_a + s + i p, |phi| = sqrt((f_a + s)^2 + p^2) = f_a + s + p^2/(2 f_a) + ...,
so at one loop the two definitions of f_a differ by delta f_a = <p^2>/(2 f_a).
The tree is M_0 = -(B_0/(2 f_a)) <S^C_10 {lambda_pi+, lambda_pi-}> in both
parametrizations, so M_0 changes by the factor f_a/(f_a + delta f_a).

Size rule (D1, D3; assumption): <p^2> is a scalar moment with one denominator on
the bare p line, whose D3 scale is m_a; it is replaced by m_a^2/(16 pi^2).  Only
the size is meaningful; the sign is not fixed by the rule.

The large-eps behaviour uses m_a^2 = eps^2 A + 2 c_0/A + ... (root of D_a,
large-eps input record) with the sources fixed.  The leading power is read on
the exact expression; nothing is truncated in the estimate.

No loop integral, scalar loop function, loop coefficient or loop amplitude is computed.
"""
import sympy as sp

import derive_sigma_pipi_cartesian_loop_large_epsilon_vertices as vertices_script
from sigma_pipi_loop_algebra import B0, EPSILON, F_PI, M_A
from sigma_pipi_cartesian_loop_algebra import ALP_MASS_COEFFICIENT, UP_QUARK_MASS

U = vertices_script.INVERSE_EPSILON
DECAY_CONSTANT = sp.Symbol('f_a', positive=True)
RADIAL, ALP_FIELD = sp.symbols('s p', real=True)
ALP_FIELD_SQUARED = sp.Symbol('<p^2>', positive=True)       # one-loop <p^2>, sized by the rule below


def modulus_expansion():
    """|phi| = sqrt((f_a + s)^2 + p^2) through second order in s and p."""
    t = sp.Symbol('t', positive=True)
    modulus = sp.sqrt((DECAY_CONSTANT + t * RADIAL) ** 2 + (t * ALP_FIELD) ** 2)
    return sp.expand(sp.series(modulus, t, 0, 3).removeO().subs(t, 1))


def vacuum_shift():
    """delta f_a = <|phi|> - <Re phi> at one loop, with <s> = 0 (Cartesian condition): the p^2 term of |phi|."""
    expansion = modulus_expansion()
    real_part = DECAY_CONSTANT + RADIAL
    difference = sp.expand(expansion - real_part)
    assert difference.coeff(RADIAL, 2) == 0 and difference.coeff(RADIAL, 1) == 0
    return sp.expand(difference.subs(ALP_FIELD ** 2, ALP_FIELD_SQUARED))


def size_of_alp_field_squared():
    """Size rule: <p^2> -> m_a^2/(16 pi^2) (one denominator, scale m_a of the bare p line, D3)."""
    return M_A ** 2 / (16 * sp.pi ** 2)


LOOP_COUNTING = sp.Symbol('h', positive=True)     # one power per loop: <p^2> is of order h


def exact_tree_factor():
    """M_0(f_a + delta f_a)/M_0(f_a) - 1 with M_0 proportional to 1/f_a, all orders in delta f_a."""
    shift = vacuum_shift().subs(ALP_FIELD_SQUARED, LOOP_COUNTING * size_of_alp_field_squared())
    return sp.simplify(DECAY_CONSTANT / (DECAY_CONSTANT + shift) - 1)


def relative_tree_change():
    """One-loop part of the relative change of M_0: the term of first order in the loop counting h
    (delta f_a is a one-loop quantity; the higher orders in delta f_a belong to two and more loops, where
    the vacuum condition is not specified).  Equal to -delta f_a/f_a."""
    return sp.simplify(sp.diff(exact_tree_factor(), LOOP_COUNTING).subs(LOOP_COUNTING, 0))


def leading_relative_change(alp_mass_squared=None):
    """(power of eps, coefficient) of the relative change of M_0 at large eps, f_a = f_pi/eps."""
    if alp_mass_squared is None:
        alp_mass_squared = vertices_script.alp_pole_root()
    expression = relative_tree_change().subs({M_A: sp.sqrt(alp_mass_squared), DECAY_CONSTANT: F_PI * U})
    series = sp.expand(sp.series(expression, U, 0, 2).removeO())
    return vertices_script.leading_term(series)


def up_only_leading_relative_change():
    """Same for the up-only sources (A = B_0 m_u)."""
    power, coefficient = leading_relative_change()
    return power, sp.factor(coefficient.subs(ALP_MASS_COEFFICIENT, B0 * UP_QUARK_MASS))


def up_only_d4_powers():
    """{label: largest power of eps over the routings} of the excluded D4 terms for the up-only sources, relative
    to M_0 (for the pi leg: dZ itself), read on the exact up-only expressions."""
    import derive_sigma_pipi_cartesian_loop_graphs as graphs
    import derive_sigma_pipi_cartesian_loop_up_only_estimate as up_only
    result = {}
    for group in graphs.GROUPS:
        for row in graphs.table(group):
            if row['D4 status'] != 'excluded (D4)':
                continue
            powers = []
            for value in row['D4 term']:
                scalar = up_only.up_only_scalar(value) * up_only.contribution_factor(group)
                leading = up_only.leading(up_only.series_in_u(sp.expand(scalar), up_only.ROOT_ORDER, 'r = 1'))
                if leading is not None:
                    powers.append(leading[0])
            result[row['label']] = max(powers) if powers else None
    return result


def ratio_scale_identity():
    """eps^4 X/(16 pi^2 f_pi^2) = (eps^2 X)/(16 pi^2 f_a^2) with f_a = f_pi/eps: the generic leading scale of R,
    eps^4 A/(16 pi^2 f_pi^2), is m_a^2/(16 pi^2 f_a^2) at leading order (m_a^2 -> eps^2 A)."""
    scale = sp.Symbol('X', positive=True)
    left = EPSILON ** 4 * scale / (16 * sp.pi ** 2 * F_PI ** 2)
    right = (EPSILON ** 2 * scale) / (16 * sp.pi ** 2 * (F_PI / EPSILON) ** 2)
    return sp.simplify(left - right) == 0


def comparison_with_ratio(power_of_ratio, coefficient_of_ratio=None):
    """(difference of powers: vacuum change minus R, ratio of the leading coefficients if R's is given)."""
    power, coefficient = leading_relative_change()
    if coefficient_of_ratio is None:
        return power - power_of_ratio, None
    return power - power_of_ratio, sp.factor(sp.Abs(coefficient) / sp.Abs(coefficient_of_ratio))


def main():
    print('|phi| through second order:', modulus_expansion())
    print('delta f_a = <|phi|> - <Re phi> (one loop):', vacuum_shift())
    print('relative change of M_0, all orders in delta f_a (h counts loops):', exact_tree_factor())
    print('one-loop relative change of M_0 (size rule <p^2> -> m_a^2/(16 pi^2)):', relative_tree_change())
    print('large eps, generic sources (power of eps, coefficient):', leading_relative_change())
    print('large eps, up-only sources:', up_only_leading_relative_change())
    print('powers above the generic leading power of R (eps^4):', comparison_with_ratio(4)[0])
    up_only_ratio = -B0 ** 2 * UP_QUARK_MASS ** 2 / (8 * sp.pi ** 2 * F_PI ** 2 * sp.Symbol('m_P', positive=True) ** 2)
    power, coefficient = up_only_leading_relative_change()
    print('up-only: powers above the leading power of R (eps^2):', power - 2,
          '; ratio of the leading coefficients:', sp.factor(sp.Abs(coefficient) / sp.Abs(up_only_ratio)))
    print('generic leading scale of R equals m^2/(16 pi^2 f_a^2):', ratio_scale_identity())
    print('up-only D4 terms, largest power relative to M_0:', up_only_d4_powers())


if __name__ == '__main__':
    main()
