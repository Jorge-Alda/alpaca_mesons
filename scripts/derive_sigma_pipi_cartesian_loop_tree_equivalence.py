"""Equivalence of the Cartesian and polar tree amplitudes of sigma -> pi+ pi- (task 2.2).

Run:
    PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_cartesian_loop_tree_equivalence.py
    PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_cartesian_loop_tree_equivalence

OpenSpec change estimate-sigma-pipi-one-loop-murayama-cartesian.  This is the
only comparison with the polar results in the derivation.  Units GeV.

Method.
 1. Relation between the coefficients.  The source M_hat is one function on
    field space.  With x = sigma/f_a, y = a/f_a, X = s/f_a, Y = p/f_a, the exact
    map 1 + X + i Y = (1 + x) e^{i y} gives x = sqrt((1+X)^2 + Y^2) - 1 and
    y = atan(Y/(1+X)).  Inserting the Taylor series of x, y into
    sum (S_{r,q} + i P_{r,q}) x^r y^q and collecting X^R Y^Q gives
    S^C_{R,Q} + i P^C_{R,Q} as real linear combinations of the polar
    S_{r,q} + i P_{r,q}; since the polar coefficients are Hermitian and the
    combinations are real, S^C (P^C) is the same combination of the S (P).
    Only the first order is needed for the tree.
 2. Cartesian on-shell vertex V(s; pi+, pi-) from
    derive_sigma_pipi_cartesian_loop_tree.py (Cartesian rules only), with the
    coefficient relation inserted, against the polar tree
        M_0 = -(B_0 eps/(2 f_pi)) <S_{1,0} {lambda_{pi+}, lambda_{pi-}}>
    taken as a stated input from murayama-sigma-pipi-loop-tree-and-legs.md
    (no polar vertex code is imported).

Remark: the charged legs have no ALP component (charge conservation, a
fact), so neither the stationarity condition nor the pole equation enters
this comparison.  No off-shell, quartic or quintic comparison is made.  No
loop is evaluated.
"""
import sympy as sp

import derive_sigma_pipi_cartesian_loop_tree as cartesian_tree
from sigma_pipi_loop_algebra import B0, EPSILON, F_PI, TraceSum

RADIAL_RATIO = sp.Symbol('X', real=True)        # s / f_a
ALP_RATIO = sp.Symbol('Y', real=True)           # p / f_a


def polar_symbol(kind, radial, alp):
    """Formal stand-in for the polar coefficient (letter 'S{r}{s}' or 'P{r}{s}')."""
    return sp.Symbol(f'{kind}{radial}{alp}', real=True)


def polar_fields_in_cartesian(order):
    """Taylor series of x = sigma/f_a and y = a/f_a in X, Y to total order `order`."""
    scale = sp.Symbol('scale', positive=True)
    radial = sp.sqrt((1 + RADIAL_RATIO) ** 2 + ALP_RATIO ** 2) - 1
    angle = sp.atan(ALP_RATIO / (1 + RADIAL_RATIO))
    result = []
    for expression in (radial, angle):
        series = sp.series(expression.subs({RADIAL_RATIO: scale * RADIAL_RATIO, ALP_RATIO: scale * ALP_RATIO}),
                           scale, 0, order + 1).removeO()
        result.append(sp.expand(series.subs(scale, 1)))
    return result


def map_check(order):
    """(1 + x) e^{i y} - (1 + X + i Y) vanishes to total order `order` in X, Y."""
    scale = sp.Symbol('scale', positive=True)
    radial, angle = polar_fields_in_cartesian(order)
    expression = (1 + radial) * sp.exp(sp.I * angle) - (1 + RADIAL_RATIO + sp.I * ALP_RATIO)
    series = sp.series(expression.subs({RADIAL_RATIO: scale * RADIAL_RATIO, ALP_RATIO: scale * ALP_RATIO}),
                       scale, 0, order + 1).removeO()
    return sp.expand(series) == 0


def cartesian_coefficients_from_polar(order):
    """{(R, Q): (S^C_{R,Q}, P^C_{R,Q})} as linear combinations of the polar symbols, R + Q <= order."""
    radial, angle = polar_fields_in_cartesian(order)
    scale = sp.Symbol('scale', positive=True)
    source = 0
    for r in range(order + 1):
        for q in range(order + 1 - r):
            source += (polar_symbol('S', r, q) + sp.I * polar_symbol('P', r, q)) * radial ** r * angle ** q
    source = sp.expand(source.subs({RADIAL_RATIO: scale * RADIAL_RATIO, ALP_RATIO: scale * ALP_RATIO}))
    source = sp.expand(sum(source.coeff(scale, k) for k in range(order + 1)))
    polynomial = sp.Poly(source, RADIAL_RATIO, ALP_RATIO)
    result = {}
    for r in range(order + 1):
        for q in range(order + 1 - r):
            coefficient = sp.expand(polynomial.coeff_monomial(RADIAL_RATIO ** r * ALP_RATIO ** q))
            real, imaginary = coefficient.as_real_imag()
            result[(r, q)] = (sp.expand(real), sp.expand(imaginary))
    return result


def letter_mapping(order=1):
    """{Cartesian letter: [(coefficient, polar letter)]} from cartesian_coefficients_from_polar."""
    mapping = {}
    for (r, q), parts in cartesian_coefficients_from_polar(order).items():
        for kind, expression in zip('SP', parts):
            terms = []
            for symbol in expression.free_symbols:
                terms.append((expression.coeff(symbol), symbol.name))
            mapping[f'{kind}C{r}{q}'] = terms
    return mapping


def substitute_letters(trace_sum, mapping):
    """Replace each mapped letter by its linear combination of letters."""
    result = TraceSum()
    for key, value in trace_sum.terms.items():
        expansions = [((), value)]
        for word in key:
            new = []
            options = [mapping.get(name, [(1, name)]) for name in word]
            partial = [((), 1)]
            for option in options:
                partial = [(prefix + (name,), factor * coefficient)
                           for prefix, factor in partial for coefficient, name in option]
            for words, factor in expansions:
                for single, coefficient in partial:
                    new.append((words + (single,), factor * coefficient))
            expansions = new
        for words, factor in expansions:
            result.add_term(words, factor)
    return result.simplified()


def polar_tree_input():
    """Polar M_0 = -(B_0 eps/(2 f_pi)) <S_{1,0} {lambda_{pi+}, lambda_{pi-}}> (stated input)."""
    return (TraceSum.single(('S10', 'l+', 'l-'), -B0 * EPSILON / (2 * F_PI))
            + TraceSum.single(('S10', 'l-', 'l+'), -B0 * EPSILON / (2 * F_PI)))


def tree_equivalence():
    """(Cartesian on-shell V(s; pi+, pi-) in polar coefficients) - (polar M_0)."""
    cartesian_value = cartesian_tree.tree_amplitude()
    return substitute_letters(cartesian_value, letter_mapping(order=1)) - polar_tree_input()


def main():
    print('Exact map checked to second order:', map_check(2))
    print('First-order relation between Cartesian and polar coefficients:')
    for (r, q), (real, imaginary) in cartesian_coefficients_from_polar(1).items():
        print(f'  S^C_{r}{q} = {real},  P^C_{r}{q} = {imaginary}')
    print('Cartesian V(s; pi+, pi-) =', cartesian_tree.TREE_LATEX)
    print('Cartesian minus polar tree, after the coefficient relation:', tree_equivalence().simplified().terms or 0)


if __name__ == '__main__':
    main()
