"""Vertex rules up to five legs for generic sources S_hat_{r,s}, P_hat_{r,s} (task 2.2).

Run:
    PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_loop_vertices.py
    PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_loop_vertices

Units GeV; conventions of sigma_pipi_loop_algebra.py (polar fields, connection-
free frame, epsilon = f_pi/f_a, all momenta incoming, d_mu -> -i k_mu, vertex =
n-th derivative of the Lagrangian).

Construction 1: direct noncommutative expansion of the Lagrangian with
Phi = sum_j phi_j x_j, a = sum_j phi_j xi_j, sigma = sum phi_j over radial legs,
and extraction of the coefficient of phi_1 ... phi_n.
Construction 2: the closed vertex rules of sigma_pipi_loop_algebra.py (word
level and tree level).
The two are compared as exact identities of cyclic words.  The script also
compares the closed formulas of murayama-radial-up-only-estimate-result.md
with the direct expansion, and checks the dimension of every vertex.  The
field ends (radial, ALP field, meson) defined here are used by
derive_sigma_pipi_loop_large_epsilon.py for the powers of epsilon of the
vertices at large epsilon.  No loop is evaluated.
"""
import math

import sympy as sp

import sigma_pipi_loop_algebra as algebra
from sigma_pipi_loop_algebra import (B0, EPSILON, F_PI, M_SIGMA, QUARTIC, TraceSum, Leg,
                                     symbolic_dot)

MAX_LEGS = 5


def generic_legs(number_of_legs, radial):
    """radial sigma legs first, then generic pseudoscalar legs (x_j, xi_j)."""
    legs = []
    for index in range(number_of_legs):
        label = f'q{index + 1}'
        if index < radial:
            legs.append(Leg(radial=True, momentum={label: 1}))
        else:
            legs.append(Leg(meson=((1, f'x{index + 1}'),),
                            alp=sp.Symbol(f'xi{index + 1}', real=True), momentum={label: 1}))
    return legs


def _noncommutative_to_trace_sum(expression, fields, markers=None, source=None):
    """Multilinear coefficient of prod(fields) of a traced noncommutative
    expression, as a TraceSum.  markers: momentum symbols kappa_j (two per
    term), replaced by the dot product of the corresponding momenta."""
    result = TraceSum()
    expression = sp.expand(expression)
    for term in sp.Add.make_args(expression):
        commuting, ordered = term.args_cnc()
        commuting = sp.Mul(*commuting)
        powers = commuting.as_powers_dict()
        if any(powers.get(field, 0) != 1 for field in fields):
            continue
        coefficient = commuting.subs({field: 1 for field in fields})
        if markers:
            present = [index for index, marker in enumerate(markers) if powers.get(marker, 0)]
            exponents = [powers.get(markers[index], 0) for index in present]
            labels = [index for index, exponent in zip(present, exponents) for _ in range(int(exponent))]
            if len(labels) != 2:
                continue
            coefficient = coefficient.subs({marker: 1 for marker in markers})
            coefficient *= symbolic_dot({f'q{labels[0] + 1}': 1}, {f'q{labels[1] + 1}': 1})
        word = []
        for factor in ordered:
            base, exponent = factor.as_base_exp()
            word += [base.name] * int(exponent)
        if source is not None:
            word = [source] + word
        result.add_term((tuple(word),), coefficient)
    return result


def direct_vertex(number_of_legs, radial, include_lambda=True):
    """Construction 1 for the pattern (radial sigma legs, generic pseudoscalars)."""
    n = number_of_legs
    fields = sp.symbols(f'phi1:{n + 1}')
    markers = sp.symbols(f'kappa1:{n + 1}')
    matrices = [sp.Symbol(f'x{j + 1}', commutative=False) for j in range(n)]
    alp_components = [sp.Symbol(f'xi{j + 1}', real=True) for j in range(n)]
    sigma = sum(fields[:radial], sp.Integer(0))
    alp = sum((alp_components[j] * fields[j] for j in range(radial, n)), sp.Integer(0))
    meson = sum((fields[j] * matrices[j] for j in range(radial, n)), sp.Integer(0))
    derivative_meson = sum((-sp.I * markers[j] * fields[j] * matrices[j] for j in range(radial, n)),
                           sp.Integer(0))
    derivative_alp = sum((-sp.I * markers[j] * fields[j] * alp_components[j] for j in range(radial, n)),
                         sp.Integer(0))
    f_a = F_PI / EPSILON
    result = TraceSum()
    # Source part: B_0 f^2 (sigma/f_a)^r (a/f_a)^s <S_rs cos(Phi/f) + P_rs sin(Phi/f)>.
    for r in range(radial + 1):
        for s in range(n - radial + 1):
            m = n - r - s
            if m < 0:
                continue
            kind = 'S' if m % 2 == 0 else 'P'
            series_coefficient = sp.Integer(-1) ** (m // 2) / (math.factorial(m) * F_PI ** m)
            prefactor = B0 * F_PI ** 2 * (sigma / f_a) ** r * (alp / f_a) ** s * series_coefficient
            matrix = meson ** m if m else sp.Integer(1)
            result = result + _noncommutative_to_trace_sum(prefactor * matrix, fields,
                                                           source=f'{kind}{r}{s}')
    # Meson kinetic part f^2/4 <dU dU^dagger>, U = exp(i Phi/f).
    if radial == 0 and n % 2 == 0:
        def derivative_of_power(power):
            if power == 0:
                return sp.Integer(0)
            return sum((meson ** p * derivative_meson * meson ** (power - 1 - p) for p in range(power)),
                       sp.Integer(0))
        kinetic = 0
        for k in range(1, n):
            l = n - k
            kinetic += ((sp.I / F_PI) ** k / math.factorial(k) * derivative_of_power(k)
                        * (-sp.I / F_PI) ** l / math.factorial(l) * derivative_of_power(l))
        result = result + _noncommutative_to_trace_sum(F_PI ** 2 / 4 * kinetic, fields, markers)
    # ALP kinetic part 1/2 (1 + sigma/f_a)^2 (d a)^2 (commuting; traced as the empty word).
    alp_kinetic = sp.expand(sp.Rational(1, 2) * (1 + sigma / f_a) ** 2 * derivative_alp ** 2)
    scalar = _scalar_multilinear(alp_kinetic, fields, markers)
    # Radial potential -lambda/4 (f_a + sigma)^4 (the quadratic part has no vertex).
    if include_lambda:
        scalar += _scalar_multilinear(sp.expand(-QUARTIC / 4 * (f_a + sigma) ** 4), fields, None)
    if scalar != 0:
        result = result + TraceSum.scalar(scalar)
    return result


def _scalar_multilinear(expression, fields, markers):
    total = 0
    for term in sp.Add.make_args(sp.expand(expression)):
        powers = term.as_powers_dict()
        if any(powers.get(field, 0) != 1 for field in fields):
            continue
        value = term.subs({field: 1 for field in fields})
        if markers:
            labels = [index for index, marker in enumerate(markers)
                      for _ in range(int(powers.get(marker, 0)))]
            if len(labels) != 2:
                continue
            value = value.subs({marker: 1 for marker in markers})
            value *= symbolic_dot({f'q{labels[0] + 1}': 1}, {f'q{labels[1] + 1}': 1})
        total += value
    return sp.expand(total)


def closed_vertex(number_of_legs, radial, eliminate_lambda=False):
    return algebra.vertex(generic_legs(number_of_legs, radial), symbolic_dot, eliminate_lambda)


def compare_constructions():
    """{(n, r): (direct == closed words, closed trees == closed words)}."""
    results = {}
    for n in range(3, MAX_LEGS + 1):
        for radial in range(n + 1):
            direct = direct_vertex(n, radial)
            closed = closed_vertex(n, radial)
            trees = algebra.compact_vertex(generic_legs(n, radial), symbolic_dot,
                                           eliminate_lambda=False).to_trace_sum()
            results[(n, radial)] = ((direct - closed).is_zero(), (trees - closed).is_zero())
    return results


# ---------------------------------------------------------------------------
# Closed formulas as written in murayama-radial-up-only-estimate-result.md
# ---------------------------------------------------------------------------

def record_potential_formula(legs):
    """V_pot of the record: B_0 sum r! s! eps^{r+s} (-1)^{m(m-1)/2}/m! ...
    written without any power of f_pi."""
    n = len(legs)
    return algebra.potential_vertex(legs).scale(F_PI ** (n - 2))


def record_meson_kinetic_formula(legs):
    """V_kin^(4) of the record: 1/48 sum (-k.k) <[,][,]>, without 1/f_pi^2."""
    return algebra.meson_kinetic_vertex(legs, symbolic_dot).scale(F_PI ** 2)


def record_alp_kinetic_formula(legs):
    """V_kin,a of the record with C(sigma) = 2 eps, C(sigma,sigma) = 2 eps^2."""
    radial = sum(1 for leg in legs if leg.radial)
    return algebra.alp_kinetic_vertex(legs, symbolic_dot).scale(F_PI ** radial)


def record_radial_cubic(leg_x, leg_y):
    """V(sigma; x, y) of the record, for generic legs x, y."""
    xi_x, xi_y = leg_x.alp, leg_y.alp
    (_, name_x), = leg_x.meson
    (_, name_y), = leg_y.meson
    dot = symbolic_dot(leg_x.momentum, leg_y.momentum)
    value = (TraceSum.single(('S10', name_x, name_y), -B0 * EPSILON / 2)
             + TraceSum.single(('S10', name_y, name_x), -B0 * EPSILON / 2)
             + TraceSum.single(('P11', name_y), B0 * EPSILON ** 2 * xi_x)
             + TraceSum.single(('P11', name_x), B0 * EPSILON ** 2 * xi_y)
             + TraceSum.single(('S12',), 2 * B0 * EPSILON ** 3 * xi_x * xi_y)
             + TraceSum.scalar(-2 * EPSILON * xi_x * xi_y * dot))
    return value.scale(1 / F_PI)


def record_radial_quartic(leg_x, leg_y):
    """V(sigma, sigma; x, y) = -2 eps^2 xi_x xi_y k_x.k_y / f_pi^2 (up-only)."""
    return TraceSum.scalar(-2 * EPSILON ** 2 * leg_x.alp * leg_y.alp
                           * symbolic_dot(leg_x.momentum, leg_y.momentum) / F_PI ** 2)


def up_only_substitution(trace_sum, highest_radial=MAX_LEGS, highest_alp=MAX_LEGS):
    """S_hat_rs + i P_hat_rs = m_u (-i)^s/s! P_u for r <= 1, (r,s) != (0,0);
    zero for r >= 2; P_hat_00 = 0; S_hat_00 unchanged."""
    up_mass = sp.Symbol('m_u', positive=True)
    mapping = {'P00': None}
    for r in range(highest_radial + 1):
        for s in range(highest_alp + 1):
            if (r, s) == (0, 0):
                continue
            if r >= 2:
                mapping[f'S{r}{s}'] = None
                mapping[f'P{r}{s}'] = None
                continue
            value = up_mass * (-sp.I) ** s / sp.factorial(s)
            mapping[f'S{r}{s}'] = (sp.re(value), 'Pu')
            mapping[f'P{r}{s}'] = (sp.im(value), 'Pu')
    result = TraceSum()
    for key, coefficient in trace_sum.terms.items():
        factor = coefficient
        words = []
        vanishes = False
        for word in key:
            new_word = []
            for name in word:
                if name in mapping:
                    if mapping[name] is None or mapping[name][0] == 0:
                        vanishes = True
                        break
                    factor *= mapping[name][0]
                    new_word.append(mapping[name][1])
                else:
                    new_word.append(name)
            if vanishes:
                break
            words.append(tuple(new_word))
        if not vanishes:
            result.add_term(tuple(words), factor)
    return result


def record_checks():
    """Comparison of the record formulas with the direct expansion."""
    checks = {}
    for n in range(3, MAX_LEGS + 1):
        for radial in range(n + 1):
            legs = generic_legs(n, radial)
            potential = algebra.potential_vertex(legs)
            record = record_potential_formula(legs)
            checks[('V_pot', n, radial)] = {
                'record equals rule': (record - potential).is_zero(),
                'record times f_pi^(2-n) equals rule': (record.scale(F_PI ** (2 - n)) - potential).is_zero()}
    legs = generic_legs(4, 0)
    kinetic = algebra.meson_kinetic_vertex(legs, symbolic_dot)
    record = record_meson_kinetic_formula(legs)
    checks[('V_kin4',)] = {'record equals rule': (record - kinetic).is_zero(),
                           'record over f_pi^2 equals rule': (record.scale(F_PI ** -2) - kinetic).is_zero()}
    for radial in (1, 2):
        legs = generic_legs(2 + radial, radial)
        kinetic = algebra.alp_kinetic_vertex(legs, symbolic_dot)
        record = record_alp_kinetic_formula(legs)
        checks[('V_kin_a', radial)] = {
            'record equals rule': (record - kinetic).is_zero(),
            'record over f_pi^r equals rule': (record.scale(F_PI ** -radial) - kinetic).is_zero()}
    legs = generic_legs(3, 1)
    rule = algebra.drop_vanishing_letters(direct_vertex(3, 1))
    record = record_radial_cubic(legs[1], legs[2])
    checks[('V(sigma;x,y)', 'generic sources')] = {'record equals direct': (record - rule).is_zero()}
    legs = generic_legs(4, 2)
    rule = up_only_substitution(direct_vertex(4, 2))
    record = record_radial_quartic(legs[2], legs[3])
    checks[('V(sigma,sigma;x,y)', 'up-only')] = {'record equals direct': (record - rule).is_zero()}
    for n in (3, 4):
        for radial in range(n + 1):
            direct = up_only_substitution(direct_vertex(n, radial))
            legs = generic_legs(n, radial)
            record = up_only_substitution(record_potential_formula(legs).scale(F_PI ** (2 - n))
                                          + algebra.meson_kinetic_vertex(legs, symbolic_dot)
                                          + algebra.alp_kinetic_vertex(legs, symbolic_dot)
                                          + algebra.radial_self_vertex(legs, eliminate_lambda=False))
            checks[('up-only', n, radial)] = {'corrected record equals direct': (record - direct).is_zero()}
    return checks


# ---------------------------------------------------------------------------
# Field ends
# ---------------------------------------------------------------------------
# A field end is a radial end (sigma), an ALP-field end (the field a, with no
# meson component) or a meson end (a bare U(3) generator).  No pole vector, no
# normalization N and no residue enters a field-end vertex.
FIELD_ENDS = ('radial', 'ALP field', 'meson')


def field_end_leg(kind, index):
    momentum = {f'q{index}': 1}
    if kind == 'radial':
        return Leg(radial=True, momentum=momentum)
    if kind == 'ALP field':
        return Leg(alp=sp.Integer(1), momentum=momentum)
    return Leg(meson=((1, f'G{index}'),), momentum=momentum)


def dimension_check():
    """V^(n) -> c^{4-n} V^(n) under X -> c^{d_X} X."""
    scale = sp.Symbol('c', positive=True)
    results = {}
    for n in range(3, MAX_LEGS + 1):
        for radial in range(n + 1):
            value = algebra.vertex(generic_legs(n, radial), symbolic_dot, eliminate_lambda=True)
            ok = True
            for key, coefficient in value.terms.items():
                dots = {symbol: scale ** 2 * symbol for symbol in sp.sympify(coefficient).free_symbols
                        if symbol.name.startswith('d_')}
                scaled = coefficient.subs({B0: scale * B0, F_PI: scale * F_PI, M_SIGMA: scale * M_SIGMA,
                                           **dots}, simultaneous=True)
                scaled *= scale ** algebra.source_count(key)
                ok &= sp.simplify(scaled - scale ** (4 - n) * coefficient) == 0
            results[(n, radial)] = ok
    return results


def closed_rules_latex():
    lines = [
        r'V^{(n)}_{\rm pot}=B_0f_\pi^{2-n}\sum_{(R,A,P)}r!\,s!\,\epsilon^{r+s}\frac{(-1)^{m(m-1)/2}}{m!}'
        r'\prod_{x\in A}\xi_x\sum_{\pi\in S_m}\langle X_{r,s}\lambda_{\pi_1}\cdots\lambda_{\pi_m}\rangle',
        r'V^{(4)}_{\rm kin}=\frac{1}{48f_\pi^2}\sum_{\pi\in S_4}(-k_{\pi_2}\cdot k_{\pi_4})'
        r'\langle[\lambda_{\pi_1},\lambda_{\pi_2}][\lambda_{\pi_3},\lambda_{\pi_4}]\rangle',
        r'V_{{\rm kin},a}(\sigma;x,y)=-\frac{2\epsilon}{f_\pi}\xi_x\xi_y\,k_x\cdot k_y,\quad '
        r'V_{{\rm kin},a}(\sigma,\sigma;x,y)=-\frac{2\epsilon^2}{f_\pi^2}\xi_x\xi_y\,k_x\cdot k_y',
        r'V(\sigma,\sigma,\sigma)=-6\lambda f_a+\frac{6B_0\epsilon^3}{f_\pi}\langle\widehat S_{3,0}\rangle'
        r'=-\frac{3\epsilon}{f_\pi}\big[m_\sigma^2-\epsilon^2B_0(\langle\widehat S_{1,0}\rangle-2\langle\widehat S_{2,0}\rangle)\big]'
        r'+\frac{6B_0\epsilon^3}{f_\pi}\langle\widehat S_{3,0}\rangle',
        r'V(\sigma,\sigma,\sigma,\sigma)=-6\lambda+\frac{24B_0\epsilon^4}{f_\pi^2}\langle\widehat S_{4,0}\rangle,\qquad '
        r'V(\sigma^5)=\frac{120B_0\epsilon^5}{f_\pi^3}\langle\widehat S_{5,0}\rangle',
    ]
    return lines


def main():
    print('Closed vertex rules (GeV; all momenta incoming):')
    for line in closed_rules_latex():
        print('  ' + line)
    print('\nConstruction 1 (direct expansion) against construction 2 (closed rules), n <= 5:')
    for (n, radial), (direct_ok, tree_ok) in compare_constructions().items():
        print(f'  n={n}, sigma legs={radial}: direct = rule: {direct_ok}; trees = words: {tree_ok}')
    print('\nRecord formulas of murayama-radial-up-only-estimate-result.md against the direct expansion:')
    for key, value in record_checks().items():
        print(f'  {key}: {value}')
    print('\nDimension check V^(n) -> c^(4-n) V^(n):', all(dimension_check().values()))


if __name__ == '__main__':
    main()
