"""Cartesian frame: invariance of the traces of the one-loop sigma -> pi+ pi- estimate under a constant axial redefinition (task 4.2).

Run:
    PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_cartesian_loop_frame_invariance.py
    PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_cartesian_loop_frame_invariance

OpenSpec change estimate-sigma-pipi-one-loop-murayama-cartesian.  Written anew
for the Cartesian Taylor variables y = s/f_a and x = p/f_a:

    M(s, p) = sum_{r,q} (S^C_{r,q} + i P^C_{r,q}) y^r x^q,   Theta(s, p) = sum_{r,q} Theta_{r,q} y^r x^q.

Redefinition (decision taken with the brainstormer; assumption on its form):
g_R = g_L^dagger = exp(i x kappa), with a constant Hermitian 3x3 matrix kappa
and [kappa, Q] = 0; s is not rotated.  Starting from the frame without
derivative couplings of s and p (C_R = C_L = 0), the new frame has the raw
source M' = g_R M g_L^dagger, the connections C'_R = kappa, C'_L = -kappa along
p, and Theta' = Theta - x <kappa>.

Dressing (transport along x at fixed s): d_x W = i C W with W(0) = 1, then
    M_hat = W_R^dagger M W_L,   Theta_hat = Theta + int_0^x <(C_R - C_L)/2>.
Redefinitions that depend on s are not covered.

Checks (exact, noncommuting symbols, no matrices):
 (1) the dressed coefficients S^C_hat_{r,q}, P^C_hat_{r,q} (r, q <= 3) agree in the two frames;
 (2) every kept T_Gamma (termwise routing maximum of stage B, D4 terms removed,
     D3 scales) is unchanged when the dressed coefficients of the new frame are
     inserted; the masses, N(m_a^2), M^{-1}_{m_a^2}(C), lambda_pi+- and the
     Piecewise conditions are held opaque, and the dressed coefficients on which
     they depend are listed and covered by (1);
 (3) negative control: inserting the raw coefficients of the new frame, to first
     order in kappa, changes the kept T_Gamma listed by the script;
 (4) Theta_hat_{0,1} is unchanged while the raw Theta_{0,1} shifts by -<kappa>.
No loop is evaluated.
"""
import functools
import sys

import sympy as sp

import sigma_pipi_cartesian_loop_algebra as cartesian
import derive_sigma_pipi_cartesian_loop_traces as traces
import derive_sigma_pipi_cartesian_loop_tree as cartesian_tree
from sigma_pipi_loop_algebra import TraceSum

ORDER = 3                       # Taylor coefficients with r, q <= ORDER
REDEFINITION_LETTER = 'kappa'   # the constant Hermitian matrix of the redefinition
TRACE_OF_REDEFINITION = sp.Symbol('<kappa>', real=True)
FIRST_ORDER = sp.Symbol('t', real=True)     # bookkeeping parameter for the first order in kappa


# ---------------------------------------------------------------------------
# Noncommutative polynomials in Hermitian letters
# ---------------------------------------------------------------------------

class Polynomial:
    """Linear combination of words of Hermitian letters: {word (tuple of names): complex coefficient}."""

    def __init__(self, terms=None):
        self.terms = {}
        for word, value in (terms or {}).items():
            self.add(word, value)

    @staticmethod
    def letter(name, value=1):
        return Polynomial({(name,): value})

    @staticmethod
    def unit(value=1):
        return Polynomial({(): value})

    def add(self, word, value):
        self.terms[word] = sp.expand(self.terms.get(word, 0) + value)
        if self.terms[word] == 0:
            del self.terms[word]

    def __add__(self, other):
        result = Polynomial(self.terms)
        for word, value in other.terms.items():
            result.add(word, value)
        return result

    def __sub__(self, other):
        return self + other.scale(-1)

    def scale(self, factor):
        return Polynomial({word: factor * value for word, value in self.terms.items()})

    def __mul__(self, other):
        result = Polynomial()
        for word_one, value_one in self.terms.items():
            for word_two, value_two in other.terms.items():
                result.add(word_one + word_two, value_one * value_two)
        return result

    def dagger(self):
        """All letters Hermitian: reverse each word and conjugate its coefficient."""
        return Polynomial({tuple(reversed(word)): sp.conjugate(value) for word, value in self.terms.items()})

    def hermitian_part(self):
        return (self + self.dagger()).scale(sp.Rational(1, 2))

    def antihermitian_part(self):
        """P in X = S + i P: (X - X^dagger)/(2i)."""
        return (self - self.dagger()).scale(1 / (2 * sp.I))

    def is_zero(self):
        return not self.terms

    def __eq__(self, other):
        return (self - other).is_zero()


class Series:
    """Truncated double series sum_{r,q <= ORDER} c_{r,q} y^r x^q with Polynomial coefficients."""

    def __init__(self, coefficients=None):
        self.coefficients = {key: value for key, value in (coefficients or {}).items()
                             if key[0] <= ORDER and key[1] <= ORDER and not value.is_zero()}

    def __add__(self, other):
        result = dict(self.coefficients)
        for key, value in other.coefficients.items():
            result[key] = result.get(key, Polynomial()) + value
        return Series(result)

    def __mul__(self, other):
        result = {}
        for (r1, q1), one in self.coefficients.items():
            for (r2, q2), two in other.coefficients.items():
                key = (r1 + r2, q1 + q2)
                if key[0] <= ORDER and key[1] <= ORDER:
                    result[key] = result.get(key, Polynomial()) + one * two
        return Series(result)

    def dagger(self):
        return Series({key: value.dagger() for key, value in self.coefficients.items()})

    def coefficient(self, r, q):
        return self.coefficients.get((r, q), Polynomial())


def source_series():
    """M(y, x) = sum (S^C_{r,q} + i P^C_{r,q}) y^r x^q with generic letters (including P^C_00, P^C_10)."""
    return Series({(r, q): Polynomial.letter(cartesian.source_letter('S', r, q))
                   + Polynomial.letter(cartesian.source_letter('P', r, q), sp.I)
                   for r in range(ORDER + 1) for q in range(ORDER + 1)})


def exponential_along_alp(matrix, sign=1):
    """exp(i sign x matrix) as a series in x."""
    power = Polynomial.unit()
    coefficients = {}
    for n in range(ORDER + 1):
        coefficients[(0, n)] = power.scale(sp.Rational(1, sp.factorial(n)))
        power = power * matrix.scale(sign * sp.I)
    return Series(coefficients)


def transport(connection):
    """W with d_x W = i C W and W(0) = 1, for a connection C(y, x) given as a Series in x (y-independent here):
    (n + 1) W_{n+1} = i sum_j C_j W_{n-j}."""
    w = {0: Polynomial.unit()}
    for n in range(ORDER):
        total = Polynomial()
        for j in range(n + 1):
            total = total + connection.coefficient(0, j) * w[n - j]
        w[n + 1] = total.scale(sp.I / (n + 1))
    return Series({(0, n): value for n, value in w.items()})


def redefinition_letter():
    return Polynomial.letter(REDEFINITION_LETTER)


def frames():
    """{'original': (M, C_R, C_L), 'new': (M', C'_R, C'_L)} as Series."""
    source = source_series()
    kappa = redefinition_letter()
    g_right = exponential_along_alp(kappa, 1)           # g_R = exp(i x kappa)
    g_left_dagger = exponential_along_alp(kappa, 1)     # g_L^dagger = (exp(-i x kappa))^dagger = exp(i x kappa)
    new_source = g_right * source * g_left_dagger
    zero = Series()
    return {'original': (source, zero, zero),
            'new': (new_source, Series({(0, 0): kappa}), Series({(0, 0): kappa.scale(-1)}))}


def dressed(source, connection_right, connection_left):
    """M_hat = W_R^dagger M W_L."""
    w_right = transport(connection_right)
    w_left = transport(connection_left)
    return w_right.dagger() * source * w_left


def coefficient_pairs(series):
    """{(kind, r, q): Polynomial} with S = Hermitian part and P = (X - X^dagger)/(2i) of each Taylor coefficient."""
    result = {}
    for r in range(ORDER + 1):
        for q in range(ORDER + 1):
            value = series.coefficient(r, q)
            result[('S', r, q)] = value.hermitian_part()
            result[('P', r, q)] = value.antihermitian_part()
    return result


@functools.lru_cache(maxsize=None)
def coefficients():
    """{'dressed original', 'dressed new', 'raw new', 'raw original'}: {(kind, r, q): Polynomial}."""
    data = frames()
    return {'raw original': coefficient_pairs(data['original'][0]),
            'raw new': coefficient_pairs(data['new'][0]),
            'dressed original': coefficient_pairs(dressed(*data['original'])),
            'dressed new': coefficient_pairs(dressed(*data['new']))}


def dressed_coefficients_agree():
    """Check (1): number of coefficients compared and whether all agree."""
    data = coefficients()
    agree = all(data['dressed new'][key] == data['dressed original'][key] for key in data['dressed original'])
    original_is_raw = all(data['dressed original'][key] == Polynomial.letter(cartesian.source_letter(*key))
                          for key in data['dressed original'])
    return len(data['dressed original']), agree and original_is_raw


def raw_first_order_shift():
    """{(kind, r, q): first-order change of the raw coefficient}; expected -{kappa, P^C_{r,q-1}} for S
    and {kappa, S^C_{r,q-1}} for P (q >= 1), zero for q = 0."""
    data = coefficients()
    shifts = {}
    for key, value in data['raw new'].items():
        difference = value - data['raw original'][key]
        first = Polynomial({word: value for word, value in difference.terms.items()
                            if word.count(REDEFINITION_LETTER) == 1})
        shifts[key] = first
    return shifts


def expected_first_order_shift(kind, r, q):
    if q == 0:
        return Polynomial()
    kappa = redefinition_letter()
    if kind == 'S':
        other = Polynomial.letter(cartesian.source_letter('P', r, q - 1))
        return (kappa * other + other * kappa).scale(-1)
    other = Polynomial.letter(cartesian.source_letter('S', r, q - 1))
    return kappa * other + other * kappa


# ---------------------------------------------------------------------------
# Insertion into the kept traces
# ---------------------------------------------------------------------------

def letter_key(name):
    kind, r, q = cartesian.source_indices(name)
    return kind, r, q


def substitute(trace_sum, replacements):
    """Replace every Cartesian source letter X by replacements[X] (a Polynomial) in every trace;
    other letters (lambda_pi+-, M^{-1}_{m_a^2}(C)) are kept."""
    result = TraceSum()
    for key, value in trace_sum.terms.items():
        products = [((), value)]   # list of (words so far, coefficient)
        for word in key:
            expanded = Polynomial.unit()
            for name in word:
                if cartesian.is_cartesian_source(name):
                    expanded = expanded * replacements[name]
                else:
                    expanded = expanded * Polynomial.letter(name)
            products = [(words + (new_word,), coefficient * new_value)
                        for words, coefficient in products for new_word, new_value in expanded.terms.items()]
        for words, coefficient in products:
            result.add_term(words, coefficient)
    return result.simplified()


def first_order_variation(trace_sum, shifts):
    """d/dt of T(X + t dX) at t = 0: each occurrence of a source letter replaced once by its shift."""
    result = TraceSum()
    for key, value in trace_sum.terms.items():
        for word_index, word in enumerate(key):
            for position, name in enumerate(word):
                if not cartesian.is_cartesian_source(name):
                    continue
                shift = shifts[name]
                for new_word, new_value in shift.terms.items():
                    replaced = word[:position] + new_word + word[position + 1:]
                    result.add_term(key[:word_index] + (replaced,) + key[word_index + 1:], value * new_value)
    return cartesian.apply_charge_rule(cartesian.apply_assumptions(result)).simplified()


def kept_values():
    """{label: T_Gamma} of every kept graph (1PI X^max and self-energy dZ^max), D4 terms removed."""
    return {entry['label']: entry['maximum']['value'] for entry in traces.entries()}


def letters_in(trace_sum):
    return {name for name in trace_sum.letters() if cartesian.is_cartesian_source(name)}


def entering_letters():
    """Source letters of the kept traces, of M_0 and of the pole data (A, kappa_sigma, C, L)."""
    letters = set()
    for value in kept_values().values():
        letters |= letters_in(value)
    letters |= letters_in(cartesian_tree.tree_amplitude())
    return letters | pole_data_letters()


def pole_data_letters():
    """Dressed coefficients on which m_a, m_sigma, N(m_a^2) and M^{-1}_{m_a^2}(C) depend:
    A = B_0(<S^C_10> - 2<S^C_02>), kappa_sigma = B_0(<S^C_10> - 2<S^C_20>), C = B_0 P^C_01,
    L(X) = B_0 {S^C_00, X} + (m_1^2/3)<X> 1 (Theta_hat_1 = 0 by assumption 6)."""
    letters = letters_in(cartesian_tree.alp_mass_coefficient()) | letters_in(cartesian_tree.kappa_sigma())
    return letters | {'PC01', 'SC00'}


def coefficients_are_scalars():
    """The coefficients of every kept T_Gamma contain no source letter: the masses, N, the size factors
    and the Piecewise conditions are scalar symbols, so only the products of traces are transformed."""
    for value in kept_values().values():
        for coefficient in value.terms.values():
            names = {symbol.name for symbol in sp.sympify(coefficient).free_symbols}
            if any(cartesian.is_cartesian_source(name) for name in names):
                return False
    return True


def within_covered_range(letters):
    return all(cartesian.source_indices(name)[1] <= ORDER and cartesian.source_indices(name)[2] <= ORDER
               for name in letters)


def replacements(frame):
    data = coefficients()[frame]
    return {cartesian.source_letter(*key): value for key, value in data.items()}


@functools.lru_cache(maxsize=None)
def trace_checks():
    """{label: (unchanged with dressed coefficients of the new frame, changed with raw coefficients
    to first order in kappa, contains a source with q >= 1, q >= 1 sources present)}."""
    dressed_new = replacements('dressed new')
    shifts = {cartesian.source_letter(*key): value for key, value in raw_first_order_shift().items()}
    result = {}
    for label, value in kept_values().items():
        unchanged = (substitute(value, dressed_new) - value).is_zero()
        changed = not first_order_variation(value, shifts).is_zero()
        alp_letters = sorted(name for name in letters_in(value) if cartesian.source_indices(name)[2] >= 1)
        result[label] = (unchanged, changed, bool(alp_letters), alp_letters)
    return result


def connection_trace():
    """<(C'_R - C'_L)/2> of the new frame: the Polynomial (C'_R - C'_L)/2, whose trace is read off
    from its single-letter word (kappa,)."""
    kappa = redefinition_letter()
    half_difference = (kappa - kappa.scale(-1)).scale(sp.Rational(1, 2))
    assert set(half_difference.terms) == {(REDEFINITION_LETTER,)}
    return half_difference.terms[(REDEFINITION_LETTER,)] * TRACE_OF_REDEFINITION


def theta_check():
    """Check (4): (raw Theta'_{0,1} - Theta_{0,1}, dressed Theta'_{0,1} - Theta_{0,1}, other dressed
    coefficients equal).  Theta' = Theta - x <kappa>; Theta_hat = Theta' + int_0^x <(C'_R - C'_L)/2>,
    the integrand being constant, so the integral adds x <(C'_R - C'_L)/2> to the x^1 coefficient."""
    theta = {(r, q): sp.Symbol(f'Theta_{r}{q}', real=True) for r in range(ORDER + 1) for q in range(ORDER + 1)}
    raw_new = dict(theta)
    raw_new[(0, 1)] = theta[(0, 1)] - TRACE_OF_REDEFINITION
    dressed_new = dict(raw_new)
    dressed_new[(0, 1)] = raw_new[(0, 1)] + connection_trace()
    others_equal = all(dressed_new[key] == theta[key] for key in theta)
    return sp.expand(raw_new[(0, 1)] - theta[(0, 1)]), sp.expand(dressed_new[(0, 1)] - theta[(0, 1)]), others_equal


def letters_with_vanishing_shift():
    """Source letters with q >= 1 whose first-order raw shift vanishes by the assumptions
    P^C_00 = P^C_10 = 0 (matrix identities): S^C_01 and S^C_11 (shift -{kappa, P^C_00}, -{kappa, P^C_10})."""
    result = set()
    for key, shift in raw_first_order_shift().items():
        if key[2] == 0:
            continue
        surviving = [word for word in shift.terms if not any(name in cartesian.VANISHING_LETTERS for name in word)]
        if not surviving:
            result.add(cartesian.source_letter(*key))
    return result


def negative_control():
    """(changed graphs, graphs with a source of q >= 1, graphs with a source of q >= 1 whose shift does not
    vanish).  The changed set equals the last one."""
    checks = trace_checks()
    vanishing = letters_with_vanishing_shift()
    changed = {label for label, value in checks.items() if value[1]}
    with_alp = {label for label, value in checks.items() if value[2]}
    with_shifted = {label for label, value in checks.items() if set(value[3]) - vanishing}
    return changed, with_alp, with_shifted


def main():
    count, agree = dressed_coefficients_agree()
    print(f'(1) dressed coefficients identical in the two frames ({count} coefficients, r, q <= {ORDER}):', agree)
    shifts = raw_first_order_shift()
    print('    raw first-order shifts equal -{kappa, P_(r,q-1)} and {kappa, S_(r,q-1)}:',
          all(shifts[key] == expected_first_order_shift(*key) for key in shifts))
    letters = entering_letters()
    print('Source letters entering (kept traces, M_0, pole data):', sorted(letters, key=cartesian.source_indices))
    print('    all within r, q <=', ORDER, ':', within_covered_range(letters))
    print('    pole data depend on:', sorted(pole_data_letters(), key=cartesian.source_indices))
    print('Coefficients of the kept traces contain no source letter:', coefficients_are_scalars())
    checks = trace_checks()
    print(f'(2) every kept T_Gamma unchanged with the dressed coefficients of the new frame ({len(checks)} graphs):',
          all(value[0] for value in checks.values()))
    changed, with_alp, with_shifted = negative_control()
    print(f'(3) changed with raw coefficients (first order in kappa): {len(changed)} graphs: {sorted(changed)}')
    print(f'    graphs with a source of q >= 1: {len(with_alp)}: {sorted(with_alp)}')
    print('    q >= 1 letters whose raw shift vanishes (P^C_00 = P^C_10 = 0):', sorted(letters_with_vanishing_shift()))
    print('    changed set = graphs with a q >= 1 source of nonvanishing shift:', changed == with_shifted)
    for label in sorted(with_alp - changed):
        print(f'    {label} unchanged; its q >= 1 sources: {checks[label][3]}')
    print('(4) (raw Theta_01 shift, dressed Theta_01 shift, other dressed Theta equal):', theta_check())


if __name__ == '__main__':
    main()
