"""External legs, kinematics, charge rule and tree amplitude of sigma -> pi+ pi- (tasks 1.1, 1.2).

Run:
    PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_loop_tree.py
    PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_loop_tree

Units GeV; conventions of sigma_pipi_loop_algebra.py.  All momenta are
incoming: p_sigma + p_1 + p_2 = 0, with p_1 the momentum of the leg
(lambda_{pi+}, 0) and p_2 that of (lambda_{pi-}, 0).  The physical decay has
outgoing momenta -p_1, -p_2; the dot products below are the same in both
conventions.

Contents:
* the exact external dot products from momentum conservation;
* the charge rule: for letters with [Q, X_k] = q_k X_k, <X_1...X_n> = 0 when
  q_1 + ... + q_n != 0; M_s maps Q-commuting matrices to Q-commuting matrices,
  so M^{-1}_{m_a^2}(C) is neutral (assumption B1: all sources commute with Q);
* the tree M_0 = V(sigma; pi+, pi-) from the vertex rule, and its agreement with
  derive_murayama_sigma_pp.py for the up-only sources, by exact substitution.
No loop is evaluated.
"""
import sympy as sp

import sigma_pipi_loop_algebra as algebra
from sigma_pipi_loop_algebra import (B0, EPSILON, F_PI, M_PI, M_SIGMA, TraceSum, CompactSum,
                                     letter, make_dot, meson_leg, radial_leg)

SIGMA_MOMENTUM, PLUS_MOMENTUM, MINUS_MOMENTUM = 'ps', 'p1', 'p2'


def external_dot_products():
    """Solve p_sigma + p_1 + p_2 = 0 with p_sigma^2 = m_sigma^2, p_i^2 = m_pi^2."""
    d12, ds1, ds2 = sp.symbols('d12 ds1 ds2')
    equations = [
        # p_sigma^2 = (p_1 + p_2)^2
        sp.Eq(M_SIGMA ** 2, 2 * M_PI ** 2 + 2 * d12),
        # p_sigma . p_1 = -(p_1 + p_2) . p_1
        sp.Eq(ds1, -(M_PI ** 2 + d12)),
        sp.Eq(ds2, -(M_PI ** 2 + d12)),
    ]
    solution = sp.solve(equations, [d12, ds1, ds2], dict=True)[0]
    return {frozenset((PLUS_MOMENTUM, MINUS_MOMENTUM)): solution[d12],
            frozenset((SIGMA_MOMENTUM, PLUS_MOMENTUM)): solution[ds1],
            frozenset((SIGMA_MOMENTUM, MINUS_MOMENTUM)): solution[ds2],
            frozenset((SIGMA_MOMENTUM,)): M_SIGMA ** 2,
            frozenset((PLUS_MOMENTUM,)): M_PI ** 2,
            frozenset((MINUS_MOMENTUM,)): M_PI ** 2}


EXTERNAL_DOTS = external_dot_products()


def external_basis_dot(first, second):
    return EXTERNAL_DOTS[frozenset((first, second))]


def external_legs():
    """sigma (radial, no meson component: P_10 = 0, <S_11> = 0), pi+ and pi-."""
    return {'sigma': radial_leg({SIGMA_MOMENTUM: 1}),
            'pi+': meson_leg('l+', {PLUS_MOMENTUM: 1}),
            'pi-': meson_leg('l-', {MINUS_MOMENTUM: 1})}


def tree_amplitude():
    """M_0 = V(sigma; pi+, pi-) from the general vertex rule, with the
    assumptions and the charge rule applied."""
    legs = external_legs()
    value = algebra.vertex([legs['sigma'], legs['pi+'], legs['pi-']], make_dot(external_basis_dot))
    return algebra.apply_charge_rule(algebra.apply_assumptions(value)).simplified()


def tree_amplitude_compact():
    """M_0 = -(B_0 eps/(2 f_pi)) <S_10 {lambda_{pi+}, lambda_{pi-}}>."""
    result = CompactSum()
    result.add_term((('P', letter('S10'), ('A', letter('l+'), letter('l-'))),), -B0 * EPSILON / (2 * F_PI))
    return result


TREE_LATEX = (r'M_0=V(\sigma;\pi^+,\pi^-)=-\frac{B_0\epsilon}{2f_\pi}'
              r'\big\langle\widehat S_{1,0}\{\lambda_{\pi^+},\lambda_{\pi^-}\}\big\rangle')


# ---------------------------------------------------------------------------
# Charge rule and neutrality of M^{-1}(C)
# ---------------------------------------------------------------------------

def _rewrite(words, rules):
    """Apply adjacent-pair rewrite rules {(x, y): [(k, replacement word)]} on
    linear (not cyclic) words until no rule applies."""
    changed = True
    while changed:
        changed = False
        new = {}
        for word, value in words.items():
            for position in range(len(word) - 1):
                pair = (word[position], word[position + 1])
                if pair in rules:
                    for factor, replacement in rules[pair]:
                        target = word[:position] + tuple(replacement) + word[position + 2:]
                        new[target] = new.get(target, 0) + factor * value
                    changed = True
                    break
            else:
                new[word] = new.get(word, 0) + value
        words = {word: sp.expand(value) for word, value in new.items() if sp.expand(value) != 0}
    return words


def charge_rule_proof(length=4):
    """Return (derivation, trace_of_commutator) for a word X_1...X_n with
    [Q, X_k] = q_k X_k:
      derivation: [Q, W] - (q_1 + ... + q_n) W written as words (must be empty);
      trace_of_commutator: <Q W> - <W Q> as cyclic words (must vanish).
    Together: (q_1 + ... + q_n) <W> = <[Q, W]> = 0."""
    names = [f'X{k + 1}' for k in range(length)]
    charges = sp.symbols(f'q1:{length + 1}')
    word = tuple(names)
    # [Q, W] = Q W - W Q; move Q through each letter with Q X_k = X_k Q + q_k X_k.
    rules = {('Q', name): [(1, (name, 'Q')), (charge, (name,))] for name, charge in zip(names, charges)}
    commutator = _rewrite({('Q',) + word: 1}, rules)
    commutator[word + ('Q',)] = commutator.get(word + ('Q',), 0) - 1
    commutator[word] = commutator.get(word, 0) - sum(charges)
    derivation = {w: v for w, v in ((w, sp.expand(v)) for w, v in commutator.items()) if v != 0}
    trace = TraceSum.single(('Q',) + word) - TraceSum.single(word + ('Q',))
    return derivation, trace.simplified().terms


def pole_operator_commutes_with_charge():
    """[Q, M_s(X)] - M_s([Q, X]) for M_s(X) = B_0{S_00, X} + m_1^2/3 <X> 1 - s X,
    with [Q, S_00] = 0.  Returns the remainder of the matrix part (words) and
    of the trace part (the coefficient of the unit matrix), both must vanish."""
    pole_mass, singlet = sp.symbols('s m_1^2')
    rules = {('Q', 'S00'): [(1, ('S00', 'Q'))]}

    def apply_pole_operator(words):
        out = {}
        for word, value in words.items():
            for new_word, factor in (((('S00',) + word), B0), (word + ('S00',), B0), (word, -pole_mass)):
                out[new_word] = out.get(new_word, 0) + factor * value
        return out

    def commutator_with_q(words):
        out = {}
        for word, value in words.items():
            out[('Q',) + word] = out.get(('Q',) + word, 0) + value
            out[word + ('Q',)] = out.get(word + ('Q',), 0) - value
        return out

    x = {('X',): 1}
    left = commutator_with_q(apply_pole_operator(x))
    right = apply_pole_operator(commutator_with_q(x))
    difference = {word: left.get(word, 0) - right.get(word, 0) for word in set(left) | set(right)}
    # Normal order: move Q to the right of S00 (they commute); compare words.
    difference = _rewrite({w: v for w, v in difference.items() if v != 0}, rules)
    # Trace part: [Q, <X> 1] = 0, and <[Q, X]> = 0 by cyclicity.
    trace_part = (TraceSum.single(('Q', 'X')) - TraceSum.single(('X', 'Q'))).scale(singlet / 3)
    return difference, trace_part.simplified().terms


# ---------------------------------------------------------------------------
# Up-only comparison with derive_murayama_sigma_pp.py
# ---------------------------------------------------------------------------

UP_ONLY_RULES = {('l+', 'l-'): [(2, ('Pu',))], ('l-', 'l+'): [(2, ('Pd',))],
                 ('Pu', 'Pu'): [(1, ('Pu',))], ('Pd', 'Pd'): [(1, ('Pd',))],
                 ('Pu', 'Pd'): [], ('Pd', 'Pu'): []}
UP_ONLY_SINGLE_TRACES = {('Pu',): 1, ('Pd',): 1}


def up_only_trace_value(trace_sum):
    """Evaluate traces with the up-only relations lambda_{pi+} lambda_{pi-} = 2 P_u,
    lambda_{pi-} lambda_{pi+} = 2 P_d, P_u^2 = P_u, P_d^2 = P_d, P_u P_d = 0,
    <P_u> = <P_d> = 1 (no explicit matrices)."""
    total = 0
    for key, value in trace_sum.terms.items():
        product = value
        for word in key:
            words = _rewrite({tuple(word): 1}, UP_ONLY_RULES)
            # cyclic closure: rewrite once more after rotation to the canonical word
            trace_value = 0
            for reduced_word, factor in words.items():
                rotated = _rewrite({algebra.canonical_word(reduced_word): 1}, UP_ONLY_RULES)
                for final_word, final_factor in rotated.items():
                    if final_word not in UP_ONLY_SINGLE_TRACES:
                        raise ValueError(f'unreduced trace {final_word}')
                    trace_value += factor * final_factor * UP_ONLY_SINGLE_TRACES[final_word]
            product *= trace_value
        total += product
    return sp.expand(total)


def up_only_tree_comparison():
    """(M_0 with S_10 = m_u P_u, -charged vertex of derive_murayama_sigma_pp)."""
    import derive_murayama_sigma_pp as reference
    (ref_b0, ref_mu, ref_fpi, ref_fphi, _), _, _, charged, _ = reference.calculate()
    up_mass = sp.Symbol('m_u', positive=True)
    tree = tree_amplitude()
    substituted = TraceSum()
    for key, value in tree.terms.items():
        words = tuple(tuple('Pu' if name == 'S10' else name for name in word) for word in key)
        factor = up_mass ** sum(word.count('S10') for word in key)
        substituted.add_term(words, factor * value)
    ours = up_only_trace_value(substituted)
    theirs = (-charged).subs({ref_fphi: F_PI / EPSILON, ref_b0: B0, ref_fpi: F_PI, ref_mu: up_mass})
    return sp.simplify(ours), sp.simplify(theirs)


ASSUMPTIONS = [
    ('No tadpoles: P_hat_00 = 0 and <S_hat_01> = 0', 'radial-decay-polar-cartesian-general-sources.md, assumption 1'),
    ('No radial mixing at quadratic order: P_hat_10 = 0 and <S_hat_11> = 0',
     'radial-decay-polar-cartesian-general-sources.md, assumption 2'),
    ('Anomalous singlet term independent of sigma and a', 'radial-decay-polar-cartesian-general-sources.md, assumption 3'),
    ('Connection-free polar frame, leading U(3) chiral order', 'radial-decay-polar-cartesian-general-sources.md, assumption 4'),
    ('Regular simple poles, M_s invertible at the ALP-like pole',
     'radial-decay-polar-cartesian-general-sources.md, assumption 5'),
    ('Theta_hat_1 = 0, so C = B_0 P_hat_01', 'murayama-radial-crossing-implicit.md (Model); frame independence in murayama-sigma-pipi-loop-frame-invariance.md'),
    ('All S_hat_rs, P_hat_rs and S_hat_00 commute with Q', 'alp-lo-charged-meson-decoupling.md'),
    ('f_pi independent of sigma (no sigma-meson kinetic coupling)', 'stated here (B4)'),
]


def main():
    print('External legs (all momenta incoming): sigma radial; pi+ = (lambda_{pi+}, 0); pi- = (lambda_{pi-}, 0).')
    print('External dot products:')
    for key, value in EXTERNAL_DOTS.items():
        print(f'  {"*".join(sorted(key)) if len(key) == 2 else sorted(key)[0] + "^2"} = {sp.latex(value)}')
    derivation, trace = charge_rule_proof()
    print('\nCharge rule: [Q, W] - (sum q_k) W =', derivation or 0, '; <QW> - <WQ> =', trace or 0)
    difference, trace_part = pole_operator_commutes_with_charge()
    print('[Q, M_s(X)] - M_s([Q, X]) =', difference or 0, '(matrix part);', trace_part or 0, '(trace part)')
    print('  => M_s preserves Q-commuting matrices; M^{-1}_{m_a^2}(C) is neutral.')
    tree = tree_amplitude()
    compact = tree_amplitude_compact()
    print('\nTree amplitude from the vertex rule:', tree.terms)
    print('Compact form equal:', (compact.to_trace_sum() - tree).is_zero())
    print('  ' + TREE_LATEX)
    ours, theirs = up_only_tree_comparison()
    print(f'\nUp-only: M_0 = {sp.latex(ours)}; derive_murayama_sigma_pp: {sp.latex(theirs)}; equal: {sp.simplify(ours - theirs) == 0}')
    print('\nAssumptions and the records that state them:')
    for text, source in ASSUMPTIONS:
        print(f'  - {text}  [{source}]')


if __name__ == '__main__':
    main()
