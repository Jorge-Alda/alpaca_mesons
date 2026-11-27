"""Cartesian vertex rules up to five legs for generic sources S^C_{r,s}, P^C_{r,s} (task 2.1).

Run:
    PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_cartesian_loop_vertices.py
    PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_cartesian_loop_vertices

OpenSpec change estimate-sigma-pipi-one-loop-murayama-cartesian.  Units GeV;
conventions of sigma_pipi_cartesian_loop_algebra.py (Cartesian fields
phi = f_a + s + i p, eps = f_pi/f_a, all momenta incoming, d_mu -> -i k_mu,
vertex = n-th derivative of the Lagrangian).

Construction 1: direct noncommutative expansion of the Cartesian Lagrangian
(source part in the exponential form <M_hat U^dagger + U M_hat^dagger>, V_phi,
kinetic terms) with s = sum phi_j over radial legs, p = sum xi_j phi_j,
Phi = sum phi_j x_j, and extraction of the coefficient of phi_1 ... phi_n.
Construction 2: the closed rules of sigma_pipi_cartesian_loop_algebra.py
(word level and compact trees).

Checks (all exact): agreement of the two constructions for every pattern
n = 3, 4, 5 with r = 0..n radial legs and generic pseudoscalar legs (each with
both lambda_x and xi_x); relabelling symmetry; dimensional scaling c^{4-n};
up-only limit (closed rules after the D-A4 substitution against the direct
expansion with M_hat = m_u P_u phi^*/f_a); no momentum dependence except in
the four-meson kinetic vertex (in particular the V_phi vertices, and no vertex
from the kinetic term of s and p).  No polar vertex is used.  No loop is
evaluated.
"""
import sympy as sp

import sigma_pipi_cartesian_loop_algebra as cartesian
from sigma_pipi_loop_algebra import B0, F_PI, M_SIGMA, TraceSum, symbolic_dot

MAX_LEGS = 5
PATTERNS = [(n, radial) for n in range(3, MAX_LEGS + 1) for radial in range(n + 1)]


def closed_vertex(number_of_legs, radial, eliminate_lambda=False):
    return cartesian.vertex(cartesian.generic_legs(number_of_legs, radial), symbolic_dot, eliminate_lambda)


def compare_constructions():
    """{(n, r): (direct == closed words, closed trees == closed words)}."""
    results = {}
    for n, radial in PATTERNS:
        direct = cartesian.direct_vertex(n, radial)
        closed = closed_vertex(n, radial)
        trees = cartesian.compact_vertex(cartesian.generic_legs(n, radial), symbolic_dot,
                                         eliminate_lambda=False).to_trace_sum()
        results[(n, radial)] = ((direct - closed).is_zero(), (trees - closed).is_zero())
    return results


# ---------------------------------------------------------------------------
# Relabelling symmetry
# ---------------------------------------------------------------------------

def relabel(trace_sum, permutation):
    """Apply the leg relabelling j -> permutation[j] (1-based) to letters x_j, components xi_j and
    dot products d_qa_qb."""
    symbol_map = {}
    for old, new in permutation.items():
        symbol_map[sp.Symbol(f'xi{old}', real=True)] = sp.Symbol(f'xi{new}', real=True)
    labels = sorted(permutation)
    for first in labels:
        for second in labels:
            if first <= second:
                old = sp.Symbol('d_' + '_'.join(sorted((f'q{first}', f'q{second}'))), real=True)
                new = sp.Symbol('d_' + '_'.join(sorted((f'q{permutation[first]}', f'q{permutation[second]}'))),
                                real=True)
                symbol_map[old] = new
    result = TraceSum()
    for key, value in trace_sum.terms.items():
        words = tuple(tuple(f'x{permutation[int(name[1:])]}' if name.startswith('x') and name[1:].isdigit()
                            else name for name in word) for word in key)
        result.add_term(words, sp.sympify(value).subs(symbol_map, simultaneous=True))
    return result


def relabelling_symmetry():
    """{(n, r): True} if the vertex is invariant under every exchange of two adjacent legs of the
    same kind (these exchanges generate all relabellings within each kind)."""
    results = {}
    for n, radial in PATTERNS:
        value = closed_vertex(n, radial, eliminate_lambda=True)
        ok = True
        for group in (range(1, radial + 1), range(radial + 1, n + 1)):
            group = list(group)
            for first, second in zip(group, group[1:]):
                permutation = {j: j for j in range(1, n + 1)}
                permutation[first], permutation[second] = second, first
                ok &= (relabel(value, permutation) - value).is_zero()
        results[(n, radial)] = ok
    return results


# ---------------------------------------------------------------------------
# Dimensions
# ---------------------------------------------------------------------------

def dimension_check():
    """V^(n) -> c^{4-n} V^(n) under X -> c^{d_X} X: d = 1 for B_0, f_pi, m_sigma and each source
    letter; d = 2 for k_i.k_j; d = 0 for eps, lambda_x, xi_x (lambda eliminated through m_sigma^2)."""
    scale = sp.Symbol('c', positive=True)
    results = {}
    for n, radial in PATTERNS:
        value = closed_vertex(n, radial, eliminate_lambda=True)
        ok = True
        for key, coefficient in value.terms.items():
            dots = {symbol: scale ** 2 * symbol for symbol in sp.sympify(coefficient).free_symbols
                    if symbol.name.startswith('d_')}
            scaled = sp.sympify(coefficient).subs({B0: scale * B0, F_PI: scale * F_PI, M_SIGMA: scale * M_SIGMA,
                                                   **dots}, simultaneous=True)
            scaled *= scale ** cartesian.source_count(key)
            ok &= sp.simplify(scaled - scale ** (4 - n) * coefficient) == 0
        results[(n, radial)] = ok
    return results


# ---------------------------------------------------------------------------
# Momentum dependence
# ---------------------------------------------------------------------------

def momentum_dependent_patterns():
    """Patterns whose vertex depends on momenta (expected: only n = 4 with no radial leg,
    through the four-meson kinetic vertex)."""
    found = []
    for n, radial in PATTERNS:
        value = closed_vertex(n, radial)
        if any(symbol.name.startswith('d_') for coefficient in value.terms.values()
               for symbol in sp.sympify(coefficient).free_symbols):
            found.append((n, radial))
    return found


def radial_potential_has_no_momentum():
    """The quartic (and cubic) V_phi vertices from the direct expansion with lambda only."""
    for n, radial in PATTERNS:
        value = cartesian.direct_vertex(n, radial).filtered(lambda key: key == ())
        for coefficient in value.terms.values():
            if any(symbol.name.startswith('d_') for symbol in sp.sympify(coefficient).free_symbols):
                return False
    return True


def sp_kinetic_term_gives_no_vertex():
    """1/2 (d s)^2 + 1/2 (d p)^2 expanded with formal legs: no term of order >= 3 (direct check)."""
    for n, radial in PATTERNS:
        fields = sp.symbols(f'phi1:{n + 1}')
        markers = sp.symbols(f'kappa1:{n + 1}')
        xis = sp.symbols(f'xi1:{n + 1}', real=True)
        ds = sum((-sp.I * markers[j] * fields[j] for j in range(radial)), sp.Integer(0))
        dp = sum((-sp.I * markers[j] * xis[j] * fields[j] for j in range(radial, n)), sp.Integer(0))
        expression = sp.expand(sp.Rational(1, 2) * (ds ** 2 + dp ** 2))
        for term in sp.Add.make_args(expression):
            powers = term.as_powers_dict()
            if all(powers.get(one_field, 0) == 1 for one_field in fields):
                return False
    return True


# ---------------------------------------------------------------------------
# Up-only limit
# ---------------------------------------------------------------------------

def up_only_checks():
    """{(n, r): closed rule after D-A4 substitution == direct expansion with M_hat = m_u P_u phi^*/f_a}."""
    results = {}
    for n, radial in PATTERNS:
        closed = cartesian.up_only_substitution(closed_vertex(n, radial))
        direct = cartesian.up_only_substitution(cartesian.direct_vertex(n, radial, up_only=True))
        results[(n, radial)] = (closed - direct).is_zero()
    return results


def up_only_source_patterns():
    """In the up-only limit only (r, q) in {(0,0), (1,0), (0,1)} survive, so a source vertex has at
    most one radial or ALP-field end.  Returns the set of (r, q) of the surviving source letters."""
    found = set()
    for n, radial in PATTERNS:
        value = cartesian.source_vertex(cartesian.generic_legs(n, radial))
        for key in value.terms:
            for word in key:
                for name in word:
                    if cartesian.is_cartesian_source(name) and cartesian.up_only_letter(name) is not None:
                        _, r, q = cartesian.source_indices(name)
                        found.add((r, q))
    return found


# ---------------------------------------------------------------------------
# Closed rules
# ---------------------------------------------------------------------------

def radial_potential_table():
    """{(r, q): vertex} of V_phi with lambda kept, legs: r radial, q bare p legs."""
    table = {}
    for n in (3, 4, 5):
        for radial in range(n + 1):
            legs = ([cartesian.radial_leg({f'q{j}': 1}) for j in range(radial)]
                    + [cartesian.alp_field_leg({f'q{j}': 1}) for j in range(radial, n)])
            value = cartesian.radial_potential_vertex(legs, eliminate_lambda=False).simplified()
            if value.terms:
                table[(radial, n - radial)] = value.terms[()]
    return table


def closed_rules_latex():
    return [
        r'V^{(n)}_{\rm src}=B_0f_\pi^{2-n}\sum_{(R,A,P)}r!\,s!\,\epsilon^{r+s}\frac{(-1)^{m(m-1)/2}}{m!}'
        r'\prod_{x\in A}\xi_x\sum_{\pi\in S_m}\langle X^C_{r,s}\lambda_{\pi_1}\cdots\lambda_{\pi_m}\rangle,'
        r'\quad X^C=\widehat S^C\ (m\ \text{even}),\ \widehat P^C\ (m\ \text{odd})',
        r'V^{(4)}_{\rm kin}=\frac{1}{48f_\pi^2}\sum_{\pi\in S_4}(-k_{\pi_2}\cdot k_{\pi_4})'
        r'\langle[\lambda_{\pi_1},\lambda_{\pi_2}][\lambda_{\pi_3},\lambda_{\pi_4}]\rangle',
        r'V_\varphi:\ V(s,s,s)=-6\lambda f_a,\ V(s;x,y)=-2\lambda f_a\,\xi_x\xi_y,\ '
        r'V(s^4)=-6\lambda,\ V(s,s;x,y)=-2\lambda\,\xi_x\xi_y,\ V(w,x,y,z)=-6\lambda\,\xi_w\xi_x\xi_y\xi_z',
        r'\lambda f_a=\frac{\epsilon}{2f_\pi}\big[m_\sigma^2-\epsilon^2B_0(\langle\widehat S^C_{1,0}\rangle'
        r'-2\langle\widehat S^C_{2,0}\rangle)\big],\qquad \lambda=\frac{\epsilon}{f_\pi}\,\lambda f_a',
    ]


def main():
    print('Closed Cartesian vertex rules (GeV; all momenta incoming):')
    for line in closed_rules_latex():
        print('  ' + line)
    print('\nV_phi vertices (r radial, q ALP-field legs; lambda kept):')
    for key, value in radial_potential_table().items():
        print(f'  {key}: {sp.latex(value)}')
    print('\nConstruction 1 (direct expansion) against construction 2 (closed rules):')
    for key, (direct_ok, tree_ok) in compare_constructions().items():
        print(f'  n={key[0]}, radial legs={key[1]}: direct = rule: {direct_ok}; trees = words: {tree_ok}')
    print('\nRelabelling symmetry:', all(relabelling_symmetry().values()))
    print('Dimension check c^(4-n):', all(dimension_check().values()))
    print('Momentum-dependent patterns:', momentum_dependent_patterns())
    print('V_phi vertices without momentum dependence:', radial_potential_has_no_momentum())
    print('No vertex from 1/2 (d s)^2 + 1/2 (d p)^2:', sp_kinetic_term_gives_no_vertex())
    print('Up-only limit, closed = direct (M_hat = m_u P_u phi^*/f_a):', all(up_only_checks().values()))
    print('Up-only surviving source indices (r, q):', sorted(up_only_source_patterns()))
    example = cartesian.up_only_substitution(closed_vertex(3, 1))
    print('Up-only V(s; x, y):', cartesian.trace_sum_latex(example))


if __name__ == '__main__':
    main()
