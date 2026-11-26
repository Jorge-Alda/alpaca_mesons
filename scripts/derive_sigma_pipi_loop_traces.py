"""Compact flavour traces of the kept one-loop graphs of sigma -> pi+ pi-, with per-term loop scales and the termwise routing maximum (tasks 3.1, 3.2, 4.4).

Run:
    PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_loop_traces.py
    PYTHONPATH=scripts .venv/bin/python scripts/derive_sigma_pipi_loop_traces.py --markdown
    PYTHONPATH=scripts .venv/bin/python -m unittest test_derive_sigma_pipi_loop_traces

Units GeV; conventions of sigma_pipi_loop_algebra.py.  The kept graphs are
those of derive_sigma_pipi_loop_graphs.py (26 one-particle-irreducible graphs,
6 sigma and 5 + 5 pi+- self-energies).  For the four meson-only graphs with an
excluded piece (decision D4), the term in which every meson line takes the
bare U(3) generator sum is removed exactly in every routing; everything else
is kept.

For every kept graph and every routing the product of vertices is built with
the generator sums of the meson lines removed by U(3) completeness,
    sum_i <A lambda_i B lambda_i> = 2 <A><B>,  sum_i <A lambda_i><lambda_i B> = 2 <A B>,
in two independent ways:
  construction 1: vertices expanded into words, completeness on words;
  construction 2: vertices as compact trees (sums over permutations,
                  commutators, anticommutators), completeness applied to the
                  trees before any expansion.
The two results are compared as exact identities of cyclic words, routing by
routing.

Loop scale (decision D3): every term (one piece per line) has its own scale
m_term, the largest of the scales of its lines: m_sigma (radial line), m_a
(ALP-like line, and the ALP-field entry, mixed entries and meson block of the
ALP-like pole of a meson line), m_P (bare generator sum of a meson line, m_P
the largest pole mass of the nine meson-like states, symbolic).  Products k.k
and k.p are set to m_term^2 and the size factor is
S_Gamma (m_term^2)^{2-r_Gamma}/(16 pi^2); for a self-energy the residue
correction dZ = -Pi'(m^2) is taken term by term at m_term.

Routing (decision D12): for every routing rho (loop momentum on each internal
line, both orientations) the size-rule value is sum_t C_t^rho t over the
independent products of traces t of words; the graph enters with
sum_t C_t^{rho*_t} t, rho*_t the routing of largest |C_t^rho| (sign kept).
Ties of equal magnitude and opposite sign take the sign of the sum of the
coefficients over all routings; if this sum vanishes, the first routing in the
canonical order is used and the product is marked (a conventional choice at a
degenerate mass relation).  Where rho*_t depends on the masses the coefficient
is a Piecewise expression in homogeneous inequalities between squared masses.
This is an upper bound per graph on the routing ambiguity of the size rule,
not on the loop.  The range of the ratios C_t^rho/C_t^{rho*} is reported as
the routing spread.  No loop integral, loop function or loop amplitude is
evaluated.
"""
import sys

import sympy as sp

import sigma_pipi_loop_algebra as algebra
import derive_sigma_pipi_loop_graphs as graphs

GROUPS = ('one_particle_irreducible', 'sigma', 'pi+', 'pi-')


def setups():
    return {'one_particle_irreducible': graphs.three_point_setup(),
            'sigma': graphs.two_point_setup('sigma'), 'pi+': graphs.two_point_setup('pi+'),
            'pi-': graphs.two_point_setup('pi-')}


def routing_trees(row, legs, basis_dot, group):
    """Size-rule value of every routing in construction 2 (compact trees)."""
    return graphs.routing_values(row, legs, basis_dot, group, 'trees', pieces=row['pieces'])


def spread_summary(maximum):
    """Diagnostic of the routing spread from graphs.termwise_maximum."""
    ratios = [ratio for values in maximum['spread'].values() for ratio in values]
    return {'products': sum(maximum['counts'].values()), 'numerical': maximum['counts']['numerical'],
            'mass-dependent': maximum['counts']['mass-dependent'],
            'minimum ratio': min(ratios) if ratios else None,
            'zero routing': any(ratio == 0 for ratio in ratios),
            'marked': len(maximum['marked']), 'marked boundary': len(maximum['marked boundary'])}


def kept_graph_traces():
    """List of dicts with the kept graphs: termwise routing maximum (words),
    its compact form where a single routing attains it, the agreement of the
    two constructions in every routing and the routing spread.  For 1PI graphs
    the value is S_Gamma (m_term^2)^{2-r}/(16 pi^2) T per term; for self-energies dZ."""
    data = graphs.summary()
    entries = []
    for group, (legs, basis_dot) in setups().items():
        for row in graphs.kept(data[group]):
            trees_per_routing = routing_trees(row, legs, basis_dot, group)
            agree = all((trees.to_trace_sum() - words).is_zero()
                        for trees, words in zip(trees_per_routing, row['routing_values']))
            common = graphs.common_routing(row['routing_values'], row['value'])
            base = common if common is not None else graphs.base_routing(row['routing_values'], row['value'])
            correction = (row['value'] - row['routing_values'][base]).simplified()
            correction.terms = {key: value for key, value in correction.terms.items()
                                if not graphs.expression_is_zero(value)}
            entry = {**row, 'group': group, 'words': row['value'], 'constructions_agree': agree,
                     'common routing': common, 'base routing': base, 'trees': trees_per_routing[base],
                     'correction': correction, 'spread': spread_summary(row['maximum'])}
            if group != 'one_particle_irreducible':
                entry['residue'] = row['value']
            entries.append(entry)
    return entries


def sources_in(entries):
    """Dressed source letters that appear in the kept T_Gamma."""
    letters = set()
    for entry in entries:
        letters |= {name for name in entry['words'].letters() if algebra.is_source(name)}
    return sorted(letters)


def no_component_notation(text):
    """Output guard: no matrix indices, no E_ij, no generator symbols."""
    forbidden = ('E_', '_{ij}', 'lambda_i', 'lambda_j', '\\lambda_{i', 'G1', 'G2', 'G3', "G'")
    return not any(token in text for token in forbidden)


# ---------------------------------------------------------------------------
# Text output
# ---------------------------------------------------------------------------
LINE_NAMES = {'S': 'radial', 'A': 'ALP-like', 'M': 'meson'}
LEG_NAMES = {'sigma': r'\sigma', 'pi+': r'\pi^+', 'pi-': r'\pi^-'}


def scales_latex(row):
    return ','.join(algebra.coefficient_latex(scale) for scale in graphs.possible_scales(row))


def lines_latex(row):
    names = [LINE_NAMES[line] for line in row['lines']]
    if len(names) == 1:
        return names[0] + ' (closed on the vertex)'
    if len(names) == 3:
        blocks = row['blocks']
        texts = ['$' + LEG_NAMES[blocks[i][0]] + '$' for i in range(3)]
        return '; '.join(f'{texts[i]}–{texts[(i + 1) % 3]}: {names[i]}' for i in range(3))
    return ', '.join(names)


def blocks_latex(row):
    return ', '.join('$(' + ','.join(LEG_NAMES[name] for name in block) + ')$' for block in row['blocks'])


def spread_text(entry):
    spread = entry['spread']
    if not spread['numerical'] and not spread['mass-dependent']:
        return 'routing independent'
    parts = [f'routing dependent in {spread["numerical"] + spread["mass-dependent"]} of {spread["products"]} '
             f'products of traces of words']
    if spread['numerical']:
        parts.append(f'{spread["numerical"]} with numerical ratios: '
                     f'$C_t^\\rho/C_t^{{\\rho^*}}$ from ${sp.latex(spread["minimum ratio"])}$ to $1$')
    if spread['mass-dependent']:
        parts.append(f'{spread["mass-dependent"]} with a mass-dependent largest routing (Piecewise)')
    if spread['marked']:
        parts.append(f'{spread["marked"]} decided by the first routing (sum over routings zero)')
    if spread['marked boundary']:
        parts.append(f'{spread["marked boundary"]} decided by the first routing at a boundary point '
                     f'(sum over routings zero)')
    if entry['common routing'] is not None:
        parts.append(f'attained by routing {entry["common routing"] + 1} for every product')
    return '; '.join(parts)


def display_coefficient(value):
    """A Piecewise coefficient written as its common factor times a Piecewise of
    the remaining factors (display only; the value is unchanged)."""
    value = sp.piecewise_fold(sp.sympify(value))
    if not isinstance(value, sp.Piecewise):
        return sp.factor(value)
    common = sp.gcd_list([piece for piece, _ in value.args])
    return sp.Mul(sp.factor(common), sp.Piecewise(*[(sp.factor(sp.cancel(piece / common)), condition)
                                                    for piece, condition in value.args]), evaluate=False)


def paired_compact(trace_sum):
    """CompactSum of a TraceSum of words in which two words W l+ l- and W l- l+
    with equal coefficients are written as <W {l+, l-}>."""
    remaining = dict(trace_sum.simplified().terms)
    result = algebra.CompactSum()
    for key in sorted(remaining, key=repr):
        if key not in remaining:
            continue
        value = remaining.pop(key)
        trees = None
        for index, word in enumerate(key):
            if word.count('l+') != 1 or word.count('l-') != 1:
                continue
            position = word.index('l+')
            rotated = word[position:] + word[:position]
            if rotated[1] != 'l-':
                continue
            rest = rotated[2:]
            others = key[:index] + key[index + 1:]
            partner = tuple(sorted(others + (algebra.canonical_word(('l-', 'l+') + rest),)))
            if partner in remaining and graphs.expression_is_zero(remaining[partner] - value):
                remaining.pop(partner)
                anticommutator = ('A', algebra.letter('l+'), algebra.letter('l-'))
                trees = tuple(('P', *(algebra.letter(name) for name in other)) for other in others)
                trees += (('P', *(algebra.letter(name) for name in rest), anticommutator),)
                break
        if trees is None:
            trees = tuple(('P', *(algebra.letter(name) for name in word)) for word in key)
        result.add_term(trees, value)
    result.terms = {key: display_coefficient(value) for key, value in result.terms.items()}
    return result


def displayed(entry):
    """Compact form of the maximum: the trees of the base routing (the one that
    attains the maximum for most products, or for all of them) plus, where
    needed, the correction sum_t (C_t^{rho*_t} - C_t^{base}) t in words with
    anticommutators.  Returns (base trees, correction CompactSum)."""
    return entry['trees'], paired_compact(entry['correction'])


def aligned(compact_sum, name):
    body = algebra.compact_latex(compact_sum, separator=' \\\\\n&+')
    return '$$\n\\begin{aligned}\n' + name + '={}&' + body + '\n\\end{aligned}\n$$'


def markdown_blocks(entries):
    """Result blocks of the two trace records (one per kept graph)."""
    blocks = {'one_particle_irreducible': [], 'self_energies': []}
    number = {'one_particle_irreducible': 0, 'sigma': 0, 'pi+': 0, 'pi-': 0}
    prefix = {'sigma': 'σ', 'pi+': 'π', 'pi-': 'π'}
    for entry in entries:
        group = entry['group']
        if group == 'pi-':
            continue
        number[group] += 1
        if group == 'one_particle_irreducible':
            opening = f'**Graph {number[group]}.** Vertices {blocks_latex(entry)}; lines'
            name = 'X_\\Gamma^{\\max}'
        else:
            opening = f'**Graph {prefix[group]}{number[group]}** (self-energy of ${LEG_NAMES[group]}$). Lines'
            name = '\\delta Z_\\Gamma^{\\max}'
        text = (f'{opening} {lines_latex(entry)}; '
                f'$S_\\Gamma={sp.latex(entry["symmetry"])}$, $r_\\Gamma={len(entry["blocks"])}$, '
                f'$m_{{\\rm term}}\\in\\{{{scales_latex(entry)}\\}}$'
                + ('; D4 term removed' if entry['excluded'] is not None else '')
                + f'; {spread_text(entry)}.\n')
        base, correction = displayed(entry)
        routing_name = f'{name}\\big|_{{\\rho={entry["base routing"] + 1}}}' if correction.terms else name
        text += aligned(base, routing_name)
        if correction.terms:
            text += ('\nThe routing $\\rho=' + str(entry['base routing'] + 1) + '$ attains the maximum for '
                     'most products; the other products are corrected by\n'
                     + aligned(correction, name + '-' + routing_name))
        key = 'one_particle_irreducible' if group == 'one_particle_irreducible' else 'self_energies'
        blocks[key].append(text)
    return blocks


def main():
    entries = kept_graph_traces()
    if '--markdown' in sys.argv:
        blocks = markdown_blocks(entries)
        for key in ('one_particle_irreducible', 'self_energies'):
            print(f'<!-- {key} -->')
            print('\n\n'.join(blocks[key]))
        return
    print('Kept graphs: termwise routing maximum of the size-rule value (1PI: S (m_term^2)^{2-r} T/(16 pi^2); '
          'self-energies: dZ).')
    for entry in entries:
        spread = entry['spread']
        print(f'\n[{entry["group"]}] {entry["label"]}  S={entry["symmetry"]}  r={len(entry["blocks"])}  '
              f'm_term in {graphs.possible_scales(entry)}  '
              f'D4 term removed: {entry["excluded"] is not None}')
        print(f'  constructions agree in every routing: {entry["constructions_agree"]}; '
              f'trace products: {len(entry["words"].terms)} (words); common routing: {entry["common routing"]}')
        print(f'  routing: {spread}')
        base, correction = displayed(entry)
        text = algebra.compact_latex(base, separator='\n  + ')
        assert no_component_notation(text)
        print(f'  routing {entry["base routing"] + 1}: ' + text)
        if correction.terms:
            text = algebra.compact_latex(correction, separator='\n  + ')
            assert no_component_notation(text)
            print('  correction to the maximum: ' + text)
    print('\nDressed sources in the kept values:', ', '.join(sources_in(entries)))


if __name__ == '__main__':
    main()
