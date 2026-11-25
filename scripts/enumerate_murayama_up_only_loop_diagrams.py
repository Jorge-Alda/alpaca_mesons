#!/usr/bin/env python3
"""Enumerate the one-loop cycle diagrams of the loop-size estimate for a -> pi+ pi- pi0 (up-only model).

NO LOOP IS COMPUTED and no vertex is evaluated.  The script only works with diagram topology:
it lists the ways of attaching the external objects to the vertices of a one-loop cycle with given
internal lines, removes the diagrams forbidden by selection rules, and compares the result with the
hand-written list DIAGRAMS of scripts/evaluate_explicit_estimate.py (used by
scripts/estimate_loop_size_dimensionless.py).

Objects attached to the loop (names as in evaluate_explicit_estimate.py):
    'A' external ALP-like pole leg, '+' pi+, '-' pi-, '0' pi0-like pole leg (pseudoscalars);
    'beta'  radial branch sigma -> (pi+, pi-) entering the loop (a scalar radial leg);
    'gamma' radial branch sigma -> (A, pi0) entering the loop (a scalar radial leg).
Admissible sets of objects (each external leg appears once):
    {A, +, -, 0},  {beta, A, 0},  {gamma, +, -}   (vertex corrections and boxes, OBJECT_SETS);
    {beta, gamma}  (self-energy of the internal sigma line of the tree exchange,
                    SIGMA_SELF_ENERGY_OBJECTS).  Only its cycle with two S lines (two sigma^3
                    vertices, diagram 6.7 of class 6) is in the estimate; its cycles with A or M lines (the
                    sigma self-energy from a pseudoscalar loop) are listed in
                    EXCLUDED_SIGMA_SELF_ENERGY_CYCLES and are not part of the decided estimate.
Internal lines: 'A' (ALP-like pole), 'M' (meson, nine U(3) generators, may carry charge),
'S' (radial sigma).  Line j joins vertex j and vertex j+1 of the cycle (cyclically).

Rules (assumptions of the decided estimate, stated explicitly):
    R0 every vertex carries at least one object (no loop-line self-energy insertions, and no
       sigma tadpoles: the sigma one-point function is removed by stationarity);
    R1 intrinsic parity: at each vertex the number of pseudoscalar fields
       (external A, +, -, 0 and internal A, M lines) is even;
    R2 radial content of the up-only Lagrangian used by evaluate_explicit_estimate.py:
       a vertex with one radial field (beta, gamma or an S line) is allowed; a vertex with two
       radial fields exists only as the 4-field kinetic vertex sigma sigma a a, so it needs exactly
       two further fields that carry an ALP component (external A or 0, internal A or M lines;
       pi+ and pi- have none); a vertex with three radial fields exists only as the sigma^3 vertex
       of the radial potential (V_sigma3 = -6 lambda f_a), so it has no further field; four or
       more radial fields: no vertex (no S_hat_2s; the sigma^4 term of the radial potential is
       not part of the vertex set);
    R3 electric charge: A, S lines and A, 0, beta, gamma are neutral, M lines carry -1, 0 or +1;
       the cycle must admit a charge assignment of the M lines conserved at every vertex;
    R4 three-or-more-line truncation of the estimate (unchanged): cycles with three or more A or S
       lines are not included (OTHER_CYCLES lists the ones that survive R0-R3 and are dropped by
       this truncation).
Diagram classes (labels of evaluate_explicit_estimate.py, diagram label "class.index"): class 1 one
A line at one vertex, class 2 A and M lines, class 3 two A lines, class 4 one S line and M lines,
class 5 one S, one A and M lines, class 6 two S lines and M lines (including the sigma self-energy
6.7); classes 7 and 8 are external-leg residues and are not cycles of this enumeration.
Symmetry factor S_d = 1/|Aut(d)|, Aut(d) = dihedral symmetries of the cycle that preserve the
line labels and the blocks (1/2 for a closed line, 1/2 for two identical lines between two
different blocks, 1 otherwise).

Run:
    .venv/bin/python scripts/enumerate_murayama_up_only_loop_diagrams.py
Tests:
    PYTHONPATH=scripts .venv/bin/python -m unittest test_enumerate_murayama_up_only_loop_diagrams
"""
import itertools
from collections import Counter

PSEUDOSCALAR_OBJECTS = {'A', '+', '-', '0'}
RADIAL_OBJECTS = {'beta', 'gamma'}
ALP_CARRYING_OBJECTS = {'A', '0'}
OBJECT_CHARGE = {'A': 0, '+': 1, '-': -1, '0': 0, 'beta': 0, 'gamma': 0}
OBJECT_SETS = (('A', '+', '-', '0'), ('beta', 'A', '0'), ('gamma', '+', '-'))
SIGMA_SELF_ENERGY_OBJECTS = (('beta', 'gamma'),)
PSEUDOSCALAR_LINES = {'A', 'M'}
ALP_CARRYING_LINES = {'A', 'M'}

# Group of diagrams of the estimate -> (class, multiset of internal lines, admissible object sets).
CYCLE_GROUPS = {
    'class 1': ('1', ('A',), OBJECT_SETS),
    'class 2': ('2', ('A', 'M'), OBJECT_SETS),
    'class 3': ('3', ('A', 'A'), OBJECT_SETS),
    'class 4 bubble': ('4', ('S', 'M'), OBJECT_SETS),
    'class 4 triangle': ('4', ('S', 'M', 'M'), OBJECT_SETS),
    'class 5 bubble': ('5', ('S', 'A'), OBJECT_SETS),
    'class 5 triangle': ('5', ('S', 'A', 'M'), OBJECT_SETS),
    'class 6 bubble': ('6', ('S', 'S'), OBJECT_SETS),
    'class 6 triangle': ('6', ('S', 'S', 'M'), OBJECT_SETS),
    'class 6 box': ('6', ('S', 'S', 'M', 'M'), OBJECT_SETS),
    'class 6 sigma self-energy': ('6', ('S', 'S'), SIGMA_SELF_ENERGY_OBJECTS),
}
# Further cycles shown to be empty under R0-R3, or dropped by the three-or-more-line truncation R4.
OTHER_CYCLES = {
    'A triangle': ('A', 'M', 'M'),
    'A box': ('A', 'M', 'M', 'M'),
    'two A lines, triangle': ('A', 'A', 'M'),
    'two A lines, box': ('A', 'A', 'M', 'M'),
    'class 4 box': ('S', 'M', 'M', 'M'),
    'class 5 box': ('S', 'A', 'M', 'M'),
    'closed S line (rule R2)': ('S',),
    'three S lines (three-or-more-line truncation)': ('S', 'S', 'S'),
    'two S lines and an A line, triangle (three-or-more-line truncation)': ('S', 'S', 'A'),
    'two S lines and an A line, box (three-or-more-line truncation)': ('S', 'S', 'A', 'M'),
    'two S lines and two A lines, box (three-or-more-line truncation)': ('S', 'S', 'A', 'A'),
    'two A lines and an S line, triangle (three-or-more-line truncation)': ('S', 'A', 'A'),
    'two A lines and an S line, box (three-or-more-line truncation)': ('S', 'A', 'A', 'M'),
}
# Cycles of the sigma self-energy {beta, gamma} with pseudoscalar lines: not in the estimate.
EXCLUDED_SIGMA_SELF_ENERGY_CYCLES = {
    'sigma self-energy, two A lines': ('A', 'A'),
    'sigma self-energy, A and M lines': ('A', 'M'),
    'sigma self-energy, two M lines': ('M', 'M'),
}


def group_of_code_diagram(label, blocks, lines):
    """Name of the entry of CYCLE_GROUPS that contains a diagram of the code list: same class, same
    multiset of internal lines and an admissible object set containing its objects."""
    from evaluate_explicit_estimate import diagram_class
    objects = sorted(o for block in blocks for o in block)
    matches = [name for name, (class_name, line_multiset, object_sets) in CYCLE_GROUPS.items()
               if class_name == diagram_class(label) and sorted(line_multiset) == sorted(lines)
               and any(sorted(candidate) == objects for candidate in object_sets)]
    assert len(matches) == 1, (label, matches)
    return matches[0]


def dihedral_maps(size):
    """(vertex map, line map) for every dihedral symmetry of a cycle with `size` vertices."""
    maps = []
    for shift in range(size):
        maps.append(([(i + shift) % size for i in range(size)], [(j + shift) % size for j in range(size)]))
        maps.append(([(shift - i) % size for i in range(size)], [(shift - j - 1) % size for j in range(size)]))
    return maps


def transform(diagram, vertex_map, line_map):
    blocks, lines = diagram
    size = len(blocks)
    new_blocks, new_lines = [None] * size, [None] * size
    for i in range(size):
        new_blocks[vertex_map[i]] = blocks[i]
        new_lines[line_map[i]] = lines[i]
    return tuple(new_blocks), tuple(new_lines)


def canonical(diagram):
    """Representative of the dihedral class of (blocks, lines); blocks are sorted tuples."""
    blocks, lines = diagram
    diagram = (tuple(tuple(sorted(block)) for block in blocks), tuple(lines))
    return min(transform(diagram, *maps) for maps in dihedral_maps(len(lines)))


def symmetry_factor(diagram):
    blocks, lines = diagram
    diagram = (tuple(tuple(sorted(block)) for block in blocks), tuple(lines))
    automorphisms = sum(transform(diagram, *maps) == diagram for maps in dihedral_maps(len(lines)))
    return 1 / automorphisms


def lines_at_vertex(lines, vertex):
    """Lines touching a vertex (a closed line counts twice)."""
    size = len(lines)
    return [lines[(vertex - 1) % size], lines[vertex]]


def parity_allowed(blocks, lines):
    return all((sum(o in PSEUDOSCALAR_OBJECTS for o in block)
                + sum(line in PSEUDOSCALAR_LINES for line in lines_at_vertex(lines, v))) % 2 == 0
               for v, block in enumerate(blocks))


def radial_allowed(blocks, lines):
    for v, block in enumerate(blocks):
        touching = lines_at_vertex(lines, v)
        radial = sum(o in RADIAL_OBJECTS for o in block) + touching.count('S')
        others = len(block) + len(touching) - radial
        if radial >= 4 or (radial == 3 and others != 0):
            return False
        if radial == 2:
            alp_carrying = (sum(o in ALP_CARRYING_OBJECTS for o in block)
                            + sum(line in ALP_CARRYING_LINES for line in touching))
            if others != 2 or alp_carrying != 2:
                return False
    return True


def charge_allowed(blocks, lines):
    size = len(lines)
    meson_lines = [j for j, line in enumerate(lines) if line == 'M']
    for charges in itertools.product((-1, 0, 1), repeat=len(meson_lines)):
        line_charge = [0] * size
        for j, q in zip(meson_lines, charges):
            line_charge[j] = q      # charge flowing from vertex j to vertex j+1
        if all(sum(OBJECT_CHARGE[o] for o in block) + line_charge[(v - 1) % size] - line_charge[v] == 0
               for v, block in enumerate(blocks)):
            return True
    return False


def enumerate_diagrams(line_multiset, rules=('parity', 'radial', 'charge'), object_sets=OBJECT_SETS):
    """Canonical diagrams with the given internal lines and object sets that pass the selected
    rules (rule R0 always; the three-or-more-line truncation R4 is applied by the choice of line
    multisets)."""
    checks = {'parity': parity_allowed, 'radial': radial_allowed, 'charge': charge_allowed}
    size = len(line_multiset)
    found = set()
    for lines in set(itertools.permutations(line_multiset)):
        for objects in object_sets:
            for assignment in itertools.product(range(size), repeat=len(objects)):
                if len(set(assignment)) != size:
                    continue
                blocks = tuple(tuple(o for o, v in zip(objects, assignment) if v == vertex) for vertex in range(size))
                if all(checks[rule](blocks, lines) for rule in rules):
                    found.add(canonical((blocks, lines)))
    return found


def code_diagrams():
    """{group: [(canonical diagram, symmetry factor, diagram label)]} for evaluate_explicit_estimate.DIAGRAMS."""
    from evaluate_explicit_estimate import DIAGRAMS
    groups = {}
    for label, weight, blocks, lines in DIAGRAMS:
        groups.setdefault(group_of_code_diagram(label, blocks, lines), []).append(
            (canonical((blocks, lines)), weight, label))
    return groups


LEG_NAMES = {'A': 'a', '+': 'pi+', '-': 'pi-', '0': 'pi0', 'beta': 'beta', 'gamma': 'gamma'}
ORDER = ['A', '+', '-', '0', 'beta', 'gamma']


def written(diagram):
    """Plain-text form V({..}) -L1- V({..}) -L2- ... (last line closes the cycle)."""
    blocks, lines = diagram
    parts = []
    for block, line in zip(blocks, lines):
        names = ','.join(LEG_NAMES[o] for o in sorted(block, key=ORDER.index))
        parts.append(f'{{{names}}} -{line}-')
    return ' '.join(parts) + ' (back to first)'


def main():
    code = code_diagrams()
    for group, (_, line_multiset, object_sets) in CYCLE_GROUPS.items():
        allowed = enumerate_diagrams(line_multiset, object_sets=object_sets)
        listed = code.get(group, [])
        listed_set = {d for d, _, _ in listed}
        print(f'== {group}: lines {line_multiset}, r = {len(line_multiset)}')
        print(f'   code list {len(listed)} (distinct {len(listed_set)}), allowed by all rules {len(allowed)}')
        for diagram, weight, label in listed:
            status = 'allowed' if diagram in allowed else 'FORBIDDEN, must not be listed'
            assert abs(weight - symmetry_factor(diagram)) < 1e-12
            print(f'   {label:5s} S = {weight:.1f}  {written(diagram)}  [{status}]')
        for diagram in sorted(allowed - listed_set):
            print(f'   MISSING from code list: {written(diagram)}')
    for name, line_multiset in OTHER_CYCLES.items():
        found = enumerate_diagrams(line_multiset)
        print(f'== other cycle {name}: {len(found)} diagrams allowed by R0-R3')
        for diagram in sorted(found):
            print(f'   dropped by the three-or-more-line truncation: {written(diagram)}')
    for name, line_multiset in EXCLUDED_SIGMA_SELF_ENERGY_CYCLES.items():
        found = enumerate_diagrams(line_multiset, object_sets=SIGMA_SELF_ENERGY_OBJECTS)
        for diagram in sorted(found):
            print(f'== excluded {name}: {written(diagram)}, S = {symmetry_factor(diagram):.1f}')


if __name__ == '__main__':
    main()
