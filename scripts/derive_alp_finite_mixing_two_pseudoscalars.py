#!/usr/bin/env python3
"""Exact three-leg pole projection; natural units, all momenta incoming.

Run: .venv/bin/python scripts/derive_alp_finite_mixing_two_pseudoscalars.py
The vertex tensor is the differentiated CCWZ Lagrangian and excludes the
Feynman-rule factor i. No truncation of neutral pole vectors is made.
"""
from itertools import product
import sympy as sp
from derive_alp_charged_u3_vector_redefinition import derive as charged_check
from derive_alp_pole_projected_amplitudes import ccwz_vertex


def projected_vertex(vertex, parent, first, second):
    """Contract an already symmetrized vertex with three pole directions."""
    return sp.expand(sum(vertex[i, j, k]*parent[i]*first[j]*second[k]
                         for i, j, k in product(range(len(parent)), repeat=3)))


def projected_ccwz_vertex(meson_basis, sources, parent, first, second, gram):
    """Evaluate the CCWZ cubic vertex on three supplied pole directions."""
    momenta = [sp.eye(3)[:, i] for i in range(3)]
    return ccwz_vertex(meson_basis, sources, [parent, first, second], momenta, gram)


def derive():
    # Two neutral meson directions and one ALP suffice to display every
    # component class; the contraction above accepts any number of fields.
    fields = sp.symbols('P_1 P_2 a')
    coefficients = {(i,j,k): sp.Symbol('V_' + ''.join(map(str,(i,j,k))))
                    for i in range(3) for j in range(i,3) for k in range(j,3)}
    lagrangian = sum(coefficients[i,j,k]*fields[i]*fields[j]*fields[k]
                     / (sp.factorial((i,j,k).count(i)) if i == k else
                        2 if i == j or j == k else 1)
                     for i,j,k in coefficients)
    tensor = sp.MutableDenseNDimArray.zeros(3,3,3)
    for i,j,k in product(range(3),repeat=3):
        tensor[i,j,k] = sp.diff(lagrangian,fields[i],fields[j],fields[k])
    parent = sp.Matrix(sp.symbols('v_a1 v_a2 v_aa'))
    first = sp.Matrix(sp.symbols('v_11 v_12 v_1a'))
    second = sp.Matrix(sp.symbols('v_21 v_22 v_2a'))
    amplitude = projected_vertex(tensor,parent,first,second)
    mixing = sp.Symbol('epsilon', real=True)
    # A nondegenerate, ALP-like parent and meson-like daughters; the
    # coefficients scale with one power per ALP source insertion.
    substitutions = {parent[0]:mixing*parent[0], parent[1]:mixing*parent[1],
                     first[2]:mixing*first[2], second[2]:mixing*second[2]}
    for (i,j,k), coefficient in coefficients.items():
        substitutions[coefficient] = (coefficient*mixing**sum(index == 2 for index in (i,j,k)) if 2 in (i,j,k) else 0)
    weak = sp.series(amplitude.subs(substitutions),mixing,0,2).removeO().expand()
    return tensor,parent,first,second,amplitude,weak


def main():
    tensor,parent,first,second,amplitude,weak = derive()
    print('Exact pole-projected three-point tensor contraction:')
    print(sp.latex(amplitude))
    print('Leading term for nondegenerate weak mixing:')
    print(sp.latex(weak))
    contact,charged,on_shell_contact = charged_check()
    print('Charged U(3) rotation contact off shell and on shell:')
    print(sp.latex(contact), sp.latex(on_shell_contact))
    print('Charged i M including the incoming pole overlap:')
    print(sp.latex(charged))


if __name__ == '__main__':
    main()
