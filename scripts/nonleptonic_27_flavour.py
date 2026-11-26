#!/usr/bin/env python3
"""Symmetric traceless SU(3) tensor and infinitesimal rotation.

Run: .venv/bin/python scripts/nonleptonic_27_flavour.py [--test]
Natural units; u,d,s ordered 0,1,2; first order in a/f_a.
The tensor indices are upper, upper, lower, lower. Its normalization is
fixed by symmetrizing E13 tensor E21, not by hadronic matching.
"""
import argparse
import unittest
import sympy as sp
from sympy.tensor.array import tensorproduct, permutedims


def calculate():
    first_current = sp.zeros(3); first_current[0, 2] = 1
    second_current = sp.zeros(3); second_current[1, 0] = 1
    # Rearrange (upper, lower, upper, lower) to (upper, upper, lower, lower).
    flavor_tensor = permutedims(tensorproduct(first_current, second_current), (0, 2, 1, 3))
    symmetric_tensor = (flavor_tensor + permutedims(flavor_tensor, (1, 0, 2, 3))
                        + permutedims(flavor_tensor, (0, 1, 3, 2))
                        + permutedims(flavor_tensor, (1, 0, 3, 2))) / 4
    # SymPy tensorcontraction contracts two specified axes.
    trace = sp.tensorcontraction(symmetric_tensor, (1, 3))
    delta = sp.eye(3)
    adjoint_subtraction = sp.ImmutableDenseNDimArray([
        (delta[i,j]*trace[k,l] + delta[i,l]*trace[k,j]
         + delta[k,j]*trace[i,l] + delta[k,l]*trace[i,j]) / 5
        for i in range(3) for k in range(3) for j in range(3) for l in range(3)
    ], (3,3,3,3))
    tensor_27 = symmetric_tensor - adjoint_subtraction
    rotation = sp.Matrix(3, 3, lambda i,j: sp.Symbol(f't_{i+1}{j+1}'))
    # Explicit contractions implement the four index transformation using SymPy matrices.
    rotation_derivative = sp.ImmutableDenseNDimArray([
        sp.I*sum(rotation[i,m]*tensor_27[m,k,j,l] + rotation[k,m]*tensor_27[i,m,j,l]
                 - tensor_27[i,k,m,l]*rotation[m,j] - tensor_27[i,k,j,m]*rotation[m,l]
                 for m in range(3))
        for i in range(3) for k in range(3) for j in range(3) for l in range(3)
    ], (3,3,3,3))
    diagonal = {rotation[i,j]: (sp.Symbol(f't_{i+1}', real=True) if i == j else 0)
                for i in range(3) for j in range(3)}
    return flavor_tensor, symmetric_tensor, trace, tensor_27, rotation, rotation_derivative, diagonal


def nonzero(tensor):
    return [((i+1,k+1,j+1,l+1), sp.simplify(tensor[i,k,j,l]))
            for i in range(3) for k in range(3) for j in range(3) for l in range(3)
            if tensor[i,k,j,l] != 0]


class Tests(unittest.TestCase):
    def test_symmetries_and_traces(self):
        _, _, _, tensor, _, _, _ = calculate()
        self.assertEqual(tensor, permutedims(tensor, (1,0,2,3)))
        self.assertEqual(tensor, permutedims(tensor, (0,1,3,2)))
        self.assertEqual(sp.Matrix(sp.tensorcontraction(tensor,(1,3))), sp.zeros(3))

    def test_diagonal_phase(self):
        _, _, _, tensor, _, derivative, substitutions = calculate()
        t2, t3 = sp.symbols('t_2 t_3', real=True)
        for i in range(3):
            for k in range(3):
                for j in range(3):
                    for l in range(3):
                        self.assertEqual(sp.simplify(derivative[i,k,j,l].subs(substitutions)
                                          - sp.I*(t2-t3)*tensor[i,k,j,l]), 0)

    def test_rotated_trace(self):
        _, _, _, _, _, derivative, _ = calculate()
        self.assertEqual(sp.Matrix(sp.tensorcontraction(derivative,(1,3))).applyfunc(sp.simplify), sp.zeros(3))


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--test', action='store_true')
    if parser.parse_args().test:
        unittest.main(argv=['test']); return
    _, symmetric, trace, tensor, _, derivative, diagonal = calculate()
    print('Natural units; first order in a/f_a; tensor indices (upper, upper, lower, lower).')
    print('Symmetric tensor nonzero entries:')
    for indices, value in nonzero(symmetric): print(indices, sp.latex(value))
    print('Single trace:'); print(sp.latex(sp.Matrix(trace)))
    print('Symmetric traceless 27-plet nonzero entries:')
    for indices, value in nonzero(tensor): print(indices, sp.latex(value))
    print('Derivative under a diagonal rotation for each nonzero entry:')
    for (i,k,j,l), _ in nonzero(tensor):
        print((i,k,j,l), sp.latex(sp.simplify(derivative[i-1,k-1,j-1,l-1].subs(diagonal))))


if __name__ == '__main__': main()
