"""Derive the invariant split of mixed ALP--meson pole overlaps.

Run with ``.venv/bin/python scripts/derive_nlo_lsz_invariant_split.py``.
The matrices describe real neutral meson fields. K and M^2 are symmetric,
K is invertible, and M^2-sK is invertible away from a meson pole. Only
terms linear in the ALP mixing are retained. K is dimensionless; M^2, s
and mu have mass dimension two; k and the shift r are dimensionless.

This script uses the quadratic convention of memory/nlo and the parent's
LSZ convention. It derives an identity for arbitrary local two-point
coefficients, without supplying omitted loop terms or fixing NLO LECs.
"""

import sympy as sp
import sys
import unittest


def show(label, expression):
    print(f"{label}:\n{sp.latex(expression)}\n")


def derive():
    dimension = sp.Symbol("n", integer=True, positive=True)
    momentum_squared = sp.Symbol("s", real=True)
    kinetic_matrix = sp.MatrixSymbol("K", dimension, dimension)
    mass_matrix = sp.MatrixSymbol("M^2", dimension, dimension)
    kinetic_mixing = sp.MatrixSymbol("k", dimension, 1)
    mass_mixing = sp.MatrixSymbol("mu", dimension, 1)
    meson_shift = sp.MatrixSymbol("r", dimension, 1)
    meson_vertex = sp.MatrixSymbol("V", dimension, 1)
    direct_vertex = sp.MatrixSymbol("V_a", 1, 1)
    inverse_meson_propagator = momentum_squared * kinetic_matrix - mass_matrix
    mixed_inverse_propagator = momentum_squared * kinetic_mixing - mass_mixing
    quadratic_matrix = sp.BlockMatrix([
        [inverse_meson_propagator, mixed_inverse_propagator],
        [mixed_inverse_propagator.T, sp.MatrixSymbol("s-m_a^2", 1, 1)],
    ])
    field_conversion = sp.BlockMatrix([
        [sp.Identity(dimension), -meson_shift],
        [sp.ZeroMatrix(1, dimension), sp.Identity(1)],
    ])
    converted_quadratic_matrix = sp.block_collapse(
        field_conversion.T * quadratic_matrix * field_conversion
    )
    converted_mixed_inverse_propagator = sp.expand(converted_quadratic_matrix.blocks[0, 1])
    converted_kinetic_mixing = kinetic_mixing - kinetic_matrix * meson_shift
    converted_mass_mixing = mass_mixing - mass_matrix * meson_shift
    pole_overlap = -inverse_meson_propagator.I * mixed_inverse_propagator
    invariant_source = mass_matrix * kinetic_matrix.I * kinetic_mixing - mass_mixing
    source_split = sp.expand(
        mixed_inverse_propagator
        - inverse_meson_propagator * kinetic_matrix.I * kinetic_mixing
    )
    split_overlap = kinetic_matrix.I * (-kinetic_mixing) - inverse_meson_propagator.I * source_split
    converted_invariant_source = sp.expand(
        mass_matrix * kinetic_matrix.I * converted_kinetic_mixing
        - converted_mass_mixing
    )
    # Keep the full propagator blocks grouped while cancelling their inverse.
    # MatrixExpr expands sums, but cannot regroup a general expanded matrix sum.
    inverse_propagator_symbol = sp.MatrixSymbol("D^{-1}_P", dimension, dimension)
    mixed_propagator_symbol = sp.MatrixSymbol("D^{-1}_{Pa}", dimension, 1)
    grouped_converted_mixing = mixed_propagator_symbol - inverse_propagator_symbol * meson_shift
    grouped_converted_overlap = sp.expand(-inverse_propagator_symbol.I * grouped_converted_mixing)
    converted_overlap = grouped_converted_overlap.xreplace({
        inverse_propagator_symbol: inverse_meson_propagator,
        mixed_propagator_symbol: mixed_inverse_propagator,
    })
    converted_vertex = sp.block_collapse(
        sp.BlockMatrix([[meson_vertex.T, direct_vertex]]) * field_conversion
    ).blocks[0, 1]
    total_vertex = direct_vertex + meson_vertex.T * pole_overlap
    converted_total_vertex = sp.expand(converted_vertex + meson_vertex.T * converted_overlap)
    return locals()


class InvariantSplitTests(unittest.TestCase):
    def test_split_solves_mixed_pole_equation(self):
        result = derive()
        replacement = {result["mass_matrix"]: result["momentum_squared"] * result["kinetic_matrix"] - result["inverse_propagator_symbol"]}
        overlap_difference = result["split_overlap"] - result["pole_overlap"]
        residual = sp.expand(overlap_difference.xreplace(replacement).doit())
        self.assertEqual(residual, sp.ZeroMatrix(result["dimension"], 1))

    def test_invariant_source(self):
        result = derive()
        self.assertEqual(result["converted_invariant_source"], result["invariant_source"])

    def test_shifted_mixing_from_quadratic_form(self):
        result = derive()
        expected = sp.expand(result["momentum_squared"] * result["converted_kinetic_mixing"] - result["converted_mass_mixing"])
        self.assertEqual(result["converted_mixed_inverse_propagator"], expected)

    def test_direct_and_pole_vertex_cancel_shift(self):
        result = derive()
        difference = sp.expand(result["converted_total_vertex"] - result["total_vertex"])
        self.assertEqual(difference, sp.ZeroMatrix(1, 1))


def main():
    result = derive()
    print("Real neutral fields; L_mix = partial(a) k^T partial(P) - a mu^T P.")
    print("Retain first order in ALP mixing; P = P' - r a; s = p^2.")
    show("Original inverse propagator", result["quadratic_matrix"])
    show("Field conversion", result["field_conversion"])
    show("Converted off-diagonal inverse propagator", result["converted_mixed_inverse_propagator"])
    show("Converted kinetic mixing", result["converted_kinetic_mixing"])
    show("Converted mass mixing", result["converted_mass_mixing"])
    show("ALP-pole meson overlap", result["pole_overlap"])
    show("Remainder after separating the kinetic term", result["source_split"])
    show("Invariant source after conversion", result["converted_invariant_source"])
    show("Overlap: one non-invariant term and invariant remainder", result["split_overlap"])
    show("Converted pole overlap", result["converted_overlap"])
    show("Converted direct one-ALP vertex", result["converted_vertex"])
    show("Direct plus meson-pole vertex", result["total_vertex"])
    show("Converted direct plus meson-pole vertex", result["converted_total_vertex"])


if __name__ == "__main__":
    if "--test" in sys.argv:
        unittest.main(argv=[sys.argv[0]])
    else:
        main()
