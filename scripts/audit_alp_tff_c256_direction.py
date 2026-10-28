"""Audit the common c2,c5,c6 shift in the retained neutral ALP TFF.

Invocation: .venv/bin/python scripts/audit_alp_tff_c256_direction.py

Inputs and conventions: memory/alp/lepton-pair/alp-doubly-virtual-tff.md
and scripts/derive_alp_doubly_virtual_tff.py (Alda:2026nrd sources,
Messerli:2025rnv photon tensor). Natural units, GeV, metric (+---),
epsilon^{0123}=+1; F=-G has dimension GeV^{-1}. U_a has dimension
GeV^{-1}, Y_a GeV, s1,s2 GeV^2, and signed vector residues GeV.

Scope: the retained tree-level original seven VJP/four VVP operators,
diagonal neutral sources, three stable neutral-vector poles, fixed contact
coefficient. No additional approximation is made. Source traces and vector
residues remain independent symbols; no source EOM, equal-spurion relation,
short-distance constraint, numerical input, or real-photon limit is imposed.
The identity is a rational-function identity, not an evaluation at a pole.
Unit tests run in this file; no output files are written.
"""

from functools import lru_cache
import unittest

import sympy as sp

import derive_alp_doubly_virtual_tff as retained


def common_shift(expression):
    """Return the exact finite variation with only c2,c5,c6 shifted."""
    shift = sp.Symbol("delta_c", real=True)
    substitutions = {
        retained.tensor.c[index]: retained.tensor.c[index] + shift
        for index in (1, 4, 5)
    }
    return sp.simplify(sp.expand(expression.subs(substitutions, simultaneous=True) - expression))


@lru_cache(maxsize=1)
def derive_variations():
    """Contract the retained operators and assemble both photon placements.

    Independent trace symbols establish the identity coefficient by
    coefficient in U_a and Y_a, without choosing their flavour entries.
    The ordered VVP sum includes all nine pairs, including diagonal pairs.
    """
    axial_trace, mass_trace = sp.symbols("U_i_gamma Y_i_gamma")
    pair_axial_trace, pair_mass_trace = sp.symbols("U_ij Y_ij")
    vjp_numerator, _ = retained.ordered_numerators(axial_trace, mass_trace)
    _, vvp_numerator = retained.ordered_numerators(pair_axial_trace, pair_mass_trace)
    vector_masses = sp.symbols("M_rho M_omega M_phi", positive=True)
    photon_residues = sp.symbols("R_rho R_omega R_phi")
    axial_traces = sp.symbols("U_rho_gamma U_omega_gamma U_phi_gamma")
    mass_traces = sp.symbols("Y_rho_gamma Y_omega_gamma Y_phi_gamma")
    vjp_numerators = [
        vjp_numerator.subs({axial_trace: axial_value, mass_trace: mass_value})
        for axial_value, mass_value in zip(axial_traces, mass_traces)
    ]
    vvp_numerators = [
        [
            vvp_numerator.subs({
                pair_axial_trace: sp.Symbol(f"U_{min(first, second)}_{max(first, second)}"),
                pair_mass_trace: sp.Symbol(f"Y_{min(first, second)}_{max(first, second)}"),
            })
            for second in range(3)
        ]
        for first in range(3)
    ]
    vjp, vvp = retained.assemble_scalar(
        vector_masses, photon_residues, vjp_numerators, vvp_numerators
    )
    contact = sp.Symbol("c_hat_gamma") / (4 * sp.pi**2 * sp.Symbol("f_a", positive=True))
    swapped_numerator = vjp_numerator.xreplace({
        retained.tensor.s1: retained.tensor.s2,
        retained.tensor.s2: retained.tensor.s1,
    })
    return {
        "VJP numerator": vjp_numerator,
        "VVP numerator": vvp_numerator,
        "delta VJP numerator": common_shift(vjp_numerator),
        "delta swapped VJP numerator": common_shift(swapped_numerator),
        "delta VVP numerator": common_shift(vvp_numerator),
        "delta F_VJP (both placements)": common_shift(vjp),
        "delta F_VVP (ordered pairs)": common_shift(vvp),
        "delta F (including fixed contact)": common_shift(contact + vjp + vvp),
    }


class CommonCoefficientDirectionTests(unittest.TestCase):
    """Exact symbolic tests, with a nonzero individual-coefficient control."""

    def test_each_placement_and_full_tff(self):
        """Require zero finite variation for arbitrary independent inputs."""
        for label, expression in derive_variations().items():
            if label.startswith("delta "):
                with self.subTest(expression=label):
                    self.assertEqual(expression, 0)

    def test_individual_coefficients_are_not_discarded(self):
        """Check each derivative and their sum against the contracted result."""
        numerator = derive_variations()["VJP numerator"]
        axial_trace = sp.Symbol("U_i_gamma")
        expected_derivatives = (
            2 * axial_trace / retained.tensor.MV
            * (retained.tensor.m**2 + retained.tensor.s1 - retained.tensor.s2),
            2 * axial_trace / retained.tensor.MV
            * (-retained.tensor.m**2 + retained.tensor.s1 + retained.tensor.s2),
            -4 * axial_trace * retained.tensor.s1 / retained.tensor.MV,
        )
        derivatives = []
        for index, expected in zip((1, 4, 5), expected_derivatives):
            derivative = sp.diff(numerator, retained.tensor.c[index])
            self.assertEqual(sp.expand(derivative - expected), 0)
            self.assertNotEqual(sp.expand(derivative), 0)
            derivatives.append(derivative)
        self.assertEqual(sp.expand(sum(derivatives)), 0)


def main():
    """Print the derivation and exit unsuccessfully if a unit test fails."""
    print("Retained tree TFF; independent U_a,Y_a and signed vector residues.", flush=True)
    for label, expression in derive_variations().items():
        print(f"{label}: {sp.latex(expression)}", flush=True)
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(CommonCoefficientDirectionTests)
    result = unittest.TextTestRunner(verbosity=1).run(suite)
    if not result.wasSuccessful():
        raise SystemExit(1)
    print("PASS: delta c2 = delta c5 = delta c6 gives identically zero delta F.")


if __name__ == "__main__":
    main()
