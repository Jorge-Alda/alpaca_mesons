"""Derive the conditional regular remainder of the retained production numerator.

Run: .venv/bin/python scripts/derive_alp_three_pion_regular_remainder.py
Input: memory/alp/lepton-pair/alp-three-pion-vector-pole-factorization.md.
Natural units; metric (+---); no width or rescattering. This is polynomial
pole subtraction of the retained VJP/VVP vertex, conditional on its off-shell
continuation. It does not determine independent local production operators,
VPPP couplings, or complementary tensor-projector contributions.
"""
import sympy as sp
from derive_alp_three_pion_vector_pole import derive_production_numerator
import derive_radiative_tff_amplitude as tensors


def derive():
    """Obtain each numerator from operator contractions and divide its difference."""
    direct, pair, _ = derive_production_numerator()
    electric_charge_resonance = sp.Symbol('e', real=True)
    vector_mass = sp.Symbol('M_i', positive=True)
    # Each term in the j sum has its own trace and W_j. Keep Sum explicit.
    intermediate = sp.Symbol('j', integer=True)
    response = sp.Function('W')(intermediate, tensors.s2)
    # Derive term by term: the intermediate trace labels are independent of s.
    direct_regular = sp.factor(sp.cancel((direct - direct.subs(tensors.s1, vector_mass**2)) / (vector_mass**2 - tensors.s1)))
    pair_regular = sp.factor(sp.cancel((pair - pair.subs(tensors.s1, vector_mass**2)) / (vector_mass**2 - tensors.s1)))
    regular = -2 * electric_charge_resonance * direct_regular + sp.sqrt(2) * electric_charge_resonance * sp.Sum(response * pair_regular, (intermediate, 1, 3))
    return dict(direct=direct, pair=pair, direct_regular=direct_regular,
                pair_regular=pair_regular, regular=regular)


def main():
    """Print retained input numerators and conditional remainder in LaTeX."""
    print('Natural units; conditional off-shell continuation of retained VJP/VVP vertices.')
    for label, expression in derive().items():
        print(f'{label}:\n{sp.latex(expression)}\n')
    print('Independent local operators and complementary tensor-projector terms remain unspecified.')


if __name__ == '__main__':
    main()
