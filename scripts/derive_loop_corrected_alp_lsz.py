"""Derive the first-order loop-corrected ALP--meson pole overlap.

Run ``.venv/bin/python scripts/derive_loop_corrected_alp_lsz.py``.

The convention is Q(s) = s K - M^2 + Sigma(s).  The meson block of
Q^(0) is invertible at the ALP pole.  Only first order in ALP--meson
mixing and first order in the NLO corrections are retained.  No matrix
commutativity is assumed.  K and the meson components of w are
dimensionless, while s, M^2, Q and Sigma have units of mass squared.
The script does not calculate the loop self energies.
"""

import sympy as sp


def show(label, expression, symbol_names):
    print(f"{label}:\n{sp.latex(expression, symbol_names=symbol_names)}\n")


def derive():
    dimension = sp.Symbol("n", integer=True, positive=True)
    momentum_squared = sp.Symbol("s", real=True)
    order = sp.Symbol("epsilon", real=True)
    pole_shift = sp.Symbol("delta_s_a", real=True)

    leading_kinetic_matrix = sp.MatrixSymbol("K_PP^(0)", dimension, dimension)
    leading_mass_matrix = sp.MatrixSymbol("M_PP^2(0)", dimension, dimension)
    leading_kinetic_mixing = sp.MatrixSymbol("K_Pa^(0)", dimension, 1)
    leading_mass_mixing = sp.MatrixSymbol("M_Pa^2(0)", dimension, 1)
    local_kinetic_matrix = sp.MatrixSymbol("delta_K_PP", dimension, dimension)
    local_mass_matrix = sp.MatrixSymbol("delta_M_PP^2", dimension, dimension)
    local_kinetic_mixing = sp.MatrixSymbol("delta_K_Pa", dimension, 1)
    local_mass_mixing = sp.MatrixSymbol("delta_M_Pa^2", dimension, 1)
    loop_meson_self_energy = sp.MatrixSymbol("Sigma_PP^loop(s)", dimension, dimension)
    loop_mixed_self_energy = sp.MatrixSymbol("Sigma_Pa^loop(s)", dimension, 1)

    leading_meson_block = sp.MatrixSymbol("Q_PP^(0)", dimension, dimension)
    leading_mixed_block = sp.MatrixSymbol("Q_Pa^(0)", dimension, 1)
    correction_meson_block = sp.MatrixSymbol("delta_Q_PP", dimension, dimension)
    correction_mixed_block = sp.MatrixSymbol("delta_Q_Pa", dimension, 1)
    leading_overlap_symbol = sp.MatrixSymbol("w_P^(a,0)", dimension, 1)
    correction_overlap_symbol = sp.MatrixSymbol("delta_w_P^(a)", dimension, 1)

    leading_meson_definition = momentum_squared * leading_kinetic_matrix - leading_mass_matrix
    leading_mixed_definition = momentum_squared * leading_kinetic_mixing - leading_mass_mixing
    correction_meson_definition = (
        momentum_squared * local_kinetic_matrix - local_mass_matrix + loop_meson_self_energy
    )
    correction_mixed_definition = (
        momentum_squared * local_kinetic_mixing - local_mass_mixing + loop_mixed_self_energy
    )

    pole_equation = (
        (leading_meson_block + order * correction_meson_block)
        * (leading_overlap_symbol + order * correction_overlap_symbol)
        + leading_mixed_block + order * correction_mixed_block
    )
    leading_equation = sp.expand(pole_equation.subs(order, 0))
    first_order_equation = sp.expand(sp.diff(pole_equation, order).subs(order, 0))
    leading_overlap = -leading_meson_block.I * leading_equation.subs(
        leading_overlap_symbol, sp.ZeroMatrix(dimension, 1)
    )
    correction_overlap = -leading_meson_block.I * first_order_equation.subs(
        correction_overlap_symbol, sp.ZeroMatrix(dimension, 1)
    )

    # Differentiate the leading pole equation implicitly in s.  Unlike a
    # component calculation, this retains the ordering of the meson matrices.
    shifted_leading_equation = (
        (leading_meson_block + order * pole_shift * leading_kinetic_matrix)
        * (leading_overlap_symbol + order * correction_overlap_symbol)
        + leading_mixed_block + order * pole_shift * leading_kinetic_mixing
    )
    first_order_shift_equation = sp.expand(
        sp.diff(shifted_leading_equation, order).subs(order, 0)
    )
    pole_shift_overlap = -leading_meson_block.I * first_order_shift_equation.subs(
        correction_overlap_symbol, sp.ZeroMatrix(dimension, 1)
    )
    return locals()


def main():
    result = derive()
    symbol_names = {
        result["leading_kinetic_matrix"]: r"K_{PP}^{(0)}",
        result["leading_mass_matrix"]: r"M_{PP}^{2(0)}",
        result["leading_kinetic_mixing"]: r"K_{P\hat a}^{(0)}",
        result["leading_mass_mixing"]: r"M_{P\hat a}^{2(0)}",
        result["local_kinetic_matrix"]: r"\delta K_{PP}^{\mathrm{local}}",
        result["local_mass_matrix"]: r"\delta M_{PP}^{2,\mathrm{local}}",
        result["local_kinetic_mixing"]: r"\delta K_{P\hat a}^{\mathrm{local}}",
        result["local_mass_mixing"]: r"\delta M_{P\hat a}^{2,\mathrm{local}}",
        result["loop_meson_self_energy"]: r"\Sigma_{PP}^{\mathrm{loop}}(s)",
        result["loop_mixed_self_energy"]: r"\Sigma_{P\hat a}^{\mathrm{loop}}(s)",
        result["leading_meson_block"]: r"\mathcal Q_{PP}^{(0)}",
        result["leading_mixed_block"]: r"\mathcal Q_{P\hat a}^{(0)}",
        result["correction_meson_block"]: r"\delta\mathcal Q_{PP}",
        result["correction_mixed_block"]: r"\delta\mathcal Q_{P\hat a}",
        result["leading_overlap_symbol"]: r"w_{P}^{(a,0)}",
        result["correction_overlap_symbol"]: r"\delta w_{P}^{(a)}",
        result["pole_shift"]: r"\delta s_a",
    }
    print("Q(s) = s K - M^2 + Sigma(s); real fields and w_a^(a) = 1.")
    print("First order in ALP mixing and NLO corrections; Q_PP^(0) is invertible.")
    print("K and w are dimensionless; s, M^2, Q and Sigma have mass dimension two.\n")
    show("Leading meson block", result["leading_meson_definition"], symbol_names)
    show("Leading ALP--meson block", result["leading_mixed_definition"], symbol_names)
    show("NLO meson block", result["correction_meson_definition"], symbol_names)
    show("NLO ALP--meson block", result["correction_mixed_definition"], symbol_names)
    show("Leading pole equation, equal to zero", result["leading_equation"], symbol_names)
    show("First-order pole equation, equal to zero", result["first_order_equation"], symbol_names)
    show("Leading meson overlap", result["leading_overlap"], symbol_names)
    show("NLO meson overlap at fixed s", result["correction_overlap"], symbol_names)
    show("Pole-shift contribution to the meson overlap", result["pole_shift_overlap"], symbol_names)


if __name__ == "__main__":
    main()
