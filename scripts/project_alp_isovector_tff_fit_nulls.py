"""Project preliminary-fit coefficient nulls onto the rho ALP TFF.

Invocation: .venv/bin/python scripts/project_alp_isovector_tff_fit_nulls.py

Inputs: memory/rcht/odd-parity/preliminary-fit-lagrangian-and-null-directions.md
and memory/alp/lepton-pair/alp-isovector-doubly-virtual-tff.md.
Natural units, GeV, metric (+---), epsilon^{0123}=+1; canonical c_i,d_i.
Tree level, original eleven odd operators, diagonal no-mixing ALP sources,
fixed even-sector parameters and source residues. No short-distance condition
or relation between the independent axial and pseudoscalar ALP sources is used.
The unrounded numerical basis is reconstructed from the stored fit using its
original exporter without writing a fit file. The symbolic dependence is
derived from the existing operator-contraction calculation.
"""

import numpy as np
import scipy.linalg as linear_algebra
import sympy as sp

import derive_alp_doubly_virtual_tff as tff
import summarize_odd_preliminary_fit as fit


def show(label, expression):
    """Print a labelled expression with the standard SymPy LaTeX printer."""
    print(f"{label}:\n{sp.latex(expression)}\n")


def main():
    """Calculate the three numerator projections and their functional rank."""
    _, null_basis = fit.main(write_output=False)
    coefficients = sp.Matrix([*tff.tensor.c, *tff.tensor.d])
    axial_trace, mass_trace = sp.symbols("A_rho B_rho", real=True)
    pair_axial_trace, pair_mass_trace = sp.symbols("C_rhorho D_rhorho", real=True)
    mixed_axial_trace, mixed_mass_trace = sp.symbols("C_rhoI D_rhoI", real=True)
    vjp, _ = tff.ordered_numerators(axial_trace, mass_trace)
    _, rho_pair = tff.ordered_numerators(pair_axial_trace, pair_mass_trace)
    _, mixed_pair = tff.ordered_numerators(mixed_axial_trace, mixed_mass_trace)
    numerators = sp.Matrix([vjp, rho_pair, mixed_pair])
    derivative_matrix = numerators.jacobian(coefficients)
    show("Ordered numerators (V_rho gamma, P_rho rho, P_rho I)", numerators)
    np.set_printoptions(precision=12, suppress=True, linewidth=140)
    print("Unrounded regenerated numerical null basis; rows c1,...,c7,d1,...,d4:")
    print(null_basis)

    # Independent polynomial coefficients determine functional dependence;
    # distinct monomials in virtualities and source traces must not be merged.
    polynomial_variables = [axial_trace, mass_trace, pair_axial_trace,
                            pair_mass_trace, mixed_axial_trace, mixed_mass_trace,
                            tff.tensor.m, tff.tensor.s1, tff.tensor.s2, tff.tensor.MV]
    polynomial_rows = []
    for numerator in numerators:
        polynomial = sp.Poly(sp.expand(numerator * tff.tensor.MV), *polynomial_variables)
        polynomial_rows.extend([sp.Matrix([coefficient]).jacobian(coefficients)
                                for _, coefficient in polynomial.terms()])
    functional_matrix = sp.Matrix.vstack(*polynomial_rows)
    numerical_projection = np.asarray(functional_matrix, dtype=float) @ null_basis
    singular_values = linear_algebra.svdvals(numerical_projection)
    numerical_rank = int(np.count_nonzero(singular_values > 1e-9))
    print("Projected independent polynomial coefficients (numerators multiplied by M_V):")
    print(numerical_projection)
    print("Functional singular values:", singular_values)
    print("Functional rank (absolute tolerance 1e-9):", numerical_rank)

    # Recover exact identities obeyed by this null space, rather than using
    # rounded SVD entries as exact symbolic numbers.
    identity_rows = np.zeros((7, 11))
    for row, column in enumerate([2, 3, 7, 8, 9]):
        identity_rows[row, column] = 1
    identity_rows[5, [0, 1, 4]] = [-1, -1, 1]
    identity_rows[6, [1, 5]] = [-1, 1]
    identity_defect = float(np.max(np.abs(identity_rows @ null_basis)))
    print("Maximum residual of exact null identities:", identity_defect)
    first_shift, second_shift, seventh_shift, fourth_pair_shift = sp.symbols(
        "delta_c1 delta_c2 delta_c7 delta_d4", real=True)
    exact_shift = sp.Matrix([first_shift, second_shift, 0, 0,
                            first_shift + second_shift, second_shift,
                            seventh_shift, 0, 0, 0, fourth_pair_shift])
    show("Exact arbitrary null shift", exact_shift)
    projected_numerators = sp.simplify(derivative_matrix * exact_shift)
    show("Exact ordered numerator shifts", projected_numerators)
    print("delta c1 in the regenerated alpha basis:", null_basis[0])
    first_response = sp.Function("W_rho")(tff.tensor.s1)
    second_response = sp.Function("W_rho")(tff.tensor.s2)
    numerator_shift = projected_numerators[0]
    scalar_shift = 2 * sp.sqrt(2) * (
        first_response * numerator_shift
        + second_response * numerator_shift.xreplace(
            {tff.tensor.s1: tff.tensor.s2, tff.tensor.s2: tff.tensor.s1}))
    show("Assembled isovector TFF shift", sp.factor(scalar_shift))
    show("TFF shift for two real photons", scalar_shift.subs(
        {tff.tensor.s1: 0, tff.tensor.s2: 0}))
    print("Scope: one visible null combination for general axial ALP sources;")
    print("rho-rho and rho-isoscalar VVP terms have zero null dependence.")


if __name__ == "__main__":
    main()
