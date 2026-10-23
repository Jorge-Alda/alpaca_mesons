"""Derive the original-eleven-operator neutral ALP two-photon scalar.

Invocation: .venv/bin/python scripts/derive_alp_doubly_virtual_tff.py

Inputs: memory/alp/photon/alp-photon-vector-mixing.md,
memory/rcht/odd-parity/radiative-tff-amplitude.md, and
memory/rcht/odd-parity/odd-alp-source-lsz-assembly.md.
Natural units, GeV, metric (+---), epsilon^{0123}=+1. Canonical c_i,d_i.
The final scalar F uses Messerli:2025rnv, Eq. (9): i e^2 R F,
where R=epsilon(epsilon1,epsilon2,q1,q2). Intermediate ordered numerators
use E=epsilon(q1,epsilon1,q2,epsilon2)=-R, so F=-G.
Assumptions: tree level, one ALP, diagonal neutral sources, original odd
operators, stable neutral vectors and their signed photon residues.
There are no hadronic cuts, widths, fitted parameter choices or UV matching.
The existing Lorentz-contraction backend is used only for normalized tensors;
their ratios are independent of the epsilon orientation. Its c_i,d_i are
never imported as physical couplings without the canonical source prefactors.
"""

import sympy as sp

import derive_radiative_tff_amplitude as tensor


def ordered_numerators(vector_axial_trace, vector_mass_trace):
    """Contract the VJP and VVP operators with independent source traces.

    For VJP the traces are tr(T_i Q U_a), tr(T_i Q Y_a).
    For VVP they are tr(T_i T_j U_a), tr(T_i T_j Y_a).
    U_a has dimension GeV^{-1}; Y_a has dimension GeV.
    """
    first = tensor.contractions(tensor.q1, tensor.q2, tensor.e1, tensor.e2)
    second = tensor.contractions(tensor.q2, tensor.q1, tensor.e2, tensor.e1)
    vjp = -2 / tensor.MV * sum(
        tensor.c[index]
        * (vector_mass_trace if index == 2 else vector_axial_trace)
        * first[index]
        for index in range(7)
    )
    vvp = -2 * (
        tensor.d[0] * vector_axial_trace * (first[0] + second[0])
        + tensor.d[1] * vector_mass_trace * (first[2] + second[2])
        + tensor.d[2] * vector_axial_trace * (first[4] + second[4])
        + tensor.d[3] * vector_axial_trace * (first[6] + second[6])
    )
    return sp.expand(vjp), sp.expand(vvp)


def photon_response():
    """Solve the quadratic tensor equation with a signed residue R_i."""
    response, residue, mass, virtuality = sp.symbols("r R_i M_i s", real=True)
    quadratic_lagrangian = (mass**2 - virtuality) * response**2 - 2 * residue * response
    solution = sp.solve(sp.diff(quadratic_lagrangian, response), response)[0]
    return quadratic_lagrangian, solution


def assemble_scalar(vector_masses, photon_residues, vjp_numerators, vvp_numerators):
    """Assemble F in the Messerli:2025rnv photon tensor convention.

    Inputs are symbolic sequences; VJP entries carry (s_1,s_2), and the
    VVP matrix is symmetric for diagonal sources. Do not interchange an
    already inserted pole: interchange the VJP numerator first. Ordered
    numerators are E-normalized; the returned terms include F=-G.
    """
    first_responses = [residue / (mass**2 - tensor.s1)
                       for mass, residue in zip(vector_masses, photon_residues)]
    second_responses = [residue / (mass**2 - tensor.s2)
                        for mass, residue in zip(vector_masses, photon_residues)]
    vjp = 2 * sp.sqrt(2) * sum(
        first_responses[index] * numerator
        + second_responses[index] * numerator.xreplace({tensor.s1: tensor.s2, tensor.s2: tensor.s1})
        for index, numerator in enumerate(vjp_numerators)
    )
    vvp = -2 * sum(
        first_responses[first_index] * second_responses[second_index] * numerator
        for first_index, row in enumerate(vvp_numerators)
        for second_index, numerator in enumerate(row)
    )
    return vjp, vvp


def show(label, expression):
    """Print deterministic LaTeX with the standard SymPy printer."""
    print(f"{label}:\n{sp.latex(expression)}\n")


def main():
    """Print the symbolic source construction without choosing UV couplings."""
    print("Tree level; original eleven odd operators; diagonal no-mixing sources.")
    print("GeV; (+---); epsilon^{0123}=+1; i M = i e^2 R F.")
    print("Messerli:2025rnv Eq. (9): R=epsilon(epsilon1,epsilon2,q1,q2).")
    print("Ordered numerators use E=epsilon(q1,epsilon1,q2,epsilon2)=-R; F=-G.")
    print("Stable vector poles; low-energy construction without hadronic cuts or UV matching.\n")
    quadratic_lagrangian, response = photon_response()
    show("Quadratic tensor Lagrangian after removing F_gamma squared/4", quadratic_lagrangian)
    show("Signed W_i(s) from its field equation", response)
    axial_trace = sp.Symbol(r"\operatorname{tr}(T_i Q U_a)")
    mass_trace = sp.Symbol(r"\operatorname{tr}(T_i Q Y_a)")
    vjp, _ = ordered_numerators(axial_trace, mass_trace)
    pair_axial_trace = sp.Symbol(r"\operatorname{tr}(T_i T_j U_a)")
    pair_mass_trace = sp.Symbol(r"\operatorname{tr}(T_i T_j Y_a)")
    _, vvp = ordered_numerators(pair_axial_trace, pair_mass_trace)
    alp_mass = sp.Symbol("m_a", positive=True)
    show("Ordered V_i gamma(s1,s2), with the vector on photon 1", sp.collect(vjp.subs(tensor.m, alp_mass), [axial_trace, mass_trace]))
    show("Ordered P_ij(s1,s2)", sp.collect(vvp.subs(tensor.m, alp_mass), [pair_axial_trace, pair_mass_trace]))
    vector_index, other_index = sp.symbols("i j", integer=True)
    response_function = sp.Function("W")
    vjp_function = sp.Function(r"V_{i\gamma}")
    vvp_function = sp.Function(r"P_{ij}")
    vjp_sum = 2 * sp.sqrt(2) * sp.Sum(
        response_function(vector_index, tensor.s1) * vjp_function(tensor.s1, tensor.s2)
        + response_function(vector_index, tensor.s2) * vjp_function(tensor.s2, tensor.s1),
        (vector_index, 1, 3),
    )
    vvp_sum = -2 * sp.Sum(
        response_function(vector_index, tensor.s1)
        * response_function(other_index, tensor.s2)
        * vvp_function(tensor.s1, tensor.s2),
        (vector_index, 1, 3), (other_index, 1, 3),
    )
    decay_constant = sp.Symbol("f_a", positive=True)
    contact_coefficient = sp.Symbol(r"\widehat c_\gamma", real=True)
    show("F_VJP from two photon placements", vjp_sum)
    show("F_VVP from ordered vector assignments", vvp_sum)
    show("Complete retained scalar F in the Messerli:2025rnv R tensor convention", contact_coefficient / (4 * sp.pi**2 * decay_constant) + vjp_sum + vvp_sum)
    show("Ordered V_i gamma(0,0)", sp.factor(vjp.subs({tensor.m: alp_mass, tensor.s1: 0, tensor.s2: 0})))
    show("Ordered P_ij(0,0)", sp.factor(vvp.subs({tensor.m: alp_mass, tensor.s1: 0, tensor.s2: 0})))


if __name__ == "__main__":
    main()
