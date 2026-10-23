"""Derive the two-photon cut and decay normalization with SymPy.

Run: .venv/bin/python scripts/derive_alp_lepton_pair_normalization.py
Inputs: Messerli:2025rnv Eqs. (6), (9), (13), and the indexed ALP TFF
record. Metric (+---); epsilon^{0123}=+1; masses in GeV, F in GeV^-1.
Leading QED order, on-shell photons, and real external s>4 m_l^2 are assumed.
In widths F_agg denotes |F(0,0)|. Overall interference phases are not derived.
Only the gamma-gamma discontinuity is evaluated; no hadronic model is used.
"""

import sympy as sp
from sympy.physics.matrices import mgamma


def derive_normalization() -> dict[str, sp.Expr]:
    """Calculate spin sums, widths, and the angular two-photon cut."""
    invariant_mass_squared, lepton_mass, alpha, form_factor = sp.symbols(
        "s m_l alpha_em F_agg", positive=True
    )
    velocity = sp.Symbol("sigma", positive=True)
    cosine = sp.Symbol("z", real=True)
    coupling = sp.Symbol("g_all", real=True)
    reduced_amplitude_squared = sp.Symbol("Abs_A_squared", nonnegative=True)
    parent_mass = sp.sqrt(invariant_mass_squared)
    lepton_energy = parent_mass / 2
    lepton_momentum = sp.sqrt(invariant_mass_squared / 4 - lepton_mass**2)
    lepton_slash = lepton_energy * mgamma(0) - lepton_momentum * mgamma(3)
    antilepton_slash = lepton_energy * mgamma(0) + lepton_momentum * mgamma(3)
    vertex = sp.I * mgamma(5)
    conjugate_vertex = mgamma(0) * vertex.conjugate().T * mgamma(0)
    spin_sum = sp.simplify(sp.trace(
        (lepton_slash + lepton_mass * sp.eye(4)) * vertex
        * (antilepton_slash - lepton_mass * sp.eye(4)) * conjugate_vertex
    ))

    # Two transverse linear polarizations on each photon, with q1 along +z.
    photon_one = sp.Matrix([parent_mass / 2, 0, 0, parent_mass / 2])
    photon_two = sp.Matrix([parent_mass / 2, 0, 0, -parent_mass / 2])
    polarizations = [sp.Matrix([0, 1, 0, 0]), sp.Matrix([0, 0, 1, 0])]
    tensor_squared_sum = sp.simplify(sum(
        sp.det(sp.Matrix.hstack(first, second, photon_one, photon_two))**2
        for first in polarizations for second in polarizations
    ))
    photon_spin_sum = (4 * sp.pi * alpha)**2 * form_factor**2 * tensor_squared_sum
    # Integrated two-body phase space: |p|/(4 pi m_a); width includes 1/(2m_a).
    photon_width = sp.simplify(photon_spin_sum / (32 * sp.pi * parent_mass))
    lepton_width = sp.simplify(coupling**2 * spin_sum * velocity / (16 * sp.pi * parent_mass))
    # Eqs. (6), (13) fix |g_all|=2 alpha_em^2 m_l |F_agg A_l|.
    # This width calculation does not determine a phase relative to a local term.
    tff_lepton_width = sp.simplify(lepton_width.subs(
        coupling**2, 4 * alpha**4 * lepton_mass**2 * form_factor**2 * reduced_amplitude_squared
    ))
    width_ratio = sp.simplify(tff_lepton_width / photon_width)

    # Cutkosky: Im A = -4/s int d^4k delta_+(k^2) delta_+((q-k)^2)
    # [s k^2-(q.k)^2]/[(p-k)^2-m_l^2]. The radial delta integral is 1/8.
    on_shell_numerator = -invariant_mass_squared**2 / 4
    lepton_denominator = -invariant_mass_squared * (1 - velocity * cosine) / 2
    angular_integrand = sp.simplify(
        -4 / invariant_mass_squared * 2 * sp.pi / 8
        * on_shell_numerator / lepton_denominator
    )
    # Integrate with real 0<sigma<1 endpoints, avoiding complex log branches.
    angular_primitive = sp.integrate(angular_integrand, cosine)
    # On this interval, log(2 sigma z-2)=log(1-sigma z)+log(2)+i pi.
    # Discard only the z-independent integration constant from SymPy's result.
    real_primitive = angular_primitive.xreplace({
        sp.log(2 * velocity * cosine - 2): sp.log(1 - velocity * cosine)
    })
    absorptive_part = sp.simplify(real_primitive.subs(cosine, 1) - real_primitive.subs(cosine, -1))
    photon_cut_ratio = sp.simplify(width_ratio.subs(reduced_amplitude_squared, absorptive_part**2))
    return {
        "Lepton spin sum divided by |g_all|^2": spin_sum,
        "Photon tensor polarization sum": tensor_squared_sum,
        "Two-photon width (identical-particle factor included)": photon_width,
        "Lepton width": lepton_width,
        "TFF lepton width": tff_lepton_width,
        "Width ratio": width_ratio,
        "Two-photon cut angular integrand": angular_integrand,
        "SymPy angular primitive": angular_primitive,
        "Real angular primitive for 0 < sigma < 1": real_primitive,
        "Two-photon absorptive reduced amplitude": absorptive_part,
        "Two-photon-cut-only width ratio": photon_cut_ratio,
        "Threshold absorptive limit": sp.limit(absorptive_part, velocity, 0, dir="+"),
    }


def main() -> None:
    """Print the derivation in deterministic LaTeX form."""
    print("Metric (+---); masses GeV; F_agg GeV^-1; leading QED; s > 4 m_l^2.")
    print("sigma = sqrt(1-4*m_l^2/s). Only the two-real-photon cut is computed.")
    print("F_agg in widths denotes |F(0,0)|; no local-term interference phase is fixed.")
    print("The full real amplitude needs ultraviolet matching; other cuts can interfere.")
    for label, expression in derive_normalization().items():
        print(f"{label}:\n{sp.latex(expression)}\n")


if __name__ == "__main__":
    main()
