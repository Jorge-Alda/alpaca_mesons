"""Construct the conditional isovector ALP TFF double-spectral part with SymPy.

Run: .venv/bin/python scripts/derive_alp_isovector_tff_spectral.py

Inputs: memory/alp/lepton-pair/alp-isovector-tff-spectral-construction.md;
Messerli:2025rnv Eqs. (5.13)-(5.14). The density prefactor is a literature
input, not an independently derived ALP unitarity normalization. Matching the
ALP P wave to that normalization and phase remains required.

Metric (+---), epsilon^{0123}=+1, i M=i e^2 R F; F has units GeV^-1.
x,y,Q_i^2,M_pi^2,Lambda^2 have units GeV^2; F_pi^V is dimensionless;
the full P wave has units GeV^-3; rho has units GeV^-1. Source coefficients
are real, dimensionless placeholders, with 1/f_a included once in each
source P wave. No independent physical source basis is selected here.
Assumptions: isospin limit, two-pion P wave, leading electromagnetic order,
spacelike photons, sharp common cutoff. No hadronic amplitudes are supplied.
Single-spectral terms and subtraction functions may also contribute to the
low-energy isovector TFF and are not determined by this calculation.
"""

import sympy as sp


def derive_spectral_construction() -> dict[str, sp.Expr]:
    """Insert real source coefficients, then build both photon placements.

    Two abstract source responses demonstrate linearity for an arbitrary
    finite real source basis. Taking Im does not commute with multiplication
    by complex coefficients; those would require a discontinuity formulation.
    The common rectangular integration domain is invariant under x<->y.
    """
    dipion_mass_squared, photon_mass_squared = sp.symbols("x y", positive=True)
    pion_mass, cutoff = sp.symbols("M_pi Lambda", positive=True)
    photon_one_spacelike, photon_two_spacelike = sp.symbols(
        "Q_1_squared Q_2_squared", nonnegative=True
    )
    source_one, source_two = sp.symbols("g_1 g_2", real=True)
    pion_vector_form_factor = sp.Function("F_pi_V")(dipion_mass_squared)
    source_one_p_wave = sp.Function(r"\mathcal{F}_1")(
        dipion_mass_squared, photon_mass_squared
    )
    source_two_p_wave = sp.Function(r"\mathcal{F}_2")(
        dipion_mass_squared, photon_mass_squared
    )
    source_p_wave = source_one * source_one_p_wave + source_two * source_two_p_wave
    pion_velocity = sp.sqrt(1 - 4 * pion_mass**2 / dipion_mass_squared)
    # This numerical coefficient is the stated literature-normalization input.
    phase_space_factor = dipion_mass_squared * pion_velocity**3 / (192 * sp.pi)
    # Differentiate the explicitly linear source amplitude to extract each
    # response. Keep Im unevaluated so the printed physical products remain
    # readable, instead of expanding into their separate real/imaginary parts.
    source_density = sp.Add(*(
        coefficient * sp.im(
            sp.conjugate(pion_vector_form_factor)
            * sp.diff(source_p_wave, coefficient), evaluate=False
        )
        for coefficient in (source_one, source_two)
    ))
    density = phase_space_factor * source_density
    ordered_denominator = (
        (dipion_mass_squared + photon_one_spacelike)
        * (photon_mass_squared + photon_two_spacelike)
    )
    ordered_integrand = density / (sp.pi**2 * ordered_denominator)
    photon_exchange = {
        photon_one_spacelike: photon_two_spacelike,
        photon_two_spacelike: photon_one_spacelike,
    }
    symmetric_integrand = ordered_integrand + ordered_integrand.xreplace(photon_exchange)
    integration_limits = (
        (dipion_mass_squared, 4 * pion_mass**2, cutoff**2),
        (photon_mass_squared, 4 * pion_mass**2, cutoff**2),
    )
    form_factor = sp.Integral(symmetric_integrand, *integration_limits)
    symmetric_density = (
        density + density.xreplace({
            dipion_mass_squared: photon_mass_squared,
            photon_mass_squared: dipion_mass_squared,
        })
    )
    # Relabel x<->y in the exchanged-photon integral. No factor 1/2 occurs.
    relabelled_form_factor = sp.Integral(
        symmetric_density / (sp.pi**2 * ordered_denominator), *integration_limits
    )
    return {
        "Source P wave": source_p_wave,
        "Two-pion velocity": pion_velocity,
        "Literature P-wave phase-space prefactor": phase_space_factor,
        "Source-linear ordered density": density,
        "Ordered dispersion integrand": ordered_integrand,
        "Bose-symmetric dispersion integrand": symmetric_integrand,
        "Low-energy isovector double-spectral contribution": form_factor,
        "Density after relabelling the exchanged integral": symmetric_density,
        "Equivalent double-spectral integral over the symmetrized density": relabelled_form_factor,
        "Real-photon integrand": sp.simplify(symmetric_integrand.subs({
            photon_one_spacelike: 0, photon_two_spacelike: 0,
        })),
        "Bose exchange difference": sp.simplify(
            symmetric_integrand - symmetric_integrand.xreplace(photon_exchange)
        ),
    }


def main() -> None:
    """Print the conditional construction in deterministic LaTeX form."""
    print("Spacelike photons; Lambda > 2 M_pi; masses GeV; F and rho GeV^-1.")
    print("Real dimensionless source coefficients; source P waves include 1/f_a once.")
    print("Isospin limit; two-pion P wave; sharp common cutoff; leading QED.")
    print("The literature prefactor requires ALP P-wave sign/normalization matching.")
    print("Only the double-spectral part is constructed; single-spectral and subtraction terms remain.")
    print("No physical P waves, pion form factor, numerical TFF, or decay width are supplied.")
    for label, expression in derive_spectral_construction().items():
        print(f"{label}:\n{sp.latex(expression)}\n")


if __name__ == "__main__":
    main()
