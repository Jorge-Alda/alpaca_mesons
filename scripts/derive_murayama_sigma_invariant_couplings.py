#!/usr/bin/env python3
"""Only-up invariant matching and radial source, LO, natural GeV units.

Run: .venv/bin/python scripts/derive_murayama_sigma_invariant_couplings.py
"""
import sympy as sp


def derive():
    mu, md, ms, f_a, sigma, x = sp.symbols('m_u m_d m_s f_a sigma x', real=True, nonzero=True)
    B0, f_pi, quartic, mass_i, mass_j = sp.symbols('B_0 f_pi lambda m_i m_j', real=True, nonzero=True)
    c_ap = sp.symbols('c_u^{AP}', real=True, nonzero=True)
    radial_coupling = sp.symbols('y_u', real=True)
    up_projector = sp.diag(1, 0, 0)
    quark_mass = sp.diag(mu, md, ms)
    mass_source = sp.diag(mu*(1+sigma/f_a)*sp.exp(-sp.I*x), md, ms)
    invariant_pseudoscalar = sp.simplify((mass_source.diff(x).subs({x: 0, sigma: 0}))/sp.I)
    invariant_coupling = sp.diag(sp.Rational(1, 2), 0, 0)
    radial_source = sp.simplify(mass_source.diff(sigma).subs({x: 0, sigma: 0}))
    mixed_source = sp.simplify(mass_source.diff(sigma, x).subs({x: 0, sigma: 0}))
    model_axion_mass = B0*mu*f_pi**2/f_a**2
    model_bracket = 2*quartic*f_a + f_a/(B0*mu*f_pi**2)*(model_axion_mass-mass_i**2)*(model_axion_mass-mass_j**2)
    invariant_bracket = sp.factor(model_bracket.subs(mu, 2*mu*c_ap))
    general_mass_source = sp.diag((mu+radial_coupling*sigma)*sp.exp(-sp.I*x), md, ms)
    general_phase_source = sp.simplify(general_mass_source.diff(x).subs({x: 0, sigma: 0}))
    general_radial_source = sp.simplify(general_mass_source.diff(sigma).subs({x: 0, sigma: 0}))
    return locals()


def main():
    result = derive()
    for label, key in [
        ('Only-up mass source', 'mass_source'), ('Invariant c_AP', 'invariant_coupling'),
        ('First phase derivative divided by i', 'invariant_pseudoscalar'),
        ('Radial derivative', 'radial_source'), ('Mixed radial and phase derivative', 'mixed_source'),
        ('Original coupling factor', 'model_bracket'), ('Coupling factor in invariant c_AP', 'invariant_bracket'),
        ('Phase derivative with independent radial coupling', 'general_phase_source'),
        ('Radial derivative with independent radial coupling', 'general_radial_source'),
    ]:
        print(f'{label}:\n{sp.latex(result[key])}\n')


if __name__ == '__main__':
    main()
