#!/usr/bin/env python3
"""LO only-up radial vertices; natural units. Run: .venv/bin/python scripts/derive_murayama_sigma_pp.py."""
import sympy as sp


def cartesian_vacuum():
    """Derive the vacuum condition and radial mass in natural units."""
    B_0, m_u, f_pi, f_phi, lam = sp.symbols('B_0 m_u f_pi f_varphi lambda', positive=True)
    mass_phi_squared = sp.symbols('m_varphi_squared', real=True)
    radial, angular = sp.symbols('s p', real=True)
    potential = (mass_phi_squared*((f_phi+radial)**2+angular**2)/2
                 +lam*((f_phi+radial)**2+angular**2)**2/4
                 -B_0*m_u*f_pi**2*(1+radial/f_phi))
    vacuum_equation = sp.diff(potential, radial).subs({radial: 0, angular: 0})
    vacuum_mass = sp.solve(vacuum_equation, mass_phi_squared)[0]
    radial_mass = sp.simplify(sp.diff(potential, radial, 2).subs(
        {radial: 0, angular: 0, mass_phi_squared: vacuum_mass}))
    angular_mass = sp.simplify(sp.diff(potential, angular, 2).subs(
        {radial: 0, angular: 0, mass_phi_squared: vacuum_mass}))
    radial_angular_vertex = sp.simplify(sp.diff(potential, radial, angular, angular).subs(
        {radial: 0, angular: 0, mass_phi_squared: vacuum_mass}))
    return potential, vacuum_equation, vacuum_mass, radial_mass, angular_mass, radial_angular_vertex


def calculate():
    B_0, m_u, f_pi, f_phi, lam = sp.symbols('B_0 m_u f_pi f_varphi lambda', positive=True)
    sigma, a, pi0, eta8, eta0, pip, pim, kp, km = sp.symbols('sigma a pi_0 eta_8 eta_0 pi_plus pi_minus K_plus K_minus', real=True)
    i = sp.I
    mesons = sp.Matrix([[pi0/sp.sqrt(2)+eta8/sp.sqrt(6)+eta0/sp.sqrt(3), pip, kp],
                        [pim, -pi0/sp.sqrt(2)+eta8/sp.sqrt(6)+eta0/sp.sqrt(3), 0],
                        [km, 0, -2*eta8/sp.sqrt(6)+eta0/sp.sqrt(3)]])
    # U = exp(i sqrt(2) Pi/f_pi), through second order in meson fields.
    U = sp.eye(3)+i*sp.sqrt(2)*mesons/f_pi-mesons**2/f_pi**2
    # Only the up mass depends on the Cartesian PQ field.
    L_mass = sp.expand(f_pi**2*B_0*m_u/2*((1+(sigma+i*a)/f_phi)*sp.conjugate(U[0,0])+(1+(sigma-i*a)/f_phi)*U[0,0]))
    L_potential = -lam/4*((f_phi+sigma)**2+a**2)**2
    L_cubic = sp.expand(L_mass+L_potential).coeff(sigma, 1)
    neutral = [pi0, eta8, eta0, a]
    mass_vertex = sp.Matrix(4, 4, lambda j,k: sp.simplify(-sp.diff(L_cubic,neutral[j],neutral[k])))
    charged_vertex = sp.simplify(-sp.diff(L_cubic,pip,pim))
    kaon_vertex = sp.simplify(-sp.diff(L_cubic,kp,km))
    return (B_0,m_u,f_pi,f_phi,lam), L_cubic, mass_vertex, charged_vertex, kaon_vertex


def polar_check():
    B_0,m_u,f_pi,f_phi,lam = sp.symbols('B_0 m_u f_pi f_varphi lambda', positive=True)
    m_sigma2, m_i2, m_j2 = sp.symbols('m_sigma_squared m_i_squared m_j_squared', real=True)
    v_ai,v_aj,c_i,c_j = sp.symbols('v_ai v_aj c_i c_j', real=True)
    # Cubic polar vertex, with 2 p_i.p_j = m_sigma² - m_i² - m_j².
    polar = B_0*m_u*f_pi**2/f_phi*(v_ai/f_phi-c_i/f_pi)*(v_aj/f_phi-c_j/f_pi)+(m_sigma2-m_i2-m_j2)*v_ai*v_aj/f_phi
    cartesian = 2*lam*f_phi*v_ai*v_aj+B_0*m_u*c_i*c_j/f_phi
    # The angular-field row of the pole equation determines each c_i.
    c_i_pole = (B_0*m_u*f_pi**2/f_phi**2-m_i2)*f_phi*v_ai/(B_0*m_u*f_pi)
    c_j_pole = (B_0*m_u*f_pi**2/f_phi**2-m_j2)*f_phi*v_aj/(B_0*m_u*f_pi)
    difference = sp.factor((polar-cartesian).subs({c_i:c_i_pole,c_j:c_j_pole,m_sigma2:B_0*m_u*f_pi**2/f_phi**2+2*lam*f_phi**2}))
    return polar, cartesian, difference


def show(label, value):
    print(f'{label}:\n{sp.latex(value)}\n')


if __name__ == '__main__':
    parameters, cubic, neutral, charged, kaon = calculate()
    polar, cartesian, difference = polar_check()
    potential, vacuum_equation, vacuum_mass, radial_mass, angular_mass, radial_angular_vertex = cartesian_vacuum()
    print('Natural units; LO chiral theory, tree level, only-up charges; on-shell pole equations for the comparison.')
    show('Cartesian scalar potential including the vacuum quark source', potential)
    show('Radial vacuum equation', vacuum_equation)
    show('Vacuum value of m_varphi squared', vacuum_mass)
    show('Radial mass squared', radial_mass)
    show('Angular diagonal mass squared', angular_mass)
    show('Radial angular angular potential derivative', radial_angular_vertex)
    show('Coefficient of sigma in the Cartesian Lagrangian', cubic)
    show('Neutral vertex matrix in (pi_0, eta_8, eta_0, a)', neutral)
    show('Charged pion vertex', charged)
    show('Charged kaon vertex', kaon)
    show('Projected polar vertex', polar)
    show('Projected Cartesian vertex', cartesian)
    show('Polar minus Cartesian on shell', difference)
