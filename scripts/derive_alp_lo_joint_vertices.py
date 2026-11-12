#!/usr/bin/env python3
"""Joint LO chiral vertices with arbitrary Taylor coefficients.

Run: .venv/bin/python scripts/derive_alp_lo_joint_vertices.py
Tests: .venv/bin/python -m unittest discover -s scripts -p test_derive_alp_lo_joint_vertices.py
All fields and decay constants use natural units; chi=2 B_0 M has mass dimension two.
Matrices are explicit flavour matrices; source arrays are derivatives at x=0,
not power-series coefficients. All momenta are incoming. A returned polynomial
uses q_i*q_j for p_i.p_j and omits the overall Feynman-rule factor i.
"""
from __future__ import annotations

import sympy as sp


def _multiply(left, right, order):
    """Multiply matrix-valued homogeneous series with a total-degree cutoff."""
    size = left[0].rows
    return [sum((left[k] * right[n-k] for k in range(n+1)), sp.zeros(size))
            for n in range(order+1)]


def _trace_product(left, right, degree):
    return sum((sp.trace(left[k] * right[degree-k]) for k in range(degree+1)), sp.S.Zero)


def vertex(meson_generators, alp_components, momentum_symbols, chi_derivatives,
           right_derivatives, left_derivatives, theta_derivatives,
           f_pi, f_a, anomaly_mass_squared, bare_alp_mass_squared):
    """Differentiate the homogeneous LO Lagrangian in labeled external legs.

    Each leg is (generator, ALP component); generators can be charged or neutral,
    and need not be Hermitian individually. q_i q_j denotes p_i dot p_j.
    The source matrices are chi derivatives and Hermitian-current derivatives.
    Theta(0)=0 and chi(0) Hermitian define the assumed aligned vacuum; an ALP
    tadpole proportional to Re tr chi'(0) must also vanish or be cancelled.
    """
    order = len(meson_generators)
    if not (len(alp_components) == len(momentum_symbols) == order):
        raise ValueError('one generator, ALP component and momentum per leg required')
    size = meson_generators[0].rows
    t = sp.symbols(f't0:{order}', real=True)
    q = [-sp.I * value for value in momentum_symbols]
    zero = sp.zeros(size)
    identity = sp.eye(size)
    phi = sum((t[i]*meson_generators[i] for i in range(order)), zero)
    dphi = sum((t[i]*q[i]*meson_generators[i] for i in range(order)), zero)
    alp = sum(t[i]*alp_components[i] for i in range(order))
    dalp = sum(t[i]*q[i]*alp_components[i] for i in range(order))
    u = [identity] + [(sp.I*phi/(2*f_pi))**n/sp.factorial(n) for n in range(1,order+1)]
    ud = [identity] + [(-sp.I*phi/(2*f_pi))**n/sp.factorial(n) for n in range(1,order+1)]
    du = [zero] + [sum(((sp.I*phi/(2*f_pi))**j*(sp.I*dphi/(2*f_pi))*
                              (sp.I*phi/(2*f_pi))**(n-j-1) for j in range(n)), zero)/sp.factorial(n)
                     for n in range(1,order+1)]
    dud = [zero] + [sum(((-sp.I*phi/(2*f_pi))**j*(-sp.I*dphi/(2*f_pi))*
                               (-sp.I*phi/(2*f_pi))**(n-j-1) for j in range(n)), zero)/sp.factorial(n)
                      for n in range(1,order+1)]
    right = [zero] + [dalp/f_a * (alp/f_a)**(n-1)/sp.factorial(n-1)*right_derivatives[n-1] for n in range(1,order+1)]
    left = [zero] + [dalp/f_a * (alp/f_a)**(n-1)/sp.factorial(n-1)*left_derivatives[n-1] for n in range(1,order+1)]
    ud_du = _multiply(ud,du,order)
    u_dud = _multiply(u,dud,order)
    ud_right_u = _multiply(_multiply(ud,right,order),u,order)
    u_left_ud = _multiply(_multiply(u,left,order),ud,order)
    u_mu = [sp.I*(ud_du[n]-u_dud[n])+ud_right_u[n]-u_left_ud[n]
            for n in range(order+1)]
    chi = [(alp/f_a)**n/sp.factorial(n)*chi_derivatives[n] for n in range(order+1)]
    chid = [(alp/f_a)**n/sp.factorial(n)*chi_derivatives[n].conjugate().T for n in range(order+1)]
    chi_plus = [a+b for a,b in zip(_multiply(_multiply(ud,chi,order),ud,order),
                                    _multiply(_multiply(u,chid,order),u,order))]
    mass = sp.trace(chi_plus[order])
    topological = -sp.trace(phi)/f_pi - 2*sum(theta_derivatives[n]*(alp/f_a)**n/sp.factorial(n) for n in range(1,order+1))
    # Extract the requested total degree before differentiating external legs.
    z = sp.Symbol('z', real=True)
    topological_degree = sp.expand((topological.subs({x:z*x for x in t}))**2).coeff(z,order)
    bare = (dalp**2-bare_alp_mass_squared*alp**2)/2 if order == 2 else 0
    lagrangian = f_pi**2/4*(_trace_product(u_mu,u_mu,order)+mass) - f_pi**2*anomaly_mass_squared/(4*size)*topological_degree + bare
    result = sp.diff(lagrangian,*t)
    return sp.expand(result.subs({x:0 for x in t}))


def quadratic_blocks(meson_generators, chi_derivatives, right_derivatives,
                     left_derivatives, theta_derivatives, f_pi, f_a,
                     anomaly_mass_squared, bare_alp_mass_squared):
    """Return (K, M²) for meson_generators followed by the ALP.

    The matrix pair can be passed directly to the spectrum solver. All Taylor
    arrays include their zeroth, first and second derivatives, even if zero.
    """
    zero = sp.zeros(meson_generators[0].rows)
    generators = list(meson_generators) + [zero]
    alp_components = [0]*len(meson_generators) + [1]
    size = len(generators)
    kinetic = sp.zeros(size)
    mass_squared = sp.zeros(size)
    for i in range(size):
        for j in range(i, size):
            common = ([generators[i], generators[j]],
                      [alp_components[i], alp_components[j]])
            parameters = (chi_derivatives, right_derivatives, left_derivatives,
                          theta_derivatives, f_pi, f_a, anomaly_mass_squared,
                          bare_alp_mass_squared)
            at_zero = vertex(*common, [0, 0], *parameters)
            at_one = vertex(*common, [1, 1], *parameters)
            kinetic[i, j] = kinetic[j, i] = sp.expand(at_zero-at_one)
            mass_squared[i, j] = mass_squared[j, i] = sp.expand(-at_zero)
    return kinetic, mass_squared


def main():
    f_pi, f_a, B_0, m_0_sq, m_a0_sq = sp.symbols('f_pi f_a B_0 m_0_sq m_a0_sq', positive=True)
    m_u, m_d, m_s, c_G, chi2 = sp.symbols('m_u m_d m_s c_G chi_2', real=True)
    q1,q2=sp.symbols('q_1 q_2', real=True)
    zero=sp.zeros(3)
    chi=[2*B_0*sp.diag(m_u,m_d,m_s),zero,sp.diag(chi2,0,0)]
    sources=[zero]*3
    out=vertex([zero,zero],[1,1],[q1,q2],chi,sources,sources,[0,c_G,0],f_pi,f_a,m_0_sq,m_a0_sq)
    print('Aligned vacuum; general mass-source second derivative; all momenta incoming; q_i q_j = p_i dot p_j.')
    print('Mass source derivatives chi(0), chi prime(0), chi double-prime(0):\n'+sp.latex(chi))
    print('Topological source through first order:\n'+sp.latex(c_G*sp.Symbol('a')/f_a))
    print('Two-ALP amputated Lagrangian vertex (omit overall i):\n'+sp.latex(out))


if __name__ == '__main__':
    main()
