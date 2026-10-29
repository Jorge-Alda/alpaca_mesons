#!/usr/bin/env python
"""Exact coefficient ranks of GomezDumm:2010ra KKpi vector form factors.

Run: .venv/bin/python scripts/audit_kkpi_coefficient_rank.py
Inputs: Appendix A C^R, C^RR, D^R, E^R and Eqs. (f31),(f32)
in taukkp_v2.tex, arXiv:0911.2640 source archive. Couplings are
dimensionless, invariants and squared masses have units GeV^2. Tree level,
exact strong isospin, ideal omega/phi mixing, nonzero fixed even couplings,
and freely varied odd couplings are assumed. No short-distance relations
are imposed. Canonical epsilon/coupling conversion reverses all odd
couplings together and preserves every rank and nullspace here.

Polynomial coefficient rows correspond to a continuum of kinematics,
not to statistically independent measurements. Zero-width denominators
are used only as rational functions for this algebraic audit.
"""

import sympy as sp


def coefficient_matrix(expressions, couplings, invariants):
    """Extract coefficient rows after multiplying each common denominator."""
    rows = []
    for expression in expressions:
        numerator = sp.cancel(expression).as_numer_denom()[0]
        polynomial = sp.Poly(numerator, *invariants)
        for coefficient in polynomial.coeffs():
            rows.append([sp.diff(coefficient, coupling) for coupling in couplings])
    return sp.Matrix(rows)


def derive_maps():
    """Construct the source polynomials and their exact coefficient maps."""
    c = sp.symbols('c_1:8')
    d = sp.symbols('d_1:5')
    g = sp.symbols('g_1:6')
    q, s, t = sp.symbols('Q2 s t')
    kaon, pion = sp.symbols('m_K2 m_pi2', nonzero=True)
    rho, omega, kstar = sp.symbols('M_rho2 M_omega2 M_Kstar2', nonzero=True)
    vector_source, vector_strong, vector_mass = sp.symbols('F_V G_V M_V', nonzero=True)
    u = q-s-t+2*kaon+pion

    def direct(x, mass1, mass2, mass3):
        return ((c[0]-c[1]+c[4])*q
                -(c[0]-c[1]-c[4]+2*c[5])*x
                +(c[0]+c[1]+8*c[2]-c[4])*mass3
                +8*c[3]*(mass1-mass2))

    def double(x, mass):
        return d[2]*(q+x)+(d[0]+8*d[1]-d[2])*mass

    local_charged = ((g[0]+2*g[1]-g[2])*(s+t)
                     -2*g[1]*(q+kaon)-(g[0]-g[2])*(3*kaon+pion)
                     +2*g[3]*(kaon+pion)+2*g[4]*kaon)
    local_neutral = (g[0]+2*g[1]-g[2])*(t-u)
    isoscalar = [g[0]-g[1]-g[2], 2*g[3]+g[4]]
    # Common normalizations -4 G_V/(M_V F^3) and
    # 2 sqrt(2) G_V/(M_V F^3), respectively, are removed.
    charged = (direct(s,kaon,kaon,pion)/(omega-s)
               +direct(t,kaon,pion,kaon)/(kstar-t)
               -2*vector_source/vector_strong*local_charged/(rho-q)
               -sp.sqrt(2)*vector_source*vector_mass/(rho-q)
               *(double(s,pion)/(omega-s)+double(t,kaon)/(kstar-t)))
    neutral = (direct(t,kaon,pion,kaon)/(kstar-t)
               -direct(u,kaon,pion,kaon)/(kstar-u)
               -2*vector_source/vector_strong*local_neutral/(rho-q)
               -sp.sqrt(2)*vector_source*vector_mass/(rho-q)
               *(double(t,kaon)/(kstar-t)-double(u,kaon)/(kstar-u)))
    # A symbolic compensation between VJP and the VPPP mass term.
    # This explicitly proves the rank loss seen in the full amplitudes.
    g4_displacement, g5_displacement = sp.symbols('delta_g4 delta_g5')
    g1_displacement = -2*((kaon+pion)*g4_displacement+kaon*g5_displacement)/(rho-2*kaon-pion)
    x_coefficient_displacement = vector_source/vector_strong*g1_displacement
    mass_coefficient_displacement = -x_coefficient_displacement*omega/pion
    compensation = {
        c[0]:mass_coefficient_displacement/2,
        c[1]:(mass_coefficient_displacement+x_coefficient_displacement)/2,
        c[3]:(-x_coefficient_displacement*kstar-mass_coefficient_displacement*kaon)/(8*(kaon-pion)),
        c[4]:x_coefficient_displacement/2,
        g[0]:g1_displacement,g[1]:-g1_displacement/2,
        g[3]:g4_displacement,g[4]:g5_displacement,
    }
    compensation_responses = [sp.factor(sum(sp.diff(amplitude,coupling)*change
                                          for coupling,change in compensation.items()))
                              for amplitude in (charged,neutral)]
    invariants = (q,s,t)
    maps = {
        'VPPP neutral': coefficient_matrix([local_neutral],g,invariants),
        'VPPP charged': coefficient_matrix([local_charged],g,invariants),
        'VPPP both': coefficient_matrix([local_neutral,local_charged],g,invariants),
        'VPPP isoscalar': coefficient_matrix(isoscalar,g,invariants),
        'VPPP neutral plus isoscalar': coefficient_matrix([local_neutral,*isoscalar],g,invariants),
        'VPPP charged plus isoscalar': coefficient_matrix([local_charged,*isoscalar],g,invariants),
        'VPPP both plus isoscalar': coefficient_matrix([local_neutral,local_charged,*isoscalar],g,invariants),
        'VJP numerator pair': coefficient_matrix([direct(s,kaon,kaon,pion),direct(t,kaon,pion,kaon)],c,invariants),
        'VVP numerator pair': coefficient_matrix([double(s,pion),double(t,kaon)],d,invariants),
    }
    # These substitutions select generic fixed inputs for a finite exact
    # certificate. They are not physical masses or a phenomenological fit.
    # Symbolic stripped maps above retain arbitrary meson masses.
    exact_inputs = {kaon:sp.Rational(1,4),pion:sp.Rational(1,50),
                    rho:sp.Rational(3,5),omega:sp.Rational(61,100),
                    kstar:sp.Rational(4,5),vector_source:sp.Rational(3,20),
                    vector_strong:sp.Rational(7,100),vector_mass:sp.Rational(77,100)}
    odd = (*c,*d,*g)
    preliminary = [c[0]+c[1]-c[4],-c[0]+c[1]+c[4]-2*c[5],c[2],c[3],d[0],d[1],d[2]]
    preliminary_matrix = coefficient_matrix(preliminary,odd,invariants)
    maps['preliminary coefficient span'] = preliminary_matrix
    for label, amplitudes in [('full charged',[charged]),('full neutral',[neutral]),('full both',[charged,neutral])]:
        maps[label] = coefficient_matrix([a.subs(exact_inputs) for a in amplitudes],odd,invariants)
        maps[label+' plus isoscalar'] = maps[label].col_join(coefficient_matrix(isoscalar,odd,invariants))
        maps[label+' plus preliminary'] = maps[label].col_join(preliminary_matrix)
        maps[label+' plus preliminary plus isoscalar'] = maps[label+' plus isoscalar'].col_join(preliminary_matrix)
    return dict(maps=maps,c=c,d=d,g=g,odd=odd,invariants=invariants,
                kaon=kaon,pion=pion,local_charged=local_charged,
                local_neutral=local_neutral,isoscalar=isoscalar,
                direct_pion=direct(s,kaon,kaon,pion),exact_inputs=exact_inputs,
                compensation=compensation,compensation_responses=compensation_responses)


def main():
    results = derive_maps()
    print(__doc__)
    for label in ('local_charged','local_neutral','direct_pion'):
        print(label+': '+sp.latex(results[label]))
    for label, matrix in results['maps'].items():
        print(label+': rank '+str(matrix.rank()))
        print('nullspace: '+sp.latex(matrix.nullspace()))
    print('Finite exact certificate inputs: '+sp.latex(results['exact_inputs']))
    print('Symbolic VJP/VPPP compensation: '+sp.latex(results['compensation']))
    print('Compensation responses (charged, neutral): '+sp.latex(results['compensation_responses']))


if __name__ == '__main__':
    main()
