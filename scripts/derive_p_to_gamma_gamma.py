"""Real-photon limit, identical-particle width and OPE sensitivity.

Run: .venv/bin/python scripts/derive_p_to_gamma_gamma.py
Tree level; common vector pole; natural units (GeV); conventions inherited
from derive_radiative_tff_amplitude. No physical eta inputs are assumed.
"""
import sympy as sp
from derive_radiative_tff_amplitude import derive, reduce, show, m, s1, s2, F, FV, MV, NC, FQ, FH, c, d


def main():
    _, _, wzw, direct, exchange = derive()
    real = [reduce(x.subs({s1: 0, s2: 0})) for x in (wzw, direct, exchange)]
    for name, value in zip(('WZW', 'VJP', 'VVP'), real):
        show(name + ' real-photon form factor', value)
    amplitude = sum(real)
    # Rest frame, two transverse polarization bases; explicit epsilon sums.
    q1, q2 = sp.Matrix([m/2, 0, 0, m/2]), sp.Matrix([m/2, 0, 0, -m/2])
    polarizations = [sp.Matrix([0, 1, 0, 0]), sp.Matrix([0, 0, 1, 0])]
    tensor_sum = sum(sum(sp.LeviCivita(a,b,r,s)*u[a]*v[b]*q1[r]*q2[s]
                         for a in range(4) for b in range(4)
                         for r in range(4) for s in range(4))**2
                     for u in polarizations for v in polarizations)
    show('Summed squared Lorentz tensor', sp.simplify(tensor_sum))
    alpha, absF = sp.symbols('alpha_em abs_F_P', positive=True)
    # Integrated massless two-body phase space=1/(8pi); initial flux=1/(2m).
    width = sp.simplify((4*sp.pi*alpha)**2 * absF**2 * tensor_sum
                        / (2*m) / (8*sp.pi) / sp.factorial(2))
    show('Width including identical-photon factor', width)
    constraints = {c[2]:-c[0]/4, c[1]:c[0]+c[4],
        c[5]:c[4]-NC*MV/(64*sp.pi**2*sp.sqrt(2)*FV),
        d[0]:-8*d[1]-NC*MV**2/(64*sp.pi**2*FV**2)+F**2/(4*FV**2),
        d[2]:-NC*MV**2/(64*sp.pi**2*FV**2)+F**2/(8*FV**2)}
    constrained = reduce(amplitude.subs(constraints, simultaneous=True))
    show('OPE-constrained real-photon form factor, independent spurion', constrained)
    equal = reduce(constrained.subs(FH, 2*m**2*FQ))
    show('Equal-spurion real-photon form factor', equal)
    show('Independent-spurion remainder', sp.factor(constrained-equal))
    pi = reduce(equal.subs({FQ:sp.sqrt(2)/3, NC:3}))
    show('Pion width', sp.factor(width.subs(absF, pi)))
    parameters = c+d
    ope = sp.Matrix([[sp.diff(x,y) for y in parameters] for x in
        (c[0]+4*c[2], c[0]-c[1]+c[4], c[4]-c[5], d[0]+8*d[1], d[2])])
    row = sp.Matrix([[sp.diff(amplitude,y) for y in parameters]])
    projection = sp.simplify(row*sp.Matrix.hstack(*ope.nullspace()))
    show('Real-photon sensitivity along OPE null directions', projection)
    # Across channels only the common coefficient times FH-2m^2 FQ varies.
    normalized = projection.applyfunc(lambda x: sp.cancel(x/(FH-2*m**2*FQ)))
    assert not normalized.has(m, FQ, FH)
    assert normalized.rank() == 1
    assert projection.subs(FH, 2*m**2*FQ) == sp.zeros(1, 6)
    assert sp.simplify(width-alpha**2*m**3*absF**2*sp.pi/4) == 0
    assert sp.simplify(pi-(1/(4*sp.pi**2*F)-F*m**2/(3*MV**4))) == 0
    show('OPE plus real-photon rank, independent spurion', sp.Integer(6))
    show('OPE plus real-photon rank, equal spurion', sp.Integer(5))


if __name__ == '__main__':
    main()
