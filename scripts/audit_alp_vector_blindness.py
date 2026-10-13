"""Audit existing tensor projections, not a complete physical-ALP amplitude.

Run: .venv/bin/python scripts/audit_alp_vector_blindness.py
Natural units, GeV, (+---), epsilon_0123=+1. Reuses the executable
published-operator contractions with independent ALP derivative/mass traces.
External vectors have canonical wedge wavefunctions; photon-induced tensors
have the same wedge form. No OPE relation, numerical fit, or ALP mixing input
is assumed. Overall field/Fourier normalization is immaterial to null tests.
"""
import sympy as sp
from derive_radiative_tff_amplitude import contractions, q1, q2, e1, e2, reduce


def audit():
    c = sp.symbols('c_1:8', real=True)
    d = sp.symbols('d_1:5', real=True)
    derivative_trace, mass_trace = sp.symbols('D M', real=True)
    first = contractions(q1, q2, e1, e2)
    second = contractions(q2, q1, e2, e1)
    # Relative derivative/mass building-block factors are independent of
    # flavor traces. Null tests do not require specifying the mass factor.
    direct = sum(c[i] * (mass_trace if i == 2 else derivative_trace)
                 * first[i] for i in range(7))
    two_vectors = (d[0]*derivative_trace*(first[0]+second[0])
                   + d[1]*mass_trace*(first[2]+second[2])
                   + d[2]*derivative_trace*(first[4]+second[4])
                   + d[3]*derivative_trace*(first[6]+second[6]))
    print('Stripped contractions, first assignment (O1,...,O7):')
    print(sp.latex(sp.Matrix(first)))
    print('Stripped contractions, second assignment (O1,...,O7):')
    print(sp.latex(sp.Matrix(second)))
    print('Stripped VJP polynomial, independent D and M:')
    print(sp.latex(reduce(direct)))
    print('Stripped VVP polynomial, independent D and M:')
    print(sp.latex(reduce(two_vectors)))
    parameters = c+d
    # Columns are c7, d4, common c2,c5,c6; rows follow c1..c7,d1..d4.
    nulls = sp.zeros(11,3)
    nulls[6,0] = 1
    nulls[10,1] = 1
    for index in (1,4,5):
        nulls[index,2] = 1
    rows = sp.Matrix([[sp.diff(poly,parameter) for parameter in parameters]
                      for poly in (direct,two_vectors)])
    projection = (rows*nulls).applyfunc(reduce)
    assert projection == sp.zeros(2,3)
    print('Blind-direction matrix N (11x3):')
    print(sp.latex(nulls))
    print('Stripped VJP/VVP sensitivity rows times N (2x3):')
    print(sp.latex(projection))
    identities = {
        'VJP c7 wedge projection': reduce(sp.diff(direct, c[6])),
        'VJP common c2,c5,c6 shift': reduce(sum(sp.diff(direct, c[i]) for i in (1,4,5))),
        'VVP d4 two-wedge projection': reduce(sp.diff(two_vectors, d[3])),
    }
    for label, result in identities.items():
        assert result == 0, (label, result)
        print(label + ': ' + sp.latex(result))
    # The transverse-frame scalar coefficient suffices for parity-odd
    # spin-zero -> two on-shell spin-one states, whose invariant is unique.
    # Also repeat the explicit c7/d4 zero with arbitrary polarization entries;
    # this establishes these operator vanishings before transversality.
    ev = sp.Matrix(sp.symbols('ev0:4'))
    er = sp.Matrix(sp.symbols('er0:4'))
    raw_first = contractions(q1, q2, ev, er, normalize=False)
    raw_second = contractions(q2, q1, er, ev, normalize=False)
    assert raw_first[6] == 0 and raw_second[6] == 0
    assert reduce(raw_first[1]+raw_first[4]+raw_first[5]) == 0
    assert reduce(raw_second[1]+raw_second[4]+raw_second[5]) == 0
    print('all three blind directions, arbitrary-polarization wedges: 0')
    return {'identities': identities, 'rows': rows, 'parameters': parameters,
            'nulls': nulls}


if __name__ == '__main__':
    audit()
