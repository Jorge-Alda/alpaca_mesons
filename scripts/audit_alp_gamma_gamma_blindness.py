"""Project real-photon ALP resonance contractions onto stage-4 blind rows.

Run: .venv/bin/python scripts/audit_alp_gamma_gamma_blindness.py
Natural units, GeV; (+---), epsilon_0123=+1, incoming p+q1+q2=0.
Scope: tree single-trace VJP/VVP action, EM-compatible neutral ALP,
f_-=0, leading EM order. This is a sensitivity audit, not a fully normalized
ALP form-factor implementation. Independent derivative and mass spurions are
retained: X=M_a^-1(C), Y=c_q^AP+X, u=-2 da X/f_a,
chi_-=-4i a {Y,B0 mq}/f_a. Their ratio implies the stripped mass weight -2M
relative to D in both operator families. Source responses may be arbitrary.
"""
import sympy as sp
from derive_radiative_tff_amplitude import contractions, q1, q2, e1, e2, s1, s2, m, c, d, reduce


def audit():
    first = [reduce(x.subs({s1: 0, s2: 0})) for x in contractions(q1, q2, e1, e2)]
    second = [reduce(x.subs({s1: 0, s2: 0})) for x in contractions(q2, q1, e2, e1)]
    # D_j=Tr({R_j(0),Q} X), M_j=Tr({R_j(0),Q}{Y,B0mq});
    # response normalizations are absorbed in arbitrary independent weights.
    D1, D2, M1, M2, D12, M12 = sp.symbols('D_1 D_2 M_1 M_2 D_12 M_12', real=True)
    vjp = sum(c[i] * (D1 * first[i] + D2 * second[i]) for i in (0, 1, 4, 5, 6))
    vjp -= 2*c[2]*(M1*first[2] + M2*second[2])
    vvp = D12*(d[0]*(first[0]+second[0]) + d[2]*(first[4]+second[4])
                 + d[3]*(first[6]+second[6]))
    vvp -= 2*M12*d[1]*(first[2]+second[2])
    theta = c+d
    blind = sp.Matrix.hstack(sp.eye(11)[:,6], sp.eye(11)[:,10],
                            sp.Matrix([0,1,0,0,1,1,0,0,0,0,0]))
    rows = sp.Matrix([[sp.diff(expr, x) for x in theta] for expr in (vjp, vvp)])
    projection = sp.simplify(rows*blind)
    assert projection == sp.zeros(2,3)
    # Each placement vanishes separately; nonideal response weights cannot
    # spoil the cancellation or resurrect an identically zero tensor.
    for values in (first, second):
        assert values[6] == 0
        assert sp.simplify(values[1]+values[4]+values[5]) == 0
    assert first[5] == second[5] == 0
    ma = sp.Symbol('m_a', positive=True)
    return {'photon_1_on_vector_stripped_tensors':sp.Matrix([first]).subs(m, ma),
            'photon_2_on_vector_stripped_tensors':sp.Matrix([second]).subs(m, ma),
            'VJP_stripped':sp.factor(vjp).subs(m, ma), 'VVP_stripped':sp.factor(vvp).subs(m, ma),
            'blind_projection':projection}


if __name__ == '__main__':
    for label, value in audit().items():
        print(f'{label}:\n{sp.latex(value)}\n')
