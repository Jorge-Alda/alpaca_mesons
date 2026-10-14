"""Exact legacy-to-canonical epsilon dictionary, metric (+---), GeV units.

Run: .venv/bin/python scripts/derive_levi_civita_convention.py
Legacy tensor lower0123=+1; physical tensor upper0123=+1.
No approximation or field rephasing is made. Legacy stored scalar matrices
remain inputs in their explicitly named source coordinates.
"""
import json
from pathlib import Path
import sympy as sp
try:
    from scripts.symbolic_utils import levi_civita_contract, legacy_epsilon_determinant
except ModuleNotFoundError:
    from symbolic_utils import levi_civita_contract, legacy_epsilon_determinant


def derive():
    metric = sp.diag(1, -1, -1, -1)
    lower0123 = sp.prod(metric[i, i] for i in range(4))
    c, d, g, C = sp.symbols('c d g C', real=True)
    A, Lc, Ld = sp.symbols('A L_c L_d', real=True)
    legacy = A + Lc*c + Ld*d
    target = sp.expand(-legacy.subs({c: -c, d: -d}, simultaneous=True))
    N, F, FV, MV = sp.symbols('N_C F F_V M_V', positive=True)
    source_rhs = dict(c5_minus_c6=N*MV/(64*sp.pi**2*sp.sqrt(2)*FV),
        d1_plus_8d2=-N*MV**2/(64*sp.pi**2*FV**2)+F**2/(4*FV**2),
        d3=-N*MV**2/(64*sp.pi**2*FV**2)+F**2/(8*FV**2))
    basis = [sp.eye(4)[:, i] for i in range(4)]
    # Explicit photon fields: both Wick assignments, Fourier derivative i*q.
    q1,q2,e1,e2 = (sp.Matrix(x) for x in ((2,0,0,2),(2,0,0,-2),(0,1,0,0),(0,0,1,0)))
    field = [sp.I*(q*e.T-e*q.T) for q,e in ((q1,e1),(q2,e2))]
    lowered = [metric*f*metric for f in field]
    ft = sum(sp.LeviCivita(mu,nu,rho,sigma)*lowered[0][mu,nu]*lowered[1][rho,sigma]
             for mu in range(4) for nu in range(4) for rho in range(4) for sigma in range(4))
    E = levi_civita_contract(q1,e1,q2,e2)
    derivative = 2*sp.I**2*E
    fa = sp.Symbol('f_a', positive=True)
    dual_coefficient_over_e2 = C/(16*sp.pi**2*fa)
    contact_over_e2 = sp.simplify(dual_coefficient_over_e2*ft/E)
    return dict(metric=metric, epsilon_lower0123=lower0123,
        legacy_basis_contraction=legacy_epsilon_determinant(*basis),
        canonical_basis_contraction=levi_civita_contract(*basis),
        coupling_dictionary=sp.Dict({c:-c,d:-d,g:-g}),
        legacy_scalar=legacy, canonical_scalar=target,
        canonical_OPE_rhs=sp.Dict({sp.Symbol(k):-v for k,v in source_rhs.items()}),
        canonical_dual_coefficient_over_e2=dual_coefficient_over_e2,
        canonical_contact_over_e2=contact_over_e2,
        photon_dual_relative_orientation=sp.Integer(1),
        canonical_Fdual_over_E=sp.simplify(ft/E),
        canonical_epsilon_dAdA_over_E=sp.simplify(derivative/E),
        action_invariance=sp.expand((-g)*(-sp.Symbol('O'))-g*sp.Symbol('O')),
        total_squared_invariance=sp.expand(target**2-legacy.subs({c:-c,d:-d}, simultaneous=True)**2))


def main():
    result = derive()
    root = Path(__file__).resolve().parents[1]
    path = root/'data/vector_sector/levi_civita_convention.json'
    path.write_text(json.dumps({k:str(v) for k,v in result.items()},indent=2)+'\n')
    transcript = root/'data/vector_sector/levi_civita_convention.md'
    transcript.write_text('# Canonical Levi-Civita dictionary\n\n'+''.join(
        f'{k}:\n\n$${sp.latex(v)}$$\n\n' for k,v in result.items()))
    print(path.relative_to(root))
    print(transcript.relative_to(root))


if __name__ == '__main__':
    main()
