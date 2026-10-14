"""Exact independent checks of source-gauge identities and trace-jet IBP.

Run .venv/bin/python scripts/verify_odd_source_gauge.py.
"""
from pathlib import Path
import json
import hashlib
import sympy as sp
import reduce_odd_covariant_higher as r

ROOT=Path(__file__).resolve().parents[1]

def zero(expression):
    return sp.expand(expression).doit()==sp.ZeroMatrix(3,3)

def main():
    A=r.jet('A');B=r.jet('B');C=r.jet('C')
    # Frechet derivative and commuting primitive partial jets.
    assert zero(r.partial(A*B,0)-r.partial(A,0)*B-A*r.partial(B,0))
    assert zero(r.partial(r.partial(A,0),1)-r.partial(r.partial(A,1),0))
    # A total derivative vanishes modulo integration by parts, including
    # cyclic stabilizers, products of traces and a scalar source leg.
    for expression in (sp.Trace(A*B),sp.Trace(A*A*B),sp.Trace(A*B)*sp.Trace(C),
                       r.Z*sp.Trace(A*B)):
        derivative=sp.expand(r.partial(expression,0).doit())
        assert r.normal_form(derivative)=={},expression
    for mu,nu in ((0,1),(1,2),(2,3)):
        curvature=sp.expand(-sp.I*(r.partial(r.jet(f'v{nu}'),mu)-r.partial(r.jet(f'v{mu}'),nu))
                            -(r.jet(f'v{mu}')*r.jet(f'v{nu}')-r.jet(f'v{nu}')*r.jet(f'v{mu}')))
        assert zero(curvature-(r.u(mu)*r.u(nu)-r.u(nu)*r.u(mu))/4+sp.I*r.fplus(mu,nu)/2)
        assert zero(r.covariant(r.u(nu),mu)-r.covariant(r.u(mu),nu)+r.fminus(mu,nu))
        assert zero(r.covariant(r.covariant(A,nu),mu)-r.covariant(r.covariant(A,mu),nu)
                    -(curvature*A-A*curvature))
    for mu,nu,rho in ((0,1,2),(1,2,3)):
        for fp,fm in ((r.fplus,r.fminus),(r.fminus,r.fplus)):
            residual=sp.ZeroMatrix(3,3)
            for i,j,k in ((mu,nu,rho),(nu,rho,mu),(rho,mu,nu)):
                residual+=r.covariant(fp(j,k),i)+sp.I*(r.u(i)*fm(j,k)-fm(j,k)*r.u(i))/2
            assert zero(residual)
    assert zero(r.Y()/2-r.Z*sp.Identity(3)-sum((r.b.G[i]*r.covariant(r.u(i),i) for i in range(4)),sp.ZeroMatrix(3,3)))
    # The faster source-polynomial path must reproduce the abstract matrix
    # definitions, including ordered second covariant derivatives.
    for mu,nu in ((0,1),(1,2),(2,3)):
        assert zero(r.matrix_from_polynomial(r.nc_fplus(mu,nu))-r.fplus(mu,nu))
        assert zero(r.matrix_from_polynomial(r.nc_fminus(mu,nu))-r.fminus(mu,nu))
        assert zero(r.matrix_from_polynomial(r.nc_covariant(r.nc_covariant(r.nc_V(0,1),nu),mu))
                    -r.covariant(r.covariant(r.V(0,1),nu),mu))
    assert zero(r.matrix_from_polynomial(r.nc_Y())-r.Y())
    output=ROOT/'data/vector_sector/odd_source_gauge_verification.json'
    paths=['scripts/verify_odd_source_gauge.py','scripts/reduce_odd_covariant_higher.py']
    output.write_text(json.dumps({'passed':True,'checks':['Frechet matrix product rule',
        'commuting primitive derivatives','cyclic-stabilizer and double-trace IBP',
        'scalar-source IBP','source-gauge curvature','covariant u curl',
        'commutator of covariant derivatives','both nonlinear Bianchi identities',
        'source-preserving Goldstone EOM','source-polynomial equivalence to matrix definitions'],
        'sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}},indent=2)+'\n')
    print(output.relative_to(ROOT))

if __name__=='__main__':main()
