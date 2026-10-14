"""Derive LO diagonal ALP source profiles and parentwise LSZ assembly.

Invocation: .venv/bin/python scripts/derive_odd_alp_lsz_assembly.py
GeV units; tree level, linear in 1/fa. The photon coefficient convention
is inherited as an explicit assumption, with its signed dual/action bridge
derived separately; neither calculation is a gauged-WZW matching result.
"""
from pathlib import Path
import contextlib
import hashlib
import json
import sympy as sp
import derive_odd_power_counted_channels as d

ROOT = Path(__file__).resolve().parents[1]


def main():
    chi = tuple(d.CHI.diagonal())
    axial = sp.symbols('cA_u cA_d cA_s', real=True)
    pseudo = sp.symbols('cP_u cP_d cP_s', real=True)
    cg, g, trace_s = sp.symbols('cG g T_S', real=True)
    cga = cg + sum(axial)
    ap = tuple(a-p for a,p in zip(axial,pseudo))
    unknowns = sp.symbols('S_u S_d S_s', real=True)
    equations = [(x-d.w)*v+g*trace_s+x*c+g*cga
                 for x,v,c in zip(chi,unknowns,ap)]
    individual = sp.solve(equations,unknowns)
    trace_solution = sp.solve(trace_s-sum(individual.values()),trace_s)[0]
    profile = sp.diag(*(sp.factor(individual[v].subs(trace_s,trace_solution)) for v in unknowns))
    A,P = sp.diag(*axial),sp.diag(*pseudo)
    C = -sp.diag(*(x*c for x,c in zip(chi,ap)))-g*cga*sp.eye(3)
    ma_image = (d.CHI*profile+profile*d.CHI)/2+g*sp.trace(profile)*sp.eye(3)-d.w*profile
    assert all(sp.cancel(x)==0 for x in ma_image-C)
    # Derive the physical sources by adding the direct source and the LSZ
    # meson replacement Phi_R/a=sqrt(2)*F*(cA+S)/fa.
    mixing = A+profile
    U_direct,U_meson = 2*A/d.fa,-2*mixing/d.fa
    Y_direct = -2*(P*d.CHI+d.CHI*P)/d.fa
    Y_meson = 2*(mixing*d.CHI+d.CHI*mixing)/d.fa
    U,Y = U_direct+U_meson,Y_direct+Y_meson
    Z_direct = -2*g*cg/d.fa
    Z_meson = -2*g*sp.trace(mixing)/d.fa
    Z = sp.factor(Z_direct+Z_meson)
    defect = Y+2*d.w*U-2*Z*sp.eye(3)
    assert all(sp.cancel(x)==0 for x in defect)
    # Use algebraically independent source entries for the expensive parent
    # identity; linearity then applies to the exact rational profile above.
    a,p,ss = [sp.diag(*sp.symbols(prefix+'_u '+prefix+'_d '+prefix+'_s'))
              for prefix in ('a','p','S')]
    ud,um = 2*a/d.fa,-2*(a+ss)/d.fa
    yd,ym = -2*(p*d.CHI+d.CHI*p)/d.fa,2*((a+ss)*d.CHI+d.CHI*(a+ss))/d.fa
    probes = [(d.RHO,d.OMEGA),(d.RHO,d.Q),(d.PHI,d.Q),(d.Q,d.Q)]
    checks = 0
    for spec in d.specifications():
        for left,right in probes:
            split = d.kernel(spec,left,right,ud,yd)+d.kernel(spec,left,right,um,ym)
            assembled = d.kernel(spec,left,right,ud+um,yd+ym)
            assert sp.expand(split-assembled)==0,(spec['id'],left,right)
            checks += 1
    nc,cgamma = sp.symbols('N_c c_gamma',real=True)
    seed = cgamma+2*nc*sp.trace(d.Q**2*mixing)
    kappa = sp.diag(*sp.symbols('kappa_u kappa_d kappa_s'))
    rotated = (cgamma-2*nc*sp.trace(d.Q**2*kappa)
               +2*nc*sp.trace(d.Q**2*(mixing+kappa)))
    assert sp.expand(rotated-seed)==0
    no_mix = sp.expand((cgamma-2*nc*sp.trace(d.Q**2*kappa)).subs(
        {kappa[i,i]:-mixing[i,i] for i in range(3)}))
    assert sp.expand(no_mix-seed)==0
    # Denominator of the LO neutral mass resolvent. The intermediate chi_i-w
    # denominators in the trace solution are removable unless this vanishes.
    mass_operator = d.CHI+g*sp.ones(3)-d.w*sp.eye(3)
    determinant = sp.factor(mass_operator.det())
    assert all(sp.cancel(v*determinant).is_polynomial() for v in profile.diagonal())
    paths = ['scripts/derive_odd_alp_lsz_assembly.py',
             'scripts/derive_odd_power_counted_channels.py',
             'memory/extracts/Alda:2026nrd/appendix-a-expansion.md',
             'memory/alp/photon/alp-photon-invariant-from-mixing.md']
    out = ROOT/'data/vector_sector/odd_alp_lsz_assembly.json'
    out.write_text(json.dumps(dict(scope='LO neutral diagonal UV source; one ALP; tree',
        g_definition='m0_squared/3',trace_profile=str(trace_solution),
        profile=[str(v) for v in profile.diagonal()],resolvent_determinant=str(determinant),
        U=[str(sp.factor(v)) for v in U.diagonal()],
        Y=[str(sp.factor(v)) for v in Y.diagonal()],Z=str(Z),
        photon_seed=str(seed),photon_status='conditional inherited anomaly convention; signed dual/action conversion derived in odd_photon_phase_matching.json',
        parent_LSZ_checks=checks,source_EOM_defect_zero=True,
        pole_domain='det(Mass-w I) nonzero and perturbative ALP-meson overlaps small',
        sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}),indent=2)+'\n')
    transcript = ROOT/'data/vector_sector/odd_alp_lsz_assembly.md'
    with transcript.open('w') as handle,contextlib.redirect_stdout(handle):
        print('# LO ALP source and LSZ assembly\n')
        for label,value in [('LO source C',C),('Component equations',sp.Matrix(equations)),
            ('Solved trace of S',trace_solution),('S = Ma inverse C',profile),
            ('LO neutral resolvent determinant',determinant),
            ('Meson replacement Phi_R/a',sp.sqrt(2)*d.F*mixing/d.fa),
            ('Direct U',U_direct),('Meson U',U_meson),('Physical U',U),
            ('Direct Y',Y_direct),('Meson Y',Y_meson),('Physical Y',Y),
            ('Topological contact Z',Z),('Source EOM defect',sp.zeros(3)),
            ('Conditional invariant photon coefficient',seed)]:d.show(label,value)
        print(f'\nAll {checks} parent/probe LSZ identities passed.\n')
        print('Photon seed uses the stated anomaly assumption; its signed dual/action conversion is derived separately in odd_photon_phase_matching.json.\n')
    print(out.relative_to(ROOT));print(transcript.relative_to(ROOT))


if __name__=='__main__':main()
