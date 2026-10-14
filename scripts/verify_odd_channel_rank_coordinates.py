"""Exact mass/response coordinate identities used by the channel rank script.

Run .venv/bin/python scripts/verify_odd_channel_rank_coordinates.py.
This derives every column's scale factor before the rank coordinate change.
Natural units, neutral cubic sector and LO EOM-compatible sources only.
"""
from pathlib import Path
import hashlib
import json
import sympy as sp
import derive_odd_power_counted_channels as d

ROOT = Path(__file__).resolve().parents[1]


def main():
    A = sp.diag(*sp.symbols('A_u A_d A_s'))
    B = sp.diag(*sp.symbols('B_u B_d B_s'))
    U = sp.diag(*sp.symbols('U_u U_d U_s'))
    x,z,ratio,rd,rs,scale = sp.symbols('x z mass_ratio r_d r_s W_scale')
    W = sp.diag(1,rd,rs)
    Y = -2*d.w*U+2*d.Z*sp.eye(3)
    scaledY = -2*x*U+2*z*sp.eye(3)
    checked = 0
    for spec in d.specifications():
        if spec['mode'] in ('fminus','fmdd'):
            continue
        mass_power = 4 if spec['mode']=='Ydd' else 2
        if spec['kind']=='VVP':
            physical = d.kernel(spec,A,B,U,Y).subs({d.s:d.M1**2,d.t:d.M1**2*ratio,
                         d.h:(d.w-d.M1**2*(1+ratio))/2})
            physical = sp.expand(physical.subs({d.w:d.M1**2*x,d.Z:d.M1**2*z},simultaneous=True))
            reduced = d.kernel(spec,A,B,U,scaledY).subs({d.s:1,d.t:ratio,d.h:(x-1-ratio)/2})
            assert sp.expand(physical-d.M1**mass_power*reduced)==0
            checked += 1
        physical = d.kernel(spec,A,d.Q if spec['kind']=='VJP' else scale*W,U,Y)
        physical = physical.subs({d.s:d.M1**2,d.t:0,d.h:(d.w-d.M1**2)/2})
        physical = sp.expand(physical.subs({d.w:d.M1**2*x,d.Z:d.M1**2*z},simultaneous=True))
        reduced = d.kernel(spec,A,d.Q if spec['kind']=='VJP' else W,U,scaledY).subs(
            {d.s:1,d.t:0,d.h:(x-1)/2})
        response_power = 0 if spec['kind']=='VJP' else 1
        assert sp.expand(physical-d.M1**mass_power*scale**response_power*reduced)==0
        checked += 1
        physical = d.kernel(spec,scale*W,d.Q if spec['kind']=='VJP' else scale*W,U,Y).subs(
            {d.s:0,d.t:0,d.h:d.w/2})
        reduced = d.kernel(spec,W,d.Q if spec['kind']=='VJP' else W,U,Y).subs(
            {d.s:0,d.t:0,d.h:d.w/2})
        response_power = 1 if spec['kind']=='VJP' else 2
        assert sp.expand(physical-scale**response_power*reduced)==0
        checked += 1
    paths = ['scripts/verify_odd_channel_rank_coordinates.py','scripts/derive_odd_power_counted_channels.py']
    output = ROOT/'data/vector_sector/odd_channel_rank_coordinate_verification.json'
    output.write_text(json.dumps(dict(passed=True,exact_column_identities=checked,
        mass_scaling='m_a2=M1^2*x; Z_a=M1^2*z; M2^2=M1^2*mass_ratio',
        response_scaling='W(0)=W_scale*diag(1,r_d,r_s)',
        scope='invertible source-family coordinates and nonzero column scales on the stated generic domain',
        sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}),indent=2)+'\n')
    transcript = ROOT/'data/vector_sector/odd_channel_rank_coordinate_verification.md'
    transcript.write_text('# Exact rank-coordinate identities\n\n'
        'All scalar sources remain LO EOM-compatible. The physical gradient normalization is preserved separately.\n\n'
        '$$\n'+sp.latex(sp.Matrix([checked]))+'\n$$\n\n'
        'Each physical column equals its ratio-coordinate column times the derived nonzero mass/response power.\n'
        'These invertible row-coordinate and column-scale changes preserve the retained-column ranks.\n')
    print(output.relative_to(ROOT));print(transcript.relative_to(ROOT))


if __name__=='__main__':
    main()
