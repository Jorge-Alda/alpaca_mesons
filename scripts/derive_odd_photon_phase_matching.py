"""Match the conditional photon dual to the bounded odd-action tensor phase.

Run: .venv/bin/python scripts/derive_odd_photon_phase_matching.py
Natural units (GeV), tree level and one ALP at O(1/fa). The common source
tensor has epsilon_lower0123=+1 and metric (+---). The old photon record
specified a dual formula but did not specify its orientation. This executable
derives the orientation required by its assumed positive WZW coefficient.
It does not derive a UV anomaly coefficient or a gauged-WZW counterterm.
"""
from pathlib import Path
from itertools import permutations
import contextlib
import hashlib
import json
import sympy as sp
import derive_radiative_tff_amplitude as sm
import derive_odd_power_counted_channels as pc

ROOT = Path(__file__).resolve().parents[1]
METRIC = sp.diag(1, -1, -1, -1)
EPS = [(indices, sp.LeviCivita(*indices))
       for indices in permutations(range(4))]


def photon_tensors(q1, q2, e1, e2, eta):
    """Explicit sums; SymPy has no Wick-assignment operation for these fields.

    eta compares the photon-dual volume form to the metric-raised common
    source volume form. Both gauge fields carry Fourier derivatives i*q.
    """
    fields = [sp.I * sm.wedge(q, e) for q, e in ((q1, e1), (q2, e2))]
    lower = [METRIC * field * METRIC for field in fields]
    duals = [sp.Matrix(4, 4, lambda mu, nu:
        eta * sum(sp.prod(METRIC[i, i] for i in indices) * sign
                  * field[indices[2], indices[3]] / 2
                  for indices, sign in EPS
                  if indices[0] == mu and indices[1] == nu))
        for field in lower]
    ft = sum(lower[0][mu, nu] * duals[1][mu, nu]
             + lower[1][mu, nu] * duals[0][mu, nu]
             for mu in range(4) for nu in range(4))
    derivatives = sum(sign * (
        sp.I*q1[mu]*e1[nu]*sp.I*q2[rho]*e2[sigma]
        + sp.I*q2[mu]*e2[nu]*sp.I*q1[rho]*e1[sigma])
        for (mu, nu, rho, sigma), sign in EPS)
    E = sum(sign*q1[mu]*e1[nu]*q2[rho]*e2[sigma]
            for (mu, nu, rho, sigma), sign in EPS)
    R = sm.reference(q1, q2, e1, e2)
    return tuple(sp.expand(x) for x in (ft, derivatives, E, R))


def resonance_gradients(U, Y, W1, W2):
    """Bounded neutral cubic action, coefficient of i*e**2*E (not of R).

    The supplied W matrices are the signed canonical spectral responses;
    two explicit photon assignments and the exported VVP Bose factor remain.
    Couplings multiply raw parent columns, which are redundant off shell.
    """
    result = {}
    # Swap the local momentum polynomial before inserting the already-labelled
    # W2(q2^2). Swapping a supplied response denominator would place it on the
    # wrong photon line. Neutral diagonal responses are the declared scope.
    slots = sp.diag(*(sp.Dummy(f'W2_slot_{i}') for i in range(3)))
    response_substitutions = {slots[i, i]: W2[i, i] for i in range(3)}
    for parent in pc.specifications():
        norm = pc.degree_normalization(parent)
        if parent['kind'] == 'VJP':
            value = -2*sp.sqrt(2)*norm * (
                pc.kernel(parent, W1, pc.Q, U, Y)
                + pc.kernel(parent, slots, pc.Q, U, Y).xreplace(
                    {pc.s: pc.t, pc.t: pc.s}).xreplace(response_substitutions))
        else:
            value = 2*norm*pc.kernel(parent, W1, W2, U, Y)
        result[parent['id']] = sp.expand(value)
    return result


def original_coefficients():
    """Read the already-derived original-coupling embedding, not fitted data."""
    artifact = json.loads((ROOT/'data/vector_sector/odd_leading_action_count.json').read_text())
    f3 = pc.b.trace_catalog(('A', 'B', 'U'))
    single = next(i for i, family in enumerate(f3)
                  if len(next(iter(family))) == 1)
    symmetric = pc.b.flavor_exchange_basis(f3, 1)
    vsingle = next(i for i, family in enumerate(symmetric)
                   if len(next(iter(family))) == 1)
    result = {}
    for coefficient, coordinates in zip(
            (sm.c[0], sm.c[1], sm.c[4], sm.c[5], sm.c[6]),
            artifact['original_VJP_derivatives_to_generic_basis']):
        for k, value in enumerate(coordinates):
            identifier = f'VJP:F{single}:L{k}'
            result[identifier] = result.get(identifier, 0) + 2*coefficient*sp.Rational(value)
    result[f'VJP:Y:F{single}'] = 2*sm.c[2]
    result['VJP:fminus'] = 2*sm.c[3]
    for coefficient, coordinates in zip(
            (sm.d[0], sm.d[2], sm.d[3]),
            artifact['original_VVP_derivatives_to_symmetric_basis']):
        for k, value in enumerate(coordinates):
            identifier = f'VVP:S:F{vsingle}:L{k}'
            result[identifier] = result.get(identifier, 0) + coefficient*sp.Rational(value)
    result[f'VVP:Y:F{vsingle}'] = 2*sm.d[1]
    return result


def main():
    eta = sp.Symbol('eta_gamma', real=True)
    ft, derivatives, E, R = photon_tensors(sm.q1, sm.q2, sm.e1, sm.e2, eta)
    ft_ratio = sp.cancel(ft/E)
    derivative_ratio = sp.cancel(derivatives/E)
    assert ft_ratio == -4*eta and derivative_ratio == -2 and sp.expand(E+R) == 0

    # Match the actual Ruiz action coefficient to the positive WZW photon
    # operator assumed by the inherited ALP coefficient record. With
    # Phi_A=sqrt(2)*Phi_R and f_pi=F these absolute coefficients coincide.
    trace_q = sp.Symbol('Tr_Q2_Phi_R', real=True)
    wzw_action = -sp.sqrt(2)*sm.NC*trace_q/(8*sp.pi**2*sm.F)
    inherited_wzw = sp.sqrt(2)*sm.NC*trace_q/(16*sp.pi**2*sm.F)
    orientation_equation = sp.Eq(wzw_action*derivative_ratio,
                                  inherited_wzw*ft_ratio)
    orientation = sp.solve(orientation_equation, eta)[0]
    assert orientation == -1

    cgamma, Cgamma = sp.symbols('c_gamma C_agammagamma', real=True)
    fa = pc.fa
    # Strip e**2. The explicit photon term has coefficient e**2/(16*pi²*fa).
    G_contact = sp.factor(Cgamma*ft_ratio.subs(eta, orientation)/(16*sp.pi**2*fa))
    G_contact_common_dual = sp.factor(Cgamma*ft_ratio.subs(eta, 1)/(16*sp.pi**2*fa))
    assert G_contact == Cgamma/(4*sp.pi**2*fa)

    TP = sp.diag(*sp.symbols('T_Pu T_Pd T_Ps', real=True))
    FQ = sp.trace((pc.Q*pc.Q+pc.Q*pc.Q)*TP)
    FH = sp.trace((pc.Q*pc.Q+pc.Q*pc.Q)*(TP*pc.CHI+pc.CHI*TP))
    _, _, sm_wzw, sm_direct, sm_exchange = sm.derive()
    replacements = {sm.FQ: FQ, sm.FH: FH, sm.m**2: pc.w,
                    sm.s1: pc.s, sm.s2: pc.t}
    sm_wzw = sm_wzw.subs(replacements, simultaneous=True)
    expected_resonance = (sm_direct+sm_exchange).subs(replacements, simultaneous=True)
    W1 = sm.FV*pc.Q/(sm.MV**2-pc.s)
    W2 = sm.FV*pc.Q/(sm.MV**2-pc.t)
    USM = -sp.sqrt(2)*TP/pc.F
    YSM = sp.sqrt(2)*(TP*pc.CHI+pc.CHI*TP)/pc.F
    gradients = resonance_gradients(USM, YSM, W1, W2)
    original = original_coefficients()
    original_resonance = sum(value*gradients[key] for key, value in original.items())
    original_resonance = sp.factor(original_resonance.subs(pc.h, (pc.w-pc.s-pc.t)/2))
    assert sp.cancel(original_resonance-expected_resonance) == 0

    # Positive SM F in the legacy R convention is a rephased meson amplitude.
    # The common-source E convention needs no meson rephasing. Its scalar
    # anomaly and resonance coefficients equal the same positive-SM scalar.
    common_wzw = sp.factor(wzw_action.subs(trace_q, sp.trace(pc.Q**2*TP))*derivative_ratio)
    assert sp.expand(common_wzw-sm_wzw) == 0
    H = sp.diag(*sp.symbols('H_u H_d H_s', real=True))
    replacement_meson = sp.sqrt(2)*pc.F*H/fa
    mixing_wzw = sp.factor(common_wzw.subs(
        {TP[i, i]: replacement_meson[i, i] for i in range(3)}))
    C_mix = sp.factor(2*sm.NC*sp.trace(pc.Q**2*H))
    assert sp.expand(mixing_wzw-C_mix/(4*sp.pi**2*fa)) == 0
    assembled_seed = cgamma+C_mix
    zero_resonance = sp.factor(cgamma/(4*sp.pi**2*fa)+mixing_wzw)
    assert sp.expand(zero_resonance-assembled_seed/(4*sp.pi**2*fa)) == 0

    # Full photon direct-plus-meson vs no-mixing assembly, for every raw parent.
    # H=cA+S, independently of the already-derived rational source profile.
    A, P, S = [sp.diag(*sp.symbols(prefix+'_u '+prefix+'_d '+prefix+'_s', real=True))
               for prefix in ('cA', 'cP', 'S')]
    Udir, Umix = 2*A/pc.fa, -2*(A+S)/pc.fa
    Ydir = -2*(P*pc.CHI+pc.CHI*P)/pc.fa
    Ymix = 2*((A+S)*pc.CHI+pc.CHI*(A+S))/pc.fa
    w1, w2 = [sp.diag(*sp.symbols(prefix+'_u '+prefix+'_d '+prefix+'_s', real=True))
              for prefix in ('W1', 'W2')]
    direct_gradients = resonance_gradients(Udir, Ydir, w1, w2)
    mixing_gradients = resonance_gradients(Umix, Ymix, w1, w2)
    total_gradients = resonance_gradients(Udir+Umix, Ydir+Ymix, w1, w2)
    bose_substitutions = {pc.s: pc.t, pc.t: pc.s,
                         **dict(zip(w1.diagonal(), w2.diagonal())),
                         **dict(zip(w2.diagonal(), w1.diagonal()))}
    lsz_checks = 0
    for identifier, total in total_gradients.items():
        assert sp.expand(direct_gradients[identifier]+mixing_gradients[identifier]-total) == 0
        assert sp.expand(total-total.xreplace(bose_substitutions)) == 0
        lsz_checks += 1
    kappa = sp.diag(*sp.symbols('kappa_u kappa_d kappa_s', real=True))
    seed = cgamma+2*sm.NC*sp.trace(pc.Q**2*(A+S))
    rotated = cgamma-2*sm.NC*sp.trace(pc.Q**2*kappa)+2*sm.NC*sp.trace(pc.Q**2*(A+S+kappa))
    assert sp.expand(rotated-seed) == 0

    # A nonzero interference check detects an isolated anomaly or resonance
    # sign reversal, which widths of the separate pieces would fail to detect.
    benchmark = {pc.s: 0, pc.t: 0, pc.w: sp.Rational(1, 4), pc.F: sp.Rational(1, 10),
                 sm.FV: sp.Rational(1, 5), pc.MV: 1, sm.NC: 3,
                 TP[0, 0]: 1/sp.sqrt(2), TP[1, 1]: -1/sp.sqrt(2), TP[2, 2]: 0,
                 **{x: 0 for x in sm.c+sm.d}}
    benchmark[sm.c[0]], benchmark[sm.d[0]] = 1, 1
    anomaly_benchmark = sp.simplify(common_wzw.subs(benchmark))
    resonance_benchmark = sp.simplify(original_resonance.subs(benchmark))
    interference = sp.simplify(2*anomaly_benchmark*resonance_benchmark)
    assert interference != 0
    assert sp.simplify((anomaly_benchmark+resonance_benchmark)**2
                      -(-anomaly_benchmark-resonance_benchmark)**2) == 0
    wrong_relative_difference = sp.simplify((anomaly_benchmark+resonance_benchmark)**2
                                            -(-anomaly_benchmark+resonance_benchmark)**2)
    assert wrong_relative_difference == 2*interference

    # Independent component contact Ward and Bose checks.
    for e1, e2 in ((sm.q1, sm.e2), (sm.e1, sm.q2)):
        ward_ft, ward_d, _, _ = photon_tensors(sm.q1, sm.q2, e1, e2, eta)
        assert ward_ft == 0 and ward_d == 0
    ft_swap, d_swap, E_swap, _ = photon_tensors(sm.q2, sm.q1, sm.e2, sm.e1, eta)
    assert all(sp.expand(x-y) == 0 for x, y in zip((ft, derivatives, E), (ft_swap, d_swap, E_swap)))

    paths = ['scripts/derive_odd_photon_phase_matching.py',
             'scripts/derive_radiative_tff_amplitude.py',
             'scripts/derive_odd_power_counted_channels.py',
             'data/vector_sector/odd_leading_action_count.json',
             'memory/extracts/Ruiz-Femenia:2003jdx/pseudoscalar-two-photon-form-factor.md',
             'memory/alp/photon/alp-photon-invariant-from-mixing.md']
    artifact = ROOT/'data/vector_sector/odd_photon_phase_matching.json'
    result = dict(scope='conditional photon convention matching; neutral cubic bounded action',
        source_epsilon_lower0123=1, source_epsilon_upper0123=-1,
        inherited_photon_dual_upper0123=1, inherited_photon_dual_relative_orientation=int(orientation),
        orientation_evidence='convention completion required by the inherited positive WZW assumption, not supplied by Alda',
        Ftilde_two_photon_over_E=str(ft_ratio), epsilon_derivative_two_photon_over_E=str(derivative_ratio),
        tensor_bridge='E=epsilon(q1,epsilon1,q2,epsilon2)=-R; R=epsilon(epsilon1,epsilon2,q1,q2)',
        G_contact_over_e2=str(G_contact), G_contact_common_raised_dual=str(G_contact_common_dual),
        coefficient_dictionary='c_gamma_common_raised=-c_gamma_inherited; C_common_raised=-C_inherited',
        amplitude_E='i M = i e^2 E [C_inherited/(4 pi^2 fa) + sum g_parent G_parent(U,Y;W1,W2)]',
        amplitude_R='i M = -i e^2 R [C_inherited/(4 pi^2 fa) + sum g_parent G_parent(U,Y;W1,W2)]',
        zero_resonance_meson_LSZ_passed=True, original_11_nonzero_resonance_phase_passed=True,
        full_photon_parent_LSZ_checks=lsz_checks, photon_parent_Bose_checks=lsz_checks,
        contact_Ward_checks=2, contact_Bose_passed=True, coefficient_rotation_passed=True,
        interference_benchmark=dict(anomaly=str(anomaly_benchmark),resonance=str(resonance_benchmark),
                                    interference=str(interference),wrong_relative_sign_squared_difference=str(wrong_relative_difference)),
        limitations=['UV photon coefficient is an external signed input in its declared dual convention',
                     'no extra anomalous resonance coefficients or finite gauged-WZW counterterms matched',
                     'tree, linear 1/fa, neutral diagonal profiles, one-nonet bounded action',
                     'perturbative ALP-meson overlaps and LO resolvent pole exclusions still apply'],
        sha256={p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths})
    artifact.write_text(json.dumps(result, indent=2)+'\n')
    transcript = ROOT/'data/vector_sector/odd_photon_phase_matching.md'
    with transcript.open('w') as handle, contextlib.redirect_stdout(handle):
        print('# Shared photon dual and odd-action phase\n')
        for label, expression in [
            ('Source metric', METRIC), ('Common tensor E', E), ('Legacy tensor R', R),
            ('E plus R', sp.expand(E+R)), ('Two Wick assignments of F Ftilde, divided by E', ft_ratio),
            ('Two Wick assignments of epsilon dA dA, divided by E', derivative_ratio),
            ('Ruiz WZW derivative coefficient, stripping e squared', wzw_action),
            ('Inherited positive WZW photon coefficient, stripping e squared', inherited_wzw),
            ('Orientation matching equation', orientation_equation),
            ('Required relative photon-dual orientation', orientation),
            ('ALP contact scalar in common E tensor', G_contact),
            ('ALP contact for common metric-raised photon dual', G_contact_common_dual),
            ('Source-meson WZW scalar in E tensor', common_wzw),
            ('Original 11-coupling resonance scalar in E tensor', original_resonance),
            ('Source meson LSZ replacement', replacement_meson),
            ('WZW meson LSZ scalar', mixing_wzw), ('Conditional seed', assembled_seed),
            ('Zero-resonance assembled scalar', zero_resonance),
            ('Nonzero benchmark anomaly', anomaly_benchmark), ('Nonzero benchmark resonance', resonance_benchmark),
            ('Anomaly-resonance interference', interference),
            ('Squared-amplitude change from an isolated relative-sign error', wrong_relative_difference)]:
            pc.show(label, expression)
        print(f'\nAll {lsz_checks} photon-parent LSZ and {lsz_checks} Bose identities passed.\n')
        print('The inherited unspecified dual orientation is completed by an explicit convention dictionary. '
              'No UV anomaly matching coefficient is inferred.\n')
    print(artifact.relative_to(ROOT))
    print(transcript.relative_to(ROOT))


if __name__ == '__main__':
    main()
