"""SM flavor algebra of the vector radiative channels V -> P gamma.

Question.  The published RChT analysis of the odd-parity sector treats exactly one
on-shell vector radiative decay, omega -> pi0 gamma (Ruiz-Femenia:2003jdx, Eq. 27).
To calibrate any other channel one needs the *operator-level* flavor invariant of that
channel.  This script computes it.

Scope correction, 2026-10-04.  V -> P gamma is a purely Standard Model process: the
final-state pseudoscalar is a real SM meson.  No ALP quantity enters.  In particular the
ALP spurion M_a^{-1}(C) and c_q^AP belong to the a V gamma vertex and must NOT be used
here.  The flavor invariants are built only from

    * the meson state flavor matrices in the (u,d,s) basis,
    * the electric charge matrix Q = diag(2/3,-1/3,-1/3), from f_+ = -2 e Q F,
    * the quark-mass matrix Mcal = diag(m_u,m_d,m_s), from the chi_- meson weight.

Field normalization (Ruiz-Femenia:2003jdx, source lines 198-232 and 291-304):
    Phi     = lambda^a phi^a / sqrt(2),   so Phi_12 = pi^+, Phi_13 = K^+, Phi_23 = K^0,
    V_{munu} = lambda^a V^a_{munu} / sqrt(2), singlet added as diag(1,1,1)/sqrt3.
The "unit field" of octet adjoint a is therefore lambda^a/sqrt(2) and of the singlet
diag(1,1,1)/sqrt(3).  States are NOT rescaled to unit Tr(T^2): doing so silently changes
the eta-eta' and omega-phi mixing angles, which this script quantifies in its last section.

The meson structures entering the flavor traces are
    Pi  = T_P                    (from u_mu = -d_mu Phi, one meson, one momentum)
    Xi  = {T_P, Mcal}            (from chi_-, whose one-meson content is linear in Mcal
                                   and Phi; expanding u = exp(i Phi/(sqrt2 F)) gives
                                   chi_- = 4 i B0 p - (2 i sqrt2 B0/F) {Mcal,Phi}).
So every operator is a single trace of a triple matrix product and the flavor factor is
    charge type   F_Q = tr({T_V, Q}    Pi)
    mass type     F_M = tr({T_V, Q}    Xi)
for the VJP operators, and the same with Q replaced by the internal octet vector matrix
for the VVP exchange operators.

Units.  Purely dimensionless flavor algebra; masses enter only through Mcal.  No ALP
quantity appears anywhere.  Deterministic, no random numbers.

Invocation:
    .venv/bin/python scripts/derive_radiative_sm_flavor_algebra.py
"""

from __future__ import annotations

import sympy as sp

SQ2, SQ3, SQ6, SQ12 = sp.sqrt(2), sp.sqrt(3), sp.sqrt(6), sp.sqrt(12)
R2, R3 = sp.sqrt(sp.Rational(2)), sp.sqrt(sp.Rational(3))

# --- unit fields -----------------------------------------------------------
E1 = sp.diag(1, 1, 1) / SQ3                       # singlet  omega_1 / sqrt(3)
E8 = sp.diag(1, 1, -2) / SQ12                     # octet    lambda^8 / sqrt(2)

PI0 = sp.diag(1, -1, 0) / SQ2                    # lambda^3/sqrt(2)
ETA8 = E8
RHO0 = PI0
PIPLUS = sp.Matrix([[0, 1, 0], [0, 0, 0], [0, 0, 0]])   # Phi_12
RHOPLUS = PIPLUS
K0 = sp.Matrix([[0, 0, 0], [0, 0, 1], [0, 0, 0]])       # Phi_23
KPLUS = sp.Matrix([[0, 0, 1], [0, 0, 0], [0, 0, 0]])     # Phi_13

# Ideal mixing, source lines 784-797, imposed on the unit fields.
OMEGA = sp.simplify(sp.sqrt(sp.Rational(2, 3)) * E1 + sp.sqrt(sp.Rational(1, 3)) * E8)
PHI = sp.simplify(-sp.sqrt(sp.Rational(1, 3)) * E1 + sp.sqrt(sp.Rational(2, 3)) * E8)

# Ideal mixing is imposed on the unit fields; verify the closed forms.
assert sp.simplify(OMEGA - sp.diag(2 * SQ2 + 1, 2 * SQ2 + 1, 2 * SQ2 - 2) / 6) == sp.zeros(3), \
    "T_omega closed form wrong"
assert sp.simplify(PHI - sp.diag(SQ2 - 2, SQ2 - 2, -2 - 2 * SQ2) / 6) == sp.zeros(3), \
    "T_phi closed form wrong"

# eta-eta' mixing: one orthogonal rotation of the unit fields, angle declared.
TH = sp.Symbol("sin_theta", real=True)
CT = sp.Symbol("cos_theta", real=True)
ETA = sp.simplify(CT * ETA8 + TH * E1)
ETAPRIME = sp.simplify(-TH * ETA8 + CT * E1)
# Aloni:2018vki value: |sin theta| = 1/3, cos theta = 2 sqrt(2)/3
ALONI_TH, ALONI_CT = sp.Rational(-1, 3), 2 * R2 / 3
assert sp.simplify(ALONI_TH**2 + ALONI_CT**2) == 1, "eta-eta' rotation must be orthogonal"

Q = sp.diag(sp.Rational(2, 3), sp.Rational(-1, 3), sp.Rational(-1, 3))
MU, MD, MS = sp.symbols("m_u m_d m_s", positive=True)
MCAL = sp.diag(MU, MD, MS)


#: the internal octet vector exchanged in the VVP diagrams is always the rho
INTERNAL_OCTET = sp.diag(1, -1, 0) / SQ2          # lambda^3 / sqrt(2)


def exchange_flavor(vector: sp.Matrix, pseudoscalar: sp.Matrix) -> tuple[sp.Expr, sp.Expr]:
    """Flavor invariants for the VVP exchange, where the second vector is the
    internal OCTET line rather than the external state.

    The photon reaches the internal line through L_2(V), so the exchanged resonance
    is the octet rho in every channel, and the flavor trace is
        F_Q^exch = tr({T_V, T_rho} T_P)          -> d1, d3
        F_M^exch = tr({T_V, T_rho} {T_P,Mcal})    -> d2
    which differs from the direct invariants because T_rho replaces Q.
    """
    a = anticommutator(vector, INTERNAL_OCTET)
    f_q = sp.simplify(sp.trace(a * pseudoscalar))
    f_m = sp.simplify(sp.trace(a * anticommutator(pseudoscalar, MCAL)))
    return f_q, f_m


def crossing(vector: sp.Matrix) -> sp.Matrix:
    """Flavor matrix of the crossed (vector) partner of a pseudoscalar: T* = T^T.

    This applies only where the two fields really are crossing partners of one another,
    i.e. to a pseudoscalar and the vector of the same quark content (kaons, charged pions).
    It does NOT apply to nonet channels, where T_V and T_P are different states that are not
    each other's crossing partners; there both matrices are carried independently.
    """
    return vector.T


def swap(vector: sp.Matrix, pseudoscalar: sp.Matrix) -> tuple[sp.Matrix, sp.Matrix]:
    """Reorder a (vector, pseudoscalar) pair into a (pseudoscalar, vector) pair.

    P -> V gamma is the same operator vertex as V -> P gamma: the VJP operators
    contain exactly one V, one f_+ and one u, and which of the two mesonic fields is
    the decaying external state is a choice of external on-shell substitution, not a
    different interaction.  P -> V gamma is therefore a first-class channel, and PDG
    tabulates B(eta' -> rho0 gamma) = (29.5 +/- 0.4) % and
    B(eta' -> omega gamma) = (2.52 +/- 0.07) %.  The two orientations share the flavor
    invariant because the trace is symmetric under the exchange for diagonal states,
    but they differ in kinematics and in phase space, so neither is redundant with
    the other: eta' -> rho0 gamma is open while rho -> eta' gamma is closed.
    """
    return pseudoscalar, vector


CHANNELS: tuple[tuple[str, sp.Matrix, sp.Matrix], ...] = (
    ("omega -> pi0 gamma", OMEGA, PI0),
    ("rho   -> pi0 gamma", RHO0, PI0),
    ("rho+  -> pi+ gamma", crossing(RHOPLUS), PIPLUS),
    ("omega -> eta  gamma", OMEGA, ETA),
    ("rho   -> eta  gamma", RHO0, ETA),
    ("phi   -> eta  gamma", PHI, ETA),
    ("omega -> eta' gamma", OMEGA, ETAPRIME),
    ("rho   -> eta' gamma", RHO0, ETAPRIME),
    ("phi   -> eta' gamma", PHI, ETAPRIME),
    ("K*0   -> K0   gamma", crossing(K0), K0),
    ("K*+   -> K+   gamma", crossing(KPLUS), KPLUS),
    # pseudoscalar-initiated channels: same vertex, roles of the external states swapped
    ("eta'  -> rho0 gamma", *swap(RHO0, ETAPRIME)),
    ("eta'  -> omega gamma", *swap(OMEGA, ETAPRIME)),
    ("eta'  -> phi   gamma", *swap(PHI, ETAPRIME)),
)


def anticommutator(a: sp.Matrix, b: sp.Matrix) -> sp.Matrix:
    return sp.simplify(a * b + b * a)


def flavor(vector: sp.Matrix, pseudoscalar: sp.Matrix) -> tuple[sp.Expr, sp.Expr]:
    """Return (F_Q, F_M) for one channel.

    F_Q = tr({T_V, Q} T_P)      charge-type flavor invariant (VJP O1,O2,O5,O6; VVP P1,P3)
    F_M = tr({T_V, Q} {T_P,Mcal})  mass-type flavor invariant (VJP O3; VVP P2)

    The ordering inside the trace follows the operators: the fields appear as
    {V, f_+} and the meson, so the charge matrix sits between the vector and the meson.
    """
    aq = anticommutator(vector, Q)
    f_q = sp.simplify(sp.trace(aq * pseudoscalar))
    f_m = sp.simplify(sp.trace(aq * anticommutator(pseudoscalar, MCAL)))
    return f_q, f_m


def main() -> None:
    print("Meson state flavor matrices, Ruiz-Femenia field normalization")
    print("  (unit Tr(T^2) is NOT imposed; it would change the mixing angles)")
    for nm, mat in (("pi0", PI0), ("eta8", ETA8), ("singlet", E1), ("eta", ETA),
                    ("eta'", ETAPRIME), ("rho0", RHO0), ("omega", OMEGA), ("phi", PHI),
                    ("pi+", PIPLUS), ("K0", K0), ("K+", KPLUS)):
        diag = [sp.simplify(x) for x in mat.diagonal()]
        tr2 = sp.simplify(sp.trace(mat * mat))
        print(f"  {nm:8s} diag = {diag}   Tr(T^2) = {tr2}")
    assert sp.simplify(sp.trace(PIPLUS * PIPLUS)) == 0, "charged pion has zero unit trace"
    assert sp.simplify(sp.trace(K0 * K0)) == 0, "kaon has zero unit trace"

    print()
    print("=" * 78)
    print("1. Flavor invariants of each channel (no ALP quantity anywhere)")
    print("=" * 78)
    print("   F_Q = tr({T_V,Q} T_P)          -> VJP c1,c2,c5,c6,c7 and VVP d1,d3,d4")
    print("   F_M = tr({T_V,Q} {T_P,Mcal})    -> VJP c3 and VVP d2")
    table: dict[str, tuple[sp.Expr, sp.Expr]] = {}
    for name, vec, pseu in CHANNELS:
        fq, fm = flavor(vec, pseu)
        table[name] = (fq, fm)
        print(f"\n  {name}")
        print(f"    F_Q = {sp.collect(sp.expand(fq), (MU, MD, MS)) if MU in fq.free_symbols else sp.expand(fq)}")
        print(f"    F_M = {sp.collect(sp.expand(fm), (MU, MD, MS)) if MU in fm.free_symbols else sp.expand(fm)}")

    print()
    print("=" * 78)
    print("1b. Exchange flavor invariants, internal octet = lambda^3/sqrt(2)")
    print("=" * 78)
    print("   F_Q^exch = tr({T_V,T_rho} T_P)         -> d1, d3, d4(vanishes on shell)")
    print("   F_M^exch = tr({T_V,T_rho} {T_P,Mcal})   -> d2")
    ex: dict[str, tuple[sp.Expr, sp.Expr]] = {}
    for name, vec, pseu in CHANNELS:
        fq, fm = exchange_flavor(vec, pseu)
        ex[name] = (fq, fm)
        print(f"  {name:22s} F_Q^exch = {sp.simplify(fq)}    F_M^exch = {sp.simplify(fm)}")
    print()
    print("  Two structural results, both assertions of the script:")
    print("   (a) for the ISOVECTOR final states the charge-type invariant vanishes,")
    print("       because {T_rho,T_rho} = 2 T_rho^2 is proportional to diag(1,1,0) and the")
    print("       isovector is orthogonal to it.  So d1 and d3 do not enter, while d2 can,")
    print("       carried only by the isospin-breaking quark-mass combination:")
    for name in ("rho   -> pi0 gamma", "rho+  -> pi+ gamma"):
        fq, fm = ex[name]
        assert sp.simplify(fq) == 0, f"{name}: F_Q^exch must vanish"
        print(f"       {name:22s} F_Q^exch = 0   F_M^exch = {sp.simplify(fm)}   ok")
    print("       rho+ -> pi+ gamma has F_M^exch = 0 identically: no exchange at all.")
    print("       rho -> pi0 gamma has F_M^exch proportional to (m_u - m_d), i.e. pure")
    print("       isospin breaking, vanishing in the exact SU(3) limit.")
    print("   (b) the kaon channels DO exchange: {T_K*,T_rho} is off diagonal but")
    print("       {E_32,E_23} contracts onto E_22, so the d component survives:")
    for name in ("K*0   -> K0   gamma", "K*+   -> K+   gamma"):
        fq, fm = ex[name]
        assert sp.simplify(fq) != 0, f"{name} must have an exchange contribution"
        print(f"       {name:22s} F_Q^exch = {sp.simplify(fq)}   ok")
    print("  So the rho-exchange is present for omega, phi, eta' and the kaons, and is")
    print("  absent only for the isovector channels rho -> pi0 gamma and")
    print("  rho+ -> pi+ gamma.  rho -> eta gamma keeps an exchange contribution because")
    print("  eta has a large isoscalar (singlet) part, so {T_rho,T_rho} does not cancel.")
    print("  These are predictions of the flavor algebra, not inputs.")

    print()
    print("=" * 78)
    print("2. Kinematic openness  (M_V^ext > m_P), PDG-scale masses in GeV")
    print("=" * 78)
    M_P = {"pi0": 0.1396, "pi+": 0.1396, "eta": 0.5479, "eta'": 0.9578,
           "K0": 0.4976, "K+": 0.4937, "rho0": 0.7753, "omega": 0.7827,
           "phi": 1.0195}
    M_V = {"omega": 0.7827, "rho0": 0.7753, "rho+": 0.7753, "phi": 1.0195,
           "K*0": 0.8956, "K*+": 0.8922, "eta'": 0.9578}
    #: label of the decaying state -> label of the final-state meson, per channel
    DECAYING, FINAL = {}, {}
    for name, _, _ in CHANNELS:
        parts = name.replace("->", " ").split()      # [decaying, final meson, "gamma"]
        assert parts[2] == "gamma", name
        DECAYING[name] = parts[0] + ("0" if parts[0] == "rho" and parts[1] != "+" else "")
        DECAYING[name] = DECAYING[name].replace("rho+0", "rho+")
        FINAL[name] = parts[1]
    OPEN: dict[str, bool] = {}
    for name, _, _ in CHANNELS:
        v, ps = DECAYING[name], FINAL[name]
        opened = M_V[v] > M_P[ps]
        OPEN[name] = opened
        print(f"  {name:22s} M_initial = {M_V[v]:.4f}   m_final = {M_P[ps]:.4f}   "
              f"{'open' if opened else 'CLOSED'}")
    for closed in ("omega -> eta' gamma", "rho   -> eta' gamma"):
        assert not OPEN[closed], f"{closed} should be kinematically closed"
    # the swapped orientations are open exactly where the direct ones are closed
    for direct, swapped in (("rho   -> eta' gamma", "eta'  -> rho0 gamma"),
                            ("omega -> eta' gamma", "eta'  -> omega gamma")):
        assert not OPEN[direct] and OPEN[swapped], \
            f"{swapped} must be open while {direct} is closed"

    print()
    print("=" * 78)
    print("3. Exact SU(3) limit: F_M becomes proportional to F_Q")
    print("=" * 78)
    iso = {MU: sp.Symbol("m"), MD: sp.Symbol("m"), MS: sp.Symbol("m")}
    for name, vec, pseu in CHANNELS:
        fq, fm = table[name]
        fqi, fmi = sp.expand(fq).subs(iso), sp.expand(fm).subs(iso)
        ratio = sp.simplify(fmi / fqi) if sp.simplify(fqi) != 0 else sp.oo
        print(f"  {name:22s} F_Q = {sp.factor(fqi)}   F_M/F_Q = {sp.simplify(ratio)}")

    print()
    print("=" * 78)
    print("4. Redundancy of the flavor algebra at general quark masses")
    print("=" * 78)
    seen: dict[tuple, list[str]] = {}
    for name, (fq, fm) in table.items():
        key = (sp.simplify(fq / (MU + 2 * MD)) if sp.simplify(fq) != 0 else sp.Integer(0),
               sp.simplify(fm))
        seen.setdefault(key, []).append(name)
    for i, (key, members) in enumerate(sorted(seen.items(), key=lambda kv: kv[1][0]), 1):
        print(f"  group {i}: {', '.join(members)}")
    print(f"\n  {len(CHANNELS)} channels fall into {len(seen)} distinct flavor classes.")

    print()
    print("=" * 78)
    print("5. Correction of the nonet matrices used in the earlier record")
    print("=" * 78)
    old = {"omega": sp.diag(1, 1, 0) / SQ2, "phi": sp.diag(0, 0, -1),
           "eta": sp.diag(1, 1, -1) / SQ3, "eta'": sp.diag(1, 1, 2) / SQ6}
    new = {"omega": OMEGA, "phi": PHI, "eta": ETA, "eta'": ETAPRIME}

    a, b = sp.symbols("a b")

    def mix(T):
        sol = sp.solve([sp.Eq(T[i, i], a * E1[i, i] + b * E8[i, i]) for i in range(3)],
                       (a, b), dict=True)[0]
        return sp.simplify(sol[a]), sp.simplify(sol[b])

    for nm in ("omega", "phi", "eta", "eta'"):
        sub = {TH: ALONI_TH, CT: ALONI_CT} if nm in ("eta", "eta'") else {}
        a_new, b_new = mix(sp.simplify(new[nm].subs(sub)))
        a_old, b_old = mix(old[nm])
        print(f"  {nm}: singlet/octet weight ratio |a/b|")
        print(f"     correct (unit fields)   : {sp.nsimplify(sp.sqrt(a_new**2 / b_new**2))}")
        print(f"     earlier unit-trace matrix: {sp.nsimplify(sp.sqrt(a_old**2 / b_old**2))}")
    print("\n  Both eta and eta' are octet-dominant / singlet-dominant respectively,")
    print("  as required by the observed mass ordering, but the mixing weights differ by a")
    print("  factor sqrt(2): rescaling each state matrix independently to unit Tr(T^2) is")
    print("  not a harmless normalization, it changes the mixing angle and hence the flavor")
    print("  invariants above.")

    print()
    print("No ALP quantity was used.  All assertions above passed.")


if __name__ == "__main__":
    main()
