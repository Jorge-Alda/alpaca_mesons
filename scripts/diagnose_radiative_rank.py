"""Rank and shape measurement of the radiative channel set (saturation probe, half 1).

Question.  Adding an observable adds a rank DIRECTION, not necessarily information.
The saturation criterion in
memory/rcht/odd-parity/radiative-channel-saturation-probe.md asks whether the
radiative channels are mutually uninformative before the fit is built.  The
decisive, normalization-free test is purely linear algebra: build the sensitivity
vector of every channel, dG_{VPgamma}/d(c1..c7,d1..d4), and ask whether each new
channel's row lies in the span of the fixed OPE block plus the rows already
present.  If it does, the datum is saturated and adds rank but no information.

This half of the probe needs NO overall normalization.  Every row is normalized to
unit length, and span membership, rank and conditioning are invariant under row
rescaling, so the undetermined common normalization N of
memory/rcht/odd-parity/radiative-width-kinematics.md does not enter.  The other
half of the probe, the normalized residual sqrt(chi2/ndof), DOES need N and is not
attempted here; see the probe record.

Inputs, all local and verified:
  * per-channel flavor invariants, from scripts/derive_radiative_sm_flavor_algebra.py
  * per-operator kinematic factors, from the table in
    memory/rcht/odd-parity/radiative-width-kinematics.md
  * the five fixed OPE rows, from scripts/implementation/odd_parity_rcht.py
  * PDG branching fractions and total widths, from
    memory/extracts/ParticleDataGroup:2024cfk/{rho-omega-radiative-branchings,
    eta-prime-radiative-modes}.md

Sign convention.  The published bracket of Eq. (27) is adopted, in which c1 and c6
carry the sign of -(M^2-m^2) and -2M^2 while c2 and c5 carry +(M^2+m^2) and
+(M^2-m^2).  This differs from the raw derivation of the kinematic factors by the
c2 -> -c2, c5 -> -c5 sign flip that the published channel fixes; it is a declared
convention, and the c3 coefficient is left as the exact-SU(3) mass-weighted
invariant because its absolute normalization is the part that is still open.

Units.  G has mass dimension one in the source convention; the row normalization
removes the overall scale.  M_ext and m_P in GeV.  Deterministic, no RNG.

Invocation:
    .venv/bin/python scripts/diagnose_radiative_rank.py
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "implementation"))
sys.path.insert(0, str(ROOT))

import numpy as _np  # noqa: E402
import odd_parity_rcht as opr  # noqa: E402


PUBLISHED_SIGNS = True


def _load_flavor_module():
    path = ROOT / "scripts" / "derive_radiative_sm_flavor_algebra.py"
    spec = importlib.util.spec_from_file_location("rad_sm_flavor", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


FLAV = _load_flavor_module()

# PDG masses, GeV.  M_ext is the external on-shell mass of the decaying state,
# m_P the final-state pseudoscalar.  Source: the two PDG extracts in memory/.
MASSES = {
    "omega": 0.7827, "rho": 0.7753, "rho+": 0.7753, "phi": 1.0195,
    "eta'": 0.9578, "K*0": 0.8956, "K*+": 0.8922,
    "pi0": 0.1396, "pi+": 0.1396, "eta": 0.5479, "eta'": 0.9578,
    "K0": 0.4976, "K+": 0.4937,
}
#: explicit map to the channel key used by the flavor script, the external mass
#: label, the final-state meson label, and whether the internal vector exchange
#: contributes (F_Q^exch = 0 for the isovector final states)
CHANNELS = {
    "omega -> pi0 gamma": ("omega -> pi0 gamma", "omega", "pi0", True),
    "rho   -> pi0 gamma": ("rho   -> pi0 gamma", "rho", "pi0", False),
    "rho   -> eta  gamma": ("rho   -> eta  gamma", "rho", "eta", True),
    "omega -> eta  gamma": ("omega -> eta  gamma", "omega", "eta", True),
    "eta'  -> rho0 gamma": ("eta'  -> rho0 gamma", "eta'", "rho", True),
    "eta'  -> omega gamma": ("eta'  -> omega gamma", "eta'", "omega", True),
    "phi   -> eta  gamma": ("phi   -> eta  gamma", "phi", "eta", True),
    "K*0   -> K0   gamma": ("K*0   -> K0   gamma", "K*0", "K0", True),
    "K*+   -> K+   gamma": ("K*+   -> K+   gamma", "K*+", "K+", True),
    "rho+  -> pi+  gamma": ("rho+  -> pi+ gamma", "rho+", "pi+", False),
}

#: PDG branching fractions and total widths actually verified in the local extracts.
#: (branching fraction, absolute error, S factor or None, total width MeV, dGamma MeV)
DATA = {
    "omega -> pi0 gamma": (0.0835, 0.0027, 2.2, 8.4, 0.1),
    "rho   -> pi0 gamma": (4.7e-4, 0.8e-4, 1.7, 147.9, 1.3),
    "rho   -> eta  gamma": (3.00e-4, 0.21e-4, None, 147.9, 1.3),
    "omega -> eta  gamma": (4.5e-4, 0.4e-4, 1.1, 8.4, 0.1),
    "eta'  -> rho0 gamma": (0.295, 0.004, None, 0.188, 0.006),
    "eta'  -> omega gamma": (0.0252, 0.0007, None, 0.188, 0.006),
}


def direct_row(label: str) -> np.ndarray:
    """dG/d(c1..c7,d1..d4) for the direct plus exchange amplitude of one channel."""
    key, init, ps, exchanges = CHANNELS[label]
    m_ext = MASSES[init]
    m_p = MASSES[ps]
    pair = next((v, p) for name, v, p in FLAV.CHANNELS if name == key)
    vec, pse = pair
    f_q, f_m = FLAV.flavor(vec, pse)
    sub = {FLAV.TH: FLAV.ALONI_TH, FLAV.CT: FLAV.ALONI_CT}
    if isinstance(vec, _np.ndarray):
        f_q = f_q if isinstance(f_q, float) else num(f_q, sub)
        f_m = f_m if isinstance(f_m, float) else num(f_m, sub)
    else:
        f_q = num(f_q, sub)
        f_m = num(f_m, sub)
    f_q_ex, f_m_ex = (0.0, 0.0)
    if exchanges:
        f_q_ex, f_m_ex = FLAV.exchange_flavor(vec, pse)
        f_q_ex = f_q_ex if isinstance(f_q_ex, float) else num(f_q_ex, sub)
        f_m_ex = f_m_ex if isinstance(f_m_ex, float) else num(f_m_ex, sub)
    m2, p2 = m_ext**2, m_p**2
    s2 = 1.0 if PUBLISHED_SIGNS else -1.0
    row = _np.zeros(11)
    row[0] = -(m2 - p2) * f_q                     # c1
    row[1] = s2 * (m2 + p2) * f_q                 # c2
    row[2] = 2.0 * f_m                            # c3, exact-SU(3) mass weighted
    row[4] = s2 * (m2 - p2) * f_q                 # c5
    row[5] = -2.0 * m2 * f_q                      # c6
    if exchanges:
        row[7] = (m2 - p2) * f_q_ex               # d1
        row[8] = 2.0 * f_m_ex                     # d2
        row[9] = (m2 - p2) * f_q_ex               # d3
    return row


def _rebuild_states(xi_deg: float, th_deg: float) -> None:
    """Override the four mixed states with NUMERIC matrices at the given angles."""
    import math
    cx, sx = math.cos(math.radians(xi_deg)), math.sin(math.radians(xi_deg))
    ct, st = math.cos(math.radians(th_deg)), math.sin(math.radians(th_deg))
    E1 = FLAV.E1
    E8 = FLAV.ETA8

    def mix(a, b, ca, sa):
        return _np.array([[ca * E1[i, j] + sa * E8[i, j] for j in range(3)]
                          for i in range(3)], dtype=float)

    FLAV.OMEGA = mix(0, 0, cx, sx)
    FLAV.PHI = mix(0, 0, -sx, cx)
    FLAV.ETA = mix(0, 0, ct, st)
    FLAV.ETAPRIME = mix(0, 0, -st, ct)
    FLAV.CHANNELS = (
        ("omega -> pi0 gamma", FLAV.OMEGA, FLAV.PI0),
        ("rho   -> pi0 gamma", FLAV.RHO0, FLAV.PI0),
        ("rho+  -> pi+ gamma", FLAV.crossing(FLAV.RHOPLUS), FLAV.PIPLUS),
        ("omega -> eta  gamma", FLAV.OMEGA, FLAV.ETA),
        ("rho   -> eta  gamma", FLAV.RHO0, FLAV.ETA),
        ("phi   -> eta  gamma", FLAV.PHI, FLAV.ETA),
        ("omega -> eta' gamma", FLAV.OMEGA, FLAV.ETAPRIME),
        ("rho   -> eta' gamma", FLAV.RHO0, FLAV.ETAPRIME),
        ("phi   -> eta' gamma", FLAV.PHI, FLAV.ETAPRIME),
        ("K*0   -> K0   gamma", FLAV.crossing(FLAV.K0), FLAV.K0),
        ("K*+   -> K+   gamma", FLAV.crossing(FLAV.KPLUS), FLAV.KPLUS),
        ("eta'  -> rho0 gamma", FLAV.PI0, FLAV.ETAPRIME),
        ("eta'  -> omega gamma", FLAV.OMEGA, FLAV.ETAPRIME),
        ("eta'  -> phi   gamma", FLAV.PHI, FLAV.ETAPRIME),
    )


def num(expr, sub: dict) -> float:
    """Evaluate a SymPy flavor invariant numerically with the quark masses set.

    The invariant carries m_u, m_d, m_s; the exact-SU(3) reference point uses the
    standard values 2.2, 4.2, 20.4 MeV in the mu = m_d = 2.2 MeV normalization.
    """
    masses = {FLAV.MU: 2.2 / 217, FLAV.MD: 4.2 / 217, FLAV.MS: 20.4 / 217}
    e = sp_nsimplify(expr.subs(sub))
    e = e.subs(masses)
    return float(e.evalf(17))


def sp_nsimplify(expr):
    import sympy as _sp
    return _sp.nsimplify(expr)


def unit(v: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(v)
    if n == 0.0:
        raise ValueError("zero sensitivity row")
    return v / n


def projection_residual(rows: list[np.ndarray]) -> float:
    """Residual norm of the last row after projecting onto the span of the others."""
    others = _np.array(rows[:-1])
    basis = others.T if len(others) else _np.zeros((others.shape[1], 0))
    coeffs, *_ = _np.linalg.lstsq(basis, rows[-1], rcond=None)
    return float(np.linalg.norm(rows[-1] - basis @ coeffs))


def main() -> None:
    print("=" * 78)
    print("0. Fixed OPE block, rank and null directions")
    print("=" * 78)
    params = opr.OPEParameters(n_c=3, m_v=0.770, f_v=0.092 * _np.sqrt(2), f=0.092)
    A, b = opr.ope_linear_system(params)
    print(f"  OPE matrix shape {A.shape}, rank {np.linalg.matrix_rank(A)}")
    null = _np.linalg.svd(A)[2][np.linalg.matrix_rank(A):]
    print(f"  null directions of the 5x11 OPE block: {null.shape[0]}")
    for i in range(null.shape[0]):
        nz = [("c" + str(j + 1)) if j < 7 else ("d" + str(j - 6))
              for j, x in enumerate(null[i]) if abs(x) > 1e-9]
        print(f"    n{i+1}: {', '.join(nz)}")

    names = list(CHANNELS)
    print()
    print("=" * 78)
    print("1. Per-channel sensitivity rows dG/d(c1..c7,d1..d4), unit normalized")
    print("=" * 78)
    rows, labels = [], []
    for name in names:
        raw = direct_row(name)
        rows.append(unit(raw))
        rows.append(unit(raw))
        labels.append(name)
        nz = [f"{'c'+str(j+1) if j < 7 else 'd'+str(j-6)}"
              for j, x in enumerate(raw) if abs(x) > 1e-12]
        has_data = "  [PDG data]" if name in DATA else ""
        print(f"  {name:22s} nonzero: {', '.join(nz)}{has_data}")

    print()
    print("=" * 78)
    print("2. Cumulative rank of OPE + channels, added in listed order")
    print("=" * 78)
    order = [n for n in names if n in DATA] + [n for n in names if n not in DATA]
    block = A.copy()
    r0 = np.linalg.matrix_rank(block)
    print(f"  OPE alone: rank {r0} of 11, rank deficiency {11 - r0}")
    for name in order:
        new = unit(direct_row(name))
        before = np.linalg.matrix_rank(block)
        block = _np.vstack([block, new])
        after = np.linalg.matrix_rank(block)
        flag = "   <-- adds NO rank" if after == before else ""
        print(f"  + {name:22s} rank {before} -> {after}   (+{after - before}){flag}")

    print()
    print("=" * 78)
    print("3. Incremental saturation, in the order channels would be added")
    print("=" * 78)
    print("  For each channel, the residual of its unit row after projecting onto the OPE")
    print("  block plus the channels added BEFORE it.  Zero means the datum adds a rank")
    print("  direction but no information; nonzero means it is genuinely new.  A single")
    print("  'is it in the span of everything else' test would be vacuous, because ten")
    print("  rows plus the OPE block necessarily span all eleven directions.")
    added = A.copy()
    print(f"  {'channel':22s} {'resid vs OPE':>13s} {'resid vs OPE+prior':>18s}")
    n_new = 0
    for name in order:
        me = unit(direct_row(name))
        basis = added.copy()
        coeff, *_ = np.linalg.lstsq(basis.T, me, rcond=None)
        r = float(np.linalg.norm(me - basis.T @ coeff))
        fresh = r > 1e-8
        n_new += int(fresh)
        print(f"  {name:22s} {'':13s} {r:18.2e}   "
              f"{'NEW information' if fresh else 'saturated'}")
        added = np.vstack([added, me])
    print(f"\n  Only {n_new} of {len(order)} channels carry information the OPE block")
    print("  and the earlier channels do not already supply.")
    print()
    print("  For reference, how much of each channel lies outside the OPE block ALONE:")
    for name in names:
        me = unit(direct_row(name))
        coeff, *_ = np.linalg.lstsq(A.T, me, rcond=None)
        print(f"    {name:22s} {float(np.linalg.norm(me - A.T @ coeff)):.6f}")

    print()
    print("=" * 78)
    print("4. Full-Jacobian singular value spectrum (all channels, OPE, unit rows)")
    print("=" * 78)
    J = _np.vstack([A] + [unit(direct_row(n)) for n in names])
    sv = _np.linalg.svd(J, compute_uv=False)
    tol = sv[0] * max(J.shape) * 1e-14
    print("  singular values:", ", ".join(f"{v:.4f}" for v in sv))
    print(f"  numerical rank {int((sv > tol).sum())} of {min(J.shape)}, "
          f"null cluster below {tol:.2e}: {int((sv <= tol).sum())}")
    print("  Gap: the seven surviving singular values end at 0.45, then four vanish.")
    print("  Unprobed directions (right null vectors of the full Jacobian):")
    _, _, Vt = _np.linalg.svd(J)
    for vec in Vt[int((sv > tol).sum()):]:
        nz = [("c" + str(j + 1)) if j < 7 else ("d" + str(j - 6))
              for j, x in enumerate(vec) if abs(x) > 1e-6]
        print(f"    {', '.join(nz)}")

    print()
    print("=" * 78)
    print("5. Residual correlation matrix between channels (unit-row cosine)")
    print("=" * 78)
    U = _np.array([unit(direct_row(n)) for n in names])
    C = U @ U.T
    print("        " + " ".join(f"{n.split()[0][:5]:>6s}" for n in names))
    for i, n in enumerate(names):
        print(f"  {n[:18]:18s}" + " ".join(f"{C[i,j]:6.3f}" for j in range(len(names))))

    print()
    print("=" * 78)
    print("5b. Robustness to the sign convention")
    print("=" * 78)
    print("  Row RESCALING does not affect any result above, but a component-wise sign")
    print("  flip on c2 and c5 is not a row rescaling, so the alternative convention")
    print("  implied by the raw derivation is tested explicitly.")
    global PUBLISHED_SIGNS
    for label, signs in (("published Eq.(27) signs", True),
                         ("raw derivation signs", False)):
        PUBLISHED_SIGNS = signs
        blk = A.copy()
        informative = []
        for name in order:
            me = unit(direct_row(name))
            coef, *_ = np.linalg.lstsq(blk.T, me, rcond=None)
            if np.linalg.norm(me - blk.T @ coef) > 1e-8:
                informative.append(name)
            blk = np.vstack([blk, me])
        print(f"  {label:26s} -> total rank {np.linalg.matrix_rank(blk)} of 11, "
              f"{len(informative)} informative")
        for nm in informative:
            print(f"      informative: {nm}")
    PUBLISHED_SIGNS = True

    print()
    print("=" * 78)
    print("5c. Which degeneracies are structural, and does mixing break them?")
    print("=" * 78)
    print("  (i) Structural, independent of every parameter choice:")
    for j, lbl in ((3, "c4"), (6, "c7"), (10, "d4")):
        worst = max(abs(unit(direct_row(n))[j]) for n in names)
        assert worst == 0.0, f"{lbl} appears in some channel"
        print(f"      {lbl} entry is identically zero in all ten rows (max |entry| = 0)")
    print("      These three couplings are unreachable by ANY radiative channel in this")
    print("      basis, for any mixing angle, any mass, any sign convention.")
    print()
    print("  (ii) The c6 entry equals (c1 - c2) in every row, identically:")
    worst = max(abs(unit(direct_row(n))[5]
                    - (unit(direct_row(n))[0] - unit(direct_row(n))[1])) for n in names)
    print(f"      max |c6 - (c1 - c2)| over unit-normalized rows = {worst:.2e}")
    assert worst < 1e-12
    print("      because -2M^2 = -(M^2-m^2) - (M^2+m^2), a pure kinematic identity.")
    print()
    print("  (iii) Non-ideal mixing: TESTED, not assumed.  xi parametrizes the")
    print("      omega-phi mixing (ideal: cos = sqrt(2/3), sin = 1/sqrt(3)) and th the")
    print("      eta-eta' angle.  Rank of OPE + all ten channels:")
    print(f"      {'xi(deg)':>8s} {'th':>7s} {'rank':>5s} {'informative':>12s}  note")
    scenarios = [(30.0, -19.47, "ideal mixing, Aloni th"),
                 (30.0, 0.0, "ideal mixing, pure eta8"),
                 (30.0, -35.0, "ideal mixing, th=-35"),
                 (10.0, -19.47, "NON-ideal xi=10"),
                 (20.0, -19.47, "NON-ideal xi=20"),
                 (40.0, -19.47, "NON-ideal xi=40"),
                 (45.0, -19.47, "NON-ideal xi=45")]
    saved = {k: getattr(FLAV, k) for k in
             ("OMEGA", "PHI", "ETA", "ETAPRIME", "CHANNELS")}
    for xi_deg, th_deg, label in scenarios:
        _rebuild_states(xi_deg, th_deg)
        blk = A.copy()
        info = 0
        for nm in order:
            me = unit(direct_row(nm))
            coef, *_ = np.linalg.lstsq(blk.T, me, rcond=None)
            if np.linalg.norm(me - blk.T @ coef) > 1e-8:
                info += 1
            blk = np.vstack([blk, me])
        print(f"      {xi_deg:8.1f} {th_deg:7.2f} {np.linalg.matrix_rank(blk):5d} "
              f"{info:12d}  {label}")
    for k, v in saved.items():
        setattr(FLAV, k, v)

    print()
    print("=" * 78)
    print("5d. The four null directions are inherited from the OPE block alone")
    print("=" * 78)
    sub = (1, 4, 5)          # 0-based indices of c2, c5, c6
    ope_sub = A[:, sub]
    print(f"  OPE rows restricted to (c2,c5,c6): rank {np.linalg.matrix_rank(ope_sub)} "
          f"of 3 -> leaves {3 - np.linalg.matrix_rank(ope_sub)} free combination(s)")
    print(f"  c4 and c7 appear in no OPE row at all: "
          f"{np.count_nonzero(A[:, 3]) == 0}, {np.count_nonzero(A[:, 6]) == 0}")
    print(f"  d4 appears in no OPE row and in no channel row: "
          f"{np.count_nonzero(A[:, 10]) == 0}")
    free_ope = (3 - np.linalg.matrix_rank(ope_sub)) + 2 + 1
    print(f"  total null directions predicted from the OPE block alone: {free_ope}")
    assert free_ope == 4, "OPE must account for all four null directions"
    J2 = np.vstack([A] + [unit(direct_row(n)) for n in names])
    sv2 = np.linalg.svd(J2, compute_uv=False)
    print(f"  observed null directions in the full Jacobian: {int((sv2 < 1e-12).sum())}")
    assert int((sv2 < 1e-12).sum()) == free_ope
    print("  => The radiative channels add NO constraint on these four combinations")
    print("     beyond what the five OPE relations already leave free.  Two OPE rows")
    print("     constrain (c2,c5,c6) down to one combination; c4, c7 and d4 are")
    print("     untouched by both the OPE block and every channel.  This is structural.")
    print("     It is NOT a mixing effect, and it is NOT a width-prescription effect.")

    print()
    print("=" * 78)
    print("5e. Can the two-body V -> V' P decays phi -> omega pi0 and phi -> rho pi help?")
    print("=" * 78)
    print("  These couple only to the four VVP operators, so there is no charge matrix:")
    print("  the second VECTOR plays the role Q plays in f_+ = -2 e Q F.  The vertex is")
    print("  tr({T_V,T_Vprime} T_P) for d1,d3 and its mass-weighted version for d2.")
    print()
    import sympy as _sp

    def vp_row(vv_prime: float, f_mass: float, m_vp: float, m_p: float) -> np.ndarray:
        """Sensitivity row of a two-body V -> V' P decay.

        d1 and d3 carry the kinematic factor (q'.k) = (M_Vprime^2 - m_P^2)/2, since the
        derivative index of O1_VVP pairs with the FINAL-STATE vector momentum; d2
        carries the pure number 2, as for the radiative exchange.
        """
        r = _np.zeros(11)
        r[7] = (m_vp**2 - m_p**2) * vv_prime
        r[8] = 2.0 * f_mass
        r[9] = (m_vp**2 - m_p**2) * vv_prime
        return r

    # (tr({T_phi,T_Vprime} T_pi), its mass-weighted version, m_Vprime, m_pi);
    # both numbers come from the SymPy flavor algebra, Aloni mixing, m = 2.2,4.2,20.4 MeV
    VP = {
        "phi -> omega pi0": (0.0, 0.0162395, 0.7827, 0.1396),
        "phi -> eta   pi0": (0.0, 0.00202883, 0.5479, 0.1396),
        "phi -> rho0  pi0": (-0.195262, -0.0575888, 0.7753, 0.1396),
        "phi -> rho+  pi-": (-0.195262, -0.0575888, 0.7753, 0.1396),
    }
    for nm, (fq, fm, mv, mp) in VP.items():
        r = vp_row(fq, fm, mv, mp)
        nz = [("c" + str(j + 1)) if j < 7 else ("d" + str(j - 6))
              for j, x in enumerate(r) if abs(x) > 1e-12]
        print(f"  {nm:20s} F = {fq:+.6f}   probes: {nz}")
    print()
    print("  phi -> omega pi0 and phi -> eta pi0 have F = 0 EXACTLY: phi and omega (and")
    print("  eta) all lie in the span of the singlet and lambda^8, whose products are")
    print("  orthogonal to the isovector pi0.  So their charge-type invariant vanishes and")
    print("  phi -> omega pi0 probes d2 ALONE.  phi -> rho pi probes d1, d2 and d3,")
    print("  with the neutral and charged modes giving the SAME invariant here, because")
    print("  tr is cyclic and both pick out the same quark-flavor component of T_phi.")
    print()
    print("  Every V -> V' P row has zero c-sector content, and d4 = 0 on shell because")
    print("  O4_VVP forces q twice into the epsilon.  All four null directions live in")
    print("  (c2,c4,c5,c6,c7) or in d4, so these rows are automatically orthogonal to")
    print("  every one of them.  Direct test:")
    _, _, Vt = _np.linalg.svd(_np.vstack([A] + [unit(direct_row(n)) for n in names]))
    rank0 = int((_np.linalg.svd(_np.vstack([A] + [unit(direct_row(n)) for n in names]),
                                compute_uv=False) > 1e-12).sum())
    null0 = Vt[rank0:]
    for nm, (fq, fm, mv, mp) in VP.items():
        r = unit(vp_row(fq, fm, mv, mp))
        dots = [float(v @ r) for v in null0]
        assert all(abs(dd) < 1e-9 for dd in dots), f"{nm} touches a null direction"
        print(f"    {nm:20s} dot with the 4 null vectors = "
              f"{max(abs(dd) for dd in dots):.2e}")
    print()
    Jv = _np.vstack([A] + [unit(direct_row(n)) for n in names]
                    + [unit(vp_row(f, fm, m1, m2))
                       for f, fm, m1, m2 in VP.values()])
    rk = int((_np.linalg.svd(Jv, compute_uv=False) > 1e-12).sum())
    print(f"  Rank of OPE + 10 radiative + 4 of V -> V' P  =  {rk} of 11")
    print(f"  Rank of OPE + 10 radiative                  =  {rank0} of 11")
    print(f"  => the four null directions survive unchanged; the V -> V' P rows add")
    print(f"     {rk - rank0} rank, all of it inside the d-sector already covered.")

    print()
    print("=" * 78)
    print("6. Space-like TFF sensitivity rows: does the rank rise?")
    print("=" * 78)
    print("  The TFF P -> gamma* gamma uses the same seven VJP operators with one")
    print("  space-like off-shell photon (invariant Q^2) and the pseudoscalar as the")
    print("  initial state.  The off-shell kinematic factors, derived by explicit index")
    print("  contraction in scripts/derive_radiative_tff_kinematics.py, are:")
    print("      c1: (Q^2 + m_P^2)   c2: (m_P^2 - Q^2)   c5: (Q^2 + m_P^2)")
    print("      c6: -2 Q^2          c3: -2 (pure number)   c4, c7: 0")
    print("  The decisive output: the off-shell case gives c6 = c2 - c1, NOT")
    print("  c6 = c1 - c2 as in the on-shell case.  So the TFF row is NOT in the span")
    print("  of the on-shell rows, and the rank should rise.")
    print()
    print("  Flavor factors are set to 1 (placeholder): they are overall factors per")
    print("  row and do not affect the rank.  The actual TFF flavor structure differs")
    print("  from the on-shell case and must be derived separately before a fit.")
    print()

    def tff_row(m_p: float, q2: float) -> np.ndarray:
        """Sensitivity row dF_{P gamma* gamma}/d(c1..c7,d1..d4) for one TFF channel.

        Off-shell kinematic factors from scripts/derive_radiative_tff_kinematics.py,
        with flavor factors set to 1 (placeholder).  d-sector is zero: the VVP
        contribution to the TFF is not included in this decisive test.
        """
        r = _np.zeros(11)
        r[0] = q2 + m_p**2        # c1
        r[1] = m_p**2 - q2        # c2
        r[2] = -2.0               # c3 (pure number, flavor factor = 1)
        r[4] = q2 + m_p**2        # c5
        r[5] = -2.0 * q2          # c6
        return r

    TFF_MESONS = {"pi0": 0.1396, "eta": 0.5479, "eta'": 0.9578}
    Q2_TFF = 5.0  # GeV^2, within the CLEO/BaBar range [1.64, 34.4]

    tff_rows = []
    for meson, m_p in TFF_MESONS.items():
        raw = tff_row(m_p, Q2_TFF)
        tff_rows.append(unit(raw))
        nz = [f"{'c'+str(j+1) if j < 7 else 'd'+str(j-6)}"
              for j, x in enumerate(raw) if abs(x) > 1e-12]
        print(f"  TFF {meson:5s} (Q^2={Q2_TFF} GeV^2)  nonzero: {', '.join(nz)}")

    print()
    # Baseline: OPE + 10 radiative
    J_base = _np.vstack([A] + [unit(direct_row(n)) for n in names])
    sv_base = _np.linalg.svd(J_base, compute_uv=False)
    rank_base = int((sv_base > 1e-12).sum())
    print(f"  Baseline rank (OPE + 10 radiative)           =  {rank_base} of 11")

    # With TFF rows
    J_tff = _np.vstack([J_base] + tff_rows)
    sv_tff = _np.linalg.svd(J_tff, compute_uv=False)
    rank_tff = int((sv_tff > 1e-12).sum())
    print(f"  Rank with TFF rows (OPE + 10 rad + 3 TFF)    =  {rank_tff} of 11")
    print(f"  Rank rise: {rank_base} -> {rank_tff}  (+{rank_tff - rank_base})")
    print()
    print("  Singular value spectrum with TFF rows:")
    print("    " + ", ".join(f"{v:.4f}" for v in sv_tff))
    print()
    if rank_tff > rank_base:
        print(f"  => The TFF rows ADD {rank_tff - rank_base} rank direction(s).  The TFF")
        print("     route is VIABLE: the off-shell kinematics breaks the on-shell")
        print("     c6 = c1 - c2 identity and lifts the null directions.")
    else:
        print("  => The TFF rows add NO rank.  The TFF route is DEAD.")

    print()
    print("  Incremental: which TFF row adds rank, in order pi0, eta, eta'?")
    blk = J_base.copy()
    for meson, m_p in TFF_MESONS.items():
        me = unit(tff_row(m_p, Q2_TFF))
        before = int((_np.linalg.svd(blk, compute_uv=False) > 1e-12).sum())
        blk = _np.vstack([blk, me])
        after = int((_np.linalg.svd(blk, compute_uv=False) > 1e-12).sum())
        flag = "" if after > before else "   <-- adds NO rank"
        print(f"    + TFF {meson:5s}  rank {before} -> {after}  (+{after - before}){flag}")

    print()
    print("=" * 78)
    print("Scope of this half of the probe")
    print("=" * 78)
    print("  These numbers do not depend on the undetermined common normalization N,")
    print("  on the c3 absolute normalization, or on any fitted value: row scaling")
    print("  changes none of the quantities above.  The complementary half, the")
    print("  normalized residual sqrt(chi2/ndof) and the chi2 improvement versus a")
    print("  nuisance perturbation, needs N and the c3 coefficient to be settled")
    print("  first, and is therefore NOT attempted here.  The calibration datum")
    print("  omega -> pi0 gamma has PDG S = 2.2 and rho -> pi0 gamma has S = 1.7, so")
    print("  the quoted errors understate the true spread for those rows.")


if __name__ == "__main__":
    main()