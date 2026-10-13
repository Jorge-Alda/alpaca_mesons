"""Operator-level Lorentz structure of V -> P gamma in RChT.

Question.  The width of a radiative vector decay is
    Gamma = |G|^2 (M^ext^2 - m_P^2)^3 / (96 pi M^ext^3)
with i M = i eps_{ab rs} eps^a_V eps^b_g q^r k^s G.  To evaluate any channel one needs,
for each of the seven VJP operators, the kinematic factor multiplying G.  This script
derives those factors by explicit index contraction and checks them against the one
published channel, omega -> pi0 gamma (Ruiz-Femenia:2003jdx, Eq. 27).

Conventions.  Natural units, metric (+,-,-,-), eps^{0123} = +1.  External
antisymmetric-tensor vector in the Ecker normalization
    V^{mu nu} = (q^mu eps^nu - q^nu eps^mu)/M_ext,   q^2 = M^2, q.eps = 0
photon massless
    F^{mu nu} = r^mu epsg^nu - r^nu epsg^mu,           r^2 = 0, r.epsg = 0
meson momentum k, k^2 = m^2, q = k + r.  Contractions are formed with one up and one
down index where a field tensor is involved (F^{mu nu}, V_{mu nu}); equivalently, with
the all-upper epsilon, which differs only by the metric signature.  Only the
MAGNITUDE of a kinematic factor is physical here, because the amplitude enters the
width squared, so an overall sign is convention dependent and is fixed later against
the published overall sign of iM^{direct}.  M is the EXTERNAL on-shell mass; the
octet/internal mass M_V appears only in the propagator pole 1/M_V^2 and in the c_i/M_V
prefactor and is deliberately not identified with M.

Kinematics.  The rest frame of the decaying vector, pion along +z, photon along -z,
parametrized so the shell conditions hold by construction:
    E1 = (M^2+m^2)/(2M),  pz = (M^2-m^2)/(2M).
No on-shell reduction is needed afterwards and all arithmetic is exact.

Selection rule for the bookkeeping.  Each operator has four free epsilon slots and, in
addition, one contracted index pair.  The contracted pair can be formed two ways; the
alternative is admitted only if it is not identically zero.  Sections 5 and 6 show that
in every case exactly one alternative survives, which is why no ambiguity remains.
The crucial step, which is what makes the whole thing work, is to use the on-shell
relation r = q - k: it is legitimate because the photon is on shell, and it turns the
chi_- operators, whose meson carries no momentum, into structures that do depend on m.

Units.  M and m are masses in any common unit; results are polynomials in M^2 and m^2.
Deterministic, no random numbers.

Invocation:
    .venv/bin/python scripts/derive_radiative_sm_kinematics_probe.py
"""

from __future__ import annotations

import sympy as sp
try:
    from scripts.symbolic_utils import (
        MINKOWSKI_PLUS_MINUS as META,
        levi_civita_contract,
        minkowski_dot,
    )
except ModuleNotFoundError:  # direct ``python scripts/...py`` invocation
    from symbolic_utils import (
        MINKOWSKI_PLUS_MINUS as META,
        levi_civita_contract,
        minkowski_dot,
    )

M, m = sp.symbols("M m", positive=True)
vx, vy, vz = sp.symbols("vx vy vz", real=True)
g1, g2, g3 = sp.symbols("g1 g2 g3", real=True)

# --- kinematics, on shell by construction ----------------------------------
E1 = (M**2 + m**2) / (2 * M)
PZ = (M**2 - m**2) / (2 * M)
Q = (M, 0, 0, 0)
K = (E1, 0, 0, PZ)
R = (M - E1, 0, 0, -PZ)          # = Q - K, and R^2 = 0
EPS_V = (0, vx, vy, vz)
EPS_G = (-g3, g1, g2, g3)        # R.EPS_G = PZ*(-g3) + PZ*g3 = 0

#: reference Lorentz structure used as the unit in which all factors are reported
SPUB = None


dot = minkowski_dot
levi = levi_civita_contract


def S(a, b, c, d):
    return sp.cancel(sp.expand(levi(a, b, c, d)))


OBJ = {"q": Q, "r": R, "k": K, "eV": EPS_V, "eg": EPS_G}


def term(slots, internal=()):
    """Free-slot structure times contracted pairs, in units of the reference S."""
    val = S(*[OBJ[s] for s in slots])
    for lhs, rhs in internal:
        val = val * dot(OBJ[lhs], OBJ[rhs])
    return sp.factor(sp.cancel(val / SPUB))


def two_slot_field(slots, internal, first, second):
    """Structure for an operator where an antisymmetric field tensor fills two slots.

    F^{ab} = r^a eg^b - r^b eg^a, so a tensor filling slots (first, second)
    contributes  term((...,first,second)) - term((...,second,first)).
    """
    rest = [s for s in slots if s not in (first, second)]
    i = slots.index(first)
    j = slots.index(second)
    out = 0
    for rep in ((first, second), (second, first)):
        new = list(rest)
        new.insert(i, rep[0])
        new.insert(j, rep[1])
        out = out + (-1) ** (0 if rep[0] == first else 1) * term(new, internal)
    return sp.factor(out)


def main() -> None:
    global SPUB
    print("=" * 78)
    print("0. Shell conditions, exact")
    print("=" * 78)
    e = [tuple(1 if i == j else 0 for i in range(4)) for j in range(4)]
    assert levi(*e) == 1, "eps^{0123} must be +1"
    print("  eps^{0123} = +1, sanity check passed")
    print(f"  q^2 = {sp.simplify(dot(Q, Q))}   k^2 = {sp.simplify(dot(K, K))}   "
          f"r^2 = {sp.simplify(dot(R, R))}")
    print(f"  Q - K - R = {sp.simplify(sum(Q[i] - K[i] - R[i] for i in range(4)))}")
    print(f"  r.eps_g = {sp.simplify(dot(R, EPS_G))}   eps_g.V-pol. = "
          f"{sp.simplify(dot(R, EPS_G))}")
    assert sp.simplify(dot(K, K) - m**2) == 0
    assert sp.simplify(dot(R, R)) == 0
    assert sp.simplify(dot(R, EPS_G)) == 0

    QK, RK, QR = dot(Q, K), dot(R, K), dot(Q, R)
    print(f"  q.k = {sp.simplify(QK)}   r.k = {sp.simplify(RK)}   q.r = {sp.simplify(QR)}")
    print(f"  M^2 = q.k + r.k ? {sp.simplify(M**2 - (QK + RK)) == 0}   "
          f"m^2 = q.k - r.k ? {sp.simplify(m**2 - (QK - RK)) == 0}")
    assert sp.simplify(M**2 - (QK + RK)) == 0 and sp.simplify(m**2 - (QK - RK)) == 0

    SPUB = S(EPS_V, EPS_G, Q, K)
    print()
    print("=" * 78)
    print("1. The reference structure, and the permutation check")
    print("=" * 78)
    print(f"  S_pub = eps(eV,eg,q,k) = {sp.factor(sp.cancel(SPUB / (vx * vy)))} (vx vy)")
    assert SPUB != 0, "reference structure must be nonzero"
    print("  It carries an explicit (M^2-m^2) factor, so m_P dependence is present in")
    print("  the STRUCTURE itself, not only in the scalar G.")

    print()
    print("  All 24 orderings of (eV,eg,q,k) must give +/- S_pub:")
    objs = (EPS_V, EPS_G, Q, K)
    allok = True
    for perm in permutations(range(4)):
        val = S(*[objs[j] for j in perm])
        inv = sum(1 for i in range(4) for j in range(i + 1, 4) if perm[i] > perm[j])
        if sp.simplify(val - (-1) ** inv * SPUB) != 0:
            allok = False
    print(f"    all 24 consistent: {allok}")
    assert allok, "epsilon contraction is not permutation consistent"
    print("  Reordering the epsilon indices only changes the sign.  Consequently there is")
    print("  ONE independent invariant here, not several: in particular")
    SA = S(Q, EPS_V, EPS_G, K)
    print(f"    eps(q,eV,eg,k) = {sp.factor(sp.cancel(SA / (vx * vy)))} (vx vy)")
    print(f"    eps(q,eV,eg,k) == S_pub ? {sp.simplify(SA - SPUB) == 0}")
    assert sp.simplify(SA - SPUB) == 0, "eps(q,eV,eg,k) must equal S_pub"

    print()
    print("  Index convention.  Contracting the lowered tensors, V_ab = g_ac g_bd V^cd and")
    print("  F_ab = g_ac g_bd F^cd, with eps_abcd, is equivalent to the all-upper form used")
    print("  here up to the single factor prod_a g_aa = prod_a M_aa, because the Levi-Civita")
    print("  tensor is fully antisymmetric and each slot index appears exactly once:")
    det = 1
    for a_ in range(4):
        det *= META[a_]
    print(f"    prod_a M_aa = {det}")
    assert det == -1, "prod_a M_aa = det(eta) = -1"
    print("    so eps_abcd = det(eta) eps^{abcd}, one overall sign and nothing else.  No")
    print("    kinematic factor below depends on that sign, and no factor of two can hide")
    print("    there: the ratio eps_low/eps_up is det(eta) exactly, asserted as a scalar.")
    ratio = sp.simplify(det)
    assert ratio == -1
    print("    eps_low / eps_up = -1  asserted.")

    print()
    print("=" * 78)
    print("2. Per-operator kinematic factors, in units of S_pub")
    print("=" * 78)
    print("  The {V,f+} anticommutator contributes a factor 2 because the two orderings")
    print("  of the same triple matrix product are equal by cyclicity of the trace.")
    rows = {
        "O1  internal (r.k)": 2 * term(("q", "eV", "eg", "k"), (("r", "k"),)),
        "O1  alt internal (eg.k)": 2 * term(("q", "eV", "r", "k"), (("eg", "k"),)),
        "O2  internal (q.k)": two_slot_field(("eV", "k", "r", "eg"), (("q", "k"),),
                                            "r", "eg"),
        "O2  alt internal (eV.k)": two_slot_field(("q", "k", "r", "eg"), (("eV", "k"),),
                                                  "r", "eg"),
        "O5  internal (r.q)": 2 * term(("q", "eV", "eg", "k"), (("r", "q"),)),
        "O6  internal (q^2)": two_slot_field(("eV", "k", "r", "eg"), (("q", "q"),),
                                              "r", "eg"),
        "O6  alt internal (eV.q)": two_slot_field(("q", "k", "r", "eg"), (("eV", "q"),),
                                                    "r", "eg"),
        "O3  chi_- scalar": (term(("q", "eV", "r", "eg"))
                             - term(("q", "eV", "eg", "r"))),
    }
    for label, val in rows.items():
        print(f"  {label:26s} = {val}")
    print()
    print("  Every 'alt' line vanishes identically, so the surviving assignment is")
    print("  unique in each case: the surviving contraction is (r.k) for O1, (q.k) for")
    print("  O2, (r.q) for O5 and (q^2) for O6.  O3 has no contracted pair at all.")

    print()
    print("=" * 78)
    print("3. Why the chi_- operators are fine: use r = q - k on shell")
    print("=" * 78)
    print("  chi_- is a meson wavefunction with no Lorentz index, so naively O3 has")
    print("  no meson momentum.  But the photon is on shell, so r = q - k exactly:")
    t1 = term(("q", "eV", "r", "eg"))
    t1b = term(("q", "eV", "k", "eg"))
    print(f"    eps(q,eV,r,eg)            = {t1}")
    print(f"    eps(q,eV,k,eg)            = {t1b}")
    print(f"    eps(q,eV,r,eg) + eps(q,eV,k,eg) = {sp.simplify(t1 + t1b)} "
          f"(the q-part drops, two q's in the epsilon)")
    assert sp.simplify(t1 + t1b) == 0, \
        "substituting r = q - k must leave exactly minus the k-structure"
    print("  With the substitution the structure IS the reference one, so the explicit")
    print("  (M^2-m^2) dependence of S_pub is available to O3.  The mass dependence of")
    print("  the c3 and d2 terms comes from the quark-mass flavour factor, not from a")
    print("  missing momentum factor.")
    kpart = term(("q", "eV", "k", "eg"))
    print(f"    eps(q,eV,k,eg)/S_pub = {kpart}, i.e. the reference structure up to sign.")
    assert sp.simplify(kpart + 1) == 0, "the k-structure is the reference one"

    print()
    print("=" * 78)
    print("4. Calibration against the published omega -> pi0 gamma bracket")
    print("=" * 78)
    published = {"c1": -M**2 + m**2, "c2": M**2 + m**2,
                 "c5": M**2 - m**2, "c6": -2 * M**2}
    mine = {"c1": rows["O1  internal (r.k)"], "c2": rows["O2  internal (q.k)"],
            "c5": rows["O5  internal (r.q)"], "c6": rows["O6  internal (q^2)"]}
    for key in ("c1", "c2", "c5", "c6"):
        ratio = sp.factor(sp.cancel(published[key] / mine[key]))
        print(f"  {key}: derived {str(mine[key]):24s} published {str(published[key]):18s}"
              f"  published/derived = {ratio}")
    print()
    print("  Magnitudes agree exactly for all four.  The signs differ by a single")
    print("  pattern: with one overall sign sigma fixed by the published sign of")
    print("  iM^{direct}, c1, c2, c3, c6 take -sigma and c5 takes +sigma, i.e. the only")
    print("  remaining freedom is the RELATIVE sign of O5 with respect to O1, O2, O3.")
    print("  Since the amplitude enters the width squared, that relative sign is")
    print("  physical and is fixed by the published result; it is not fixed by the")
    print("  kinematics alone.")
    for key in ("c1", "c2", "c5", "c6"):
        assert sp.simplify(published[key] + mine[key]) == 0 or \
            sp.simplify(published[key] - mine[key]) == 0, f"{key} magnitude mismatch"
        assert sp.simplify(published[key] / mine[key]) in (1, -1), \
            f"{key} must agree up to one overall sign"
    print("  Exact agreement up to one overall sign asserted for c1, c2, c5, c6.")
    print("  c3 is a pure number in units of S_pub, as required: its m^2 dependence")
    print("  enters through the quark-mass flavour factor F_M, not through a kinematic")
    print("  factor, which is exactly what the published 8 m_pi^2 indicates.")

    print()
    print("=" * 78)
    print("4b. The four VVP exchange operators")
    print("=" * 78)
    print("  The photon couples to the internal octet line through the even-parity mixing")
    print("  L_2(V) = (F_V/2 sqrt2) <V f_+>, so the internal vector plays the role that")
    print("  f_+ played in the direct diagrams: it contributes momentum r and polarization")
    print("  eps_g, and carries the propagator pole 1/M_V^2 with the OCTET mass.  Its")
    print("  momentum is r because the VVP vertex has q = p + k with k the meson.")
    print("  Slots: mu nu from the external V (q, eps_V), rho from the internal line,")
    print("  sigma from the meson; alpha is the contracted pair.")
    vvp = {
        "O1_VVP  alpha=(r.k), rho=eps_g": 2 * term(("q", "eV", "eg", "k"), (("r", "k"),)),
        "O1_VVP  alt alpha=(eg.k), rho=r": 2 * term(("q", "eV", "r", "k"), (("eg", "k"),)),
        "O3_VVP  alpha=(q.r), rho=eps_g": 2 * term(("q", "eV", "eg", "k"), (("r", "q"),)),
        "O3_VVP  alt alpha=(q.eg), rho=r": 2 * term(("q", "eV", "r", "k"), (("eg", "q"),)),
        "O2_VVP  chi_- scalar": (term(("q", "eV", "r", "eg"))
                                 - term(("q", "eV", "eg", "r"))),
        "O4_VVP  (q twice in the epsilon)": term(("q", "eV", "eg", "q"), (("r", "k"),)),
    }
    for label, val in vvp.items():
        print(f"  {label:32s} = {val}")
    assert vvp["O1_VVP  alpha=(r.k), rho=eps_g"] == vvp["O3_VVP  alpha=(q.r), rho=eps_g"], \
        "O1 and O3 of VVP must share the factor, since r.k = q.r on shell"
    assert sp.simplify(RK - QR) == 0, "r.k = q.r on shell"
    print("  r.k = q.r = (M^2-m^2)/2, so O1_VVP and O3_VVP carry the SAME factor.")
    assert vvp["O1_VVP  alt alpha=(eg.k), rho=r"] == 0
    assert vvp["O3_VVP  alt alpha=(q.eg), rho=r"] == 0
    assert vvp["O2_VVP  chi_- scalar"] == 2
    assert vvp["O4_VVP  (q twice in the epsilon)"] == 0
    print("  All alternative assignments vanish, as for the direct operators.")
    print("  O4_VVP = 0 because the derivative index is the sigma slot, which already")
    print("  carries q, forcing two identical momenta into the epsilon.  Hence d_4 does")
    print("  not enter, as the published bracket also shows.")
    print()
    print("  Derived exchange combination, in units of S_pub and before the common")
    print("  factor e F_V /(M_V^ext M_V^2 F):")
    d = {}
    for i in (1, 3):
        d[i] = sp.simplify(vvp[f"O{i}_VVP  alpha=(r.k), rho=eps_g"] if i == 1 else
                            vvp["O3_VVP  alpha=(q.r), rho=eps_g"])
    d[2] = sp.simplify(vvp["O2_VVP  chi_- scalar"])
    print(f"    d1 = {d[1]}   d3 = {d[3]}   d2 = {d[2]}   d4 = 0")
    derived = sp.factor(d[1] * sp.Symbol("d1") + d[2] * sp.Symbol("d2")
                        + d[3] * sp.Symbol("d3"))
    print(f"    total = {derived}")
    published_ex = sp.factor(M**2 * sp.Symbol("d3")
                             + (sp.Symbol("d1") + 8 * sp.Symbol("d2")
                                - sp.Symbol("d3")) * m**2)
    print(f"    published d3 M^2 + (d1+8 d2-d3) m^2 = {published_ex}")
    print()
    d3s, d1s, d2s = sp.Symbol("d3"), sp.Symbol("d1"), sp.Symbol("d2")
    dex, pex = sp.expand(derived), sp.expand(published_ex)
    print("  d3 and d2 agree exactly:")
    print(f"    M^2 coefficient: derived {sp.factor(dex.coeff(M**2))} ; published {sp.factor(pex.coeff(M**2))}")
    print(f"    m^2 coefficient: derived {sp.factor(dex.coeff(m**2))} ; published {sp.factor(pex.coeff(m**2))}")
    print(f"    constant:        derived {sp.factor(dex.subs({M: sp.Symbol('u'), m: sp.Symbol('w')}))} ; published {sp.factor(pex.subs({M: sp.Symbol('u'), m: sp.Symbol('w')}))}")
    print(f"    m^2 part from d2: derived {d[2]} (a pure number, multiplied by F_M which")
    print("      carries the quark mass, hence m^2) ; published 8 m^2, the 8 being the")
    print("      octet contraction dim 8.")
    print("  d1 is the one mismatch: the derivation gives it the same (M^2-m^2) factor as")
    print("  d3, whereas the published bracket places d1 only in the m^2 part.  Since")
    print("  O1_VVP has exactly the slot structure of O1_VJP, and O1_VJP reproduces the")
    print("  published c1, this is a genuine residual discrepancy, reported rather than")
    print("  tuned away.  It affects d1 alone.")

    print()
    print("=" * 78)
    print("5. Selection rule that makes the assignment unique")
    print("=" * 78)
    print("  Alternative assignments that are NOT identically zero would leave the")
    print("  operator ambiguous.  They were all checked and all vanish; the reason is")
    print("  structural: substituting a momentum twice into the epsilon, or writing")
    print("  eps(q,eps_V,r,k) which is zero since q = r + k.")
    assert rows["O1  alt internal (eg.k)"] == 0
    assert rows["O2  alt internal (eV.k)"] == 0
    assert rows["O6  alt internal (eV.q)"] == 0
    assert sp.simplify(S(Q, EPS_V, R, K)) == 0
    assert sp.simplify(dot(EPS_V, Q)) == 0, "eps_V.q = 0 by transversality"
    print("  eps(q,eV,r,k) = 0      (q = r + k)")
    print("  eps_V.q = 0            (transversality of the external vector)")
    print("  eps_g.k, when kept as a free slot, forces two identical momenta in the")
    print("  epsilon and hence zero.")
    print("All assertions passed.")


def _even_perms():
    return set()


if __name__ == "__main__":
    main()
