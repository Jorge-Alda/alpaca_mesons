"""Operator-level Lorentz structure of P -> gamma* gamma (space-like TFF) in RChT.

Question.  The transition form factor F_{P gamma* gamma}(Q^2) is the same VJP vertex
with one photon off shell (space-like, invariant Q^2) and the pseudoscalar as the
initial state.  To evaluate any TFF channel one needs, for each of the seven VJP
operators, the kinematic factor multiplying the scalar G.  This script derives those
factors by explicit index contraction, exactly as
scripts/derive_radiative_sm_kinematics_probe.py did for the on-shell case.

The decisive question: does the off-shell kinematics break the on-shell fixed-ratio
degeneracy in (c2, c5, c6)?  The on-shell rows all satisfy the identity
    c6 = c1 - c2          (in the c-sector, published sign convention)
because  -2 M^2 = -(M^2-m^2) - (M^2+m^2).  If the off-shell factors do NOT satisfy
this identity, the TFF sensitivity row is not in the span of the on-shell rows, and
the rank of the Jacobian rises.

Conventions.  Natural units, metric (+,-,-,-), eps^{0123} = +1.  The "vector" in the
VJP operators is the OFF-SHELL photon (momentum ell, polarization eps_ell,
ell^2 = -Q^2); the "photon" f_+ is the REAL photon (momentum k, polarization eps_k,
k^2 = 0); the "pseudoscalar" u is P (momentum p, p^2 = m_P^2).  The off-shell photon
field tensor is V^{mu nu} = ell^mu eps_ell^nu - ell^nu eps_ell^mu (NO 1/M factor,
because the off-shell photon is not a resonance with a pole mass).  The real photon
field tensor is F^{mu nu} = k^mu eps_k^nu - k^nu eps_k^mu.

Kinematics.  Pseudoscalar rest frame:
    p = (m_P, 0, 0, 0)
    k = (E, 0, 0, E),   E = (m_P^2 + Q^2)/(2 m_P)
    ell = p - k = ((Q^2-m_P^2)/(2 m_P), 0, 0, -(m_P^2+Q^2)/(2 m_P))
    ell^2 = -Q^2  (space-like)
The on-shell relation r = q - k does NOT hold: the off-shell photon has ell^2 != 0.

Reference structure.  The TFF amplitude is proportional to
    eps(ell, k, eps_ell, eps_k)  =  eps(p, k, eps_ell, eps_k)
(the two are equal because ell = p - k and eps(k,k,.,.) = 0).  All factors are
reported in units of this reference structure.

Units.  m_P and Q^2 are masses / mass-squared in any common unit; results are
polynomials in m_P^2 and Q^2.  Deterministic, no random numbers.

Invocation:
    .venv/bin/python scripts/derive_radiative_tff_kinematics.py
"""

from __future__ import annotations

import sympy as sp
try:
    from scripts.symbolic_utils import (
        levi_civita_contract,
        minkowski_dot,
    )
except ModuleNotFoundError:  # direct ``python scripts/...py`` invocation
    from symbolic_utils import (
        levi_civita_contract,
        minkowski_dot,
    )

mP, Q2 = sp.symbols("m_P Q^2", positive=True)
vx, vy = sp.symbols("vx vy", real=True)
g1, g2, g3 = sp.symbols("g1 g2 g3", real=True)

# --- off-shell kinematics, exact --------------------------------------------
E = (mP**2 + Q2) / (2 * mP)
P = (mP, 0, 0, 0)                              # pseudoscalar momentum, p^2 = m_P^2
K = (E, 0, 0, E)                               # real photon, k^2 = 0
ELL = (sp.simplify(mP - E), 0, 0, sp.simplify(-E))   # off-shell photon, ell = p - k
# ell^2 = (m_P - E)^2 - E^2 = m_P^2 - 2 m_P E = m_P^2 - (m_P^2 + Q^2) = -Q^2

# Polarizations.
# Real photon: k.eps_k = E*eps_k^0 - E*eps_k^3 = 0  =>  eps_k^0 = eps_k^3.
EPS_K = (g3, g1, g2, g3)
# Off-shell photon: ell.eps_ell = 0.  Parameterize with free symbols and impose
# the transversality constraint by solving for the time component.
v0, v3 = sp.symbols("v0 v3", real=True)
EPS_L = (v0, vx, vy, v3)
# ell.eps_ell = ELL[0]*v0 - ELL[3]*v3 = (m_P-E)*v0 + E*v3 = 0
# => v0 = -E/(m_P - E) * v3 = -(m_P^2+Q^2)/(Q^2-m_P^2) * v3
v0_sol = sp.simplify(-E / (mP - E) * v3)
EPS_L = (v0_sol, vx, vy, v3)

#: reference Lorentz structure, in units of which all factors are reported
STFF = None


dot = minkowski_dot
levi = levi_civita_contract


def S(a, b, c, d):
    return sp.cancel(sp.expand(levi(a, b, c, d)))


OBJ = {"q": ELL, "r": K, "k": P, "eV": EPS_L, "eg": EPS_K}


def term(slots, internal=()):
    """Free-slot structure times contracted pairs, in units of the reference S."""
    val = S(*[OBJ[s] for s in slots])
    for lhs, rhs in internal:
        val = val * dot(OBJ[lhs], OBJ[rhs])
    return sp.factor(sp.cancel(val / STFF))


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
    global STFF
    print("=" * 78)
    print("0. Shell conditions, exact")
    print("=" * 78)
    e = [tuple(1 if i == j else 0 for i in range(4)) for j in range(4)]
    assert levi(*e) == 1, "eps^{0123} must be +1"
    print("  eps^{0123} = +1, sanity check passed")
    print(f"  p^2 = {sp.simplify(dot(P, P))}   k^2 = {sp.simplify(dot(K, K))}   "
          f"ell^2 = {sp.simplify(dot(ELL, ELL))}")
    assert sp.simplify(dot(P, P) - mP**2) == 0
    assert sp.simplify(dot(K, K)) == 0
    assert sp.simplify(dot(ELL, ELL) + Q2) == 0, "ell^2 must equal -Q^2"
    print(f"  p - k - ell = {sp.simplify(sum(P[i] - K[i] - ELL[i] for i in range(4)))}")
    print(f"  ell.eps_ell = {sp.simplify(dot(ELL, EPS_L))}   "
          f"k.eps_k = {sp.simplify(dot(K, EPS_K))}")
    assert sp.simplify(dot(ELL, EPS_L)) == 0
    assert sp.simplify(dot(K, EPS_K)) == 0

    PK, ELKP, ELK = dot(P, K), dot(ELL, P), dot(ELL, K)
    print(f"  p.k = {sp.simplify(PK)}   ell.p = {sp.simplify(ELKP)}   "
          f"ell.k = {sp.simplify(ELK)}")
    print(f"  ell^2 = {sp.simplify(dot(ELL, ELL))}   (must be -Q^2)")

    STFF = S(P, K, EPS_L, EPS_K)
    print()
    print("=" * 78)
    print("1. The reference structure, and the permutation check")
    print("=" * 78)
    print(f"  S_TFF = eps(p,k,eps_ell,eps_k) = "
          f"{sp.factor(sp.cancel(STFF / (vx * vy)))} (vx vy)")
    assert STFF != 0, "reference structure must be nonzero"
    print("  It carries an explicit (m_P^2 - Q^2) factor (from the off-shell photon")
    print("  polarization constraint), so Q^2 dependence is present in the STRUCTURE.")

    print()
    print("  All 24 orderings of (p,k,eps_ell,eps_k) must give +/- S_TFF:")
    objs = (P, K, EPS_L, EPS_K)
    allok = True
    for perm in permutations(range(4)):
        val = S(*[objs[j] for j in perm])
        inv = sum(1 for i in range(4) for j in range(i + 1, 4) if perm[i] > perm[j])
        if sp.simplify(val - (-1) ** inv * STFF) != 0:
            allok = False
    print(f"    all 24 consistent: {allok}")
    assert allok, "epsilon contraction is not permutation consistent"

    # Verify eps(ell, k, eps_ell, eps_k) == eps(p, k, eps_ell, eps_k)
    SA = S(ELL, K, EPS_L, EPS_K)
    print(f"    eps(ell,k,eps_ell,eps_k) == S_TFF ? "
          f"{sp.simplify(SA - STFF) == 0}")
    assert sp.simplify(SA - STFF) == 0, "eps(ell,k,..) must equal eps(p,k,..)"

    print()
    print("=" * 78)
    print("2. Per-operator kinematic factors, in units of S_TFF")
    print("=" * 78)
    print("  The {V,f+} anticommutator contributes a factor 2 because the two orderings")
    print("  of the same triple matrix product are equal by cyclicity of the trace.")
    print("  Slot mapping: q->ell (off-shell photon), eV->eps_ell, r->k (real photon),")
    print("  eg->eps_k, k->p (pseudoscalar).")
    rows = {
        "O1  internal (r.k)": 2 * term(("q", "eV", "eg", "k"), (("r", "k"),)),
        "O1  alt internal (eg.k)": 2 * term(("q", "eV", "r", "k"), (("eg", "k"),)),
        "O2  internal (q.k)": two_slot_field(("eV", "k", "r", "eg"), (("q", "k"),),
                                             "r", "eg"),
        "O2  alt internal (eV.k)": two_slot_field(("q", "k", "r", "eg"), (("eV", "k"),),
                                                   "r", "eg"),
        "O5  internal (r.q)": 2 * term(("q", "eV", "eg", "k"), (("r", "q"),)),
        "O5  alt internal (eg.q)": 2 * term(("q", "eV", "r", "k"), (("eg", "q"),)),
        "O6  internal (q^2)": two_slot_field(("eV", "k", "r", "eg"), (("q", "q"),),
                                               "r", "eg"),
        "O6  alt internal (eV.q)": two_slot_field(("q", "k", "r", "eg"), (("eV", "q"),),
                                                     "r", "eg"),
        "O3  chi_- scalar": (term(("q", "eV", "r", "eg"))
                              - term(("q", "eV", "eg", "r"))),
    }
    for label, val in rows.items():
        print(f"  {label:28s} = {val}")
    print()

    # The decisive question: does the off-shell kinematics break the on-shell
    # identity c6 = c1 - c2?
    print("=" * 78)
    print("3. Decisive test: does the off-shell case break c6 = c1 - c2?")
    print("=" * 78)
    c1 = rows["O1  internal (r.k)"]
    c2 = rows["O2  internal (q.k)"]
    c5 = rows["O5  internal (r.q)"]
    c6 = rows["O6  internal (q^2)"]
    print(f"  c1 factor = {sp.simplify(c1)}")
    print(f"  c2 factor = {sp.simplify(c2)}")
    print(f"  c5 factor = {sp.simplify(c5)}")
    print(f"  c6 factor = {sp.simplify(c6)}")
    print()
    lhs = sp.simplify(c6)
    rhs = sp.simplify(c1 - c2)
    print(f"  c1 - c2 = {rhs}")
    print(f"  c6      = {lhs}")
    print(f"  c6 == c1 - c2 ?  {sp.simplify(lhs - rhs) == 0}")
    print(f"  c6 == c2 - c1 ?  {sp.simplify(lhs + rhs) == 0}")
    print()
    if sp.simplify(lhs - rhs) != 0:
        print("  => The off-shell kinematics BREAKS the on-shell identity c6 = c1 - c2.")
        print("     The TFF sensitivity row is NOT in the span of the on-shell rows.")
        print("     The rank of the Jacobian should RISE.")
    else:
        print("  => The off-shell kinematics PRESERVES the on-shell identity c6 = c1 - c2.")
        print("     The TFF sensitivity row MAY be in the span of the on-shell rows.")
        print("     The rank of the Jacobian may NOT rise.")

    print()
    print("=" * 78)
    print("4. Selection rule: alternative assignments")
    print("=" * 78)
    print("  Alternative assignments that are NOT identically zero would leave the")
    print("  operator ambiguous.  They were all checked:")
    for alt in ("O1  alt internal (eg.k)", "O2  alt internal (eV.k)",
                "O5  alt internal (eg.q)", "O6  alt internal (eV.q)"):
        val = rows[alt]
        status = "VANISHES" if sp.simplify(val) == 0 else f"= {sp.simplify(val)}"
        print(f"    {alt:30s} {status}")
    print()
    print("  The on-shell selection rule relied on eps(q,eV,r,k)=0 (q=r+k) and")
    print("  eps_V.q=0 (transversality).  For the off-shell case these do NOT hold,")
    print("  so the selection rule may be different.  The surviving contracted pairs")
    print("  are: (r.k) for O1, (q.k) for O2, (r.q) for O5, (q^2) for O6.")
    print("All assertions passed.")


if __name__ == "__main__":
    main()
