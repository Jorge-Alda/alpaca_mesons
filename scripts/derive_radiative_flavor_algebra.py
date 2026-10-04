"""SymPy derivation of the flavor algebra for V -> P gamma radiative channels.

Question.  Ruiz-Femenia:2003jdx treats exactly one on-shell vector radiative
decay, omega -> pi0 gamma.  To use the rho and omega channels in the coupled fit
we need the operator-level amplitude of each.  This script establishes the
flavor half: which ALP spurion combinations each channel can reach.

Flavor algebra (see the odd-parity numerical kernel record).  With source
matrices T_A, T_B,

    D = tr( {T_A,T_B} Ma )                          Ma = M_a^{-1}(C)
    M = tr( {T_A,T_B} {cq + Ma, B0mq} )

    VJP: (c1..c7) ~ ( 4eD,  4eD, -8eM,  0,  4eD,  4eD,  4eD)/f_a
    VVP: (d1..d4) ~ (-2D,  4M, -2D, -2D)/f_a

so c1,c2,c5,c6,c7 and d1,d3,d4 share D, c3 and d2 share M, and c4 vanishes.
A channel therefore reaches at most two independent spurion combinations.

Conventions, taken from the on-shell radiative extract of Ruiz-Femenia:2003jdx:
octet fields enter as Phi = lambda^a phi^a / sqrt(2), i.e. the matrix of an octet
generator a is lambda^a/sqrt(2); the singlet is added as omega_1/sqrt(3) with
omega_1 = diag(1,1,1); ideal mixing is
    omega = sqrt(2/3) omega_1 + sqrt(1/3) omega_8,
    phi   = -sqrt(1/3) omega_1 + sqrt(2/3) omega_8.
Every state is then normalized to unit Tr(V^2), which is the nonet requirement.
In the (u,d,s) basis that gives, exactly,

    pi0 = diag(1,-1,0)/sqrt(2)      rho   = diag(1,-1,0)/sqrt(2)
    eta = diag(1,1,-1)/sqrt(3)       omega = diag(1,1, 0)/sqrt(2)
    phi = diag(0,0,-1)

The eta shape diag(1,1,-1) is the deliberately simplified choice; the final fit
must substitute accurate pseudoscalar matrices, and most likely accurate vector
matrices, from one common nonet parametrization. Note that the phi state carries
no 1/sqrt(2): normalizing all states by a single common factor would leave
Tr(phi^2) = 1 against Tr(rho^2) = Tr(omega^2) = 2, which is not a nonet.

No isospin limit is imposed anywhere: B0mq = diag(m_u, m_d, m_s) is kept fully
general and the spurion M_a^{-1}(C) is a general Hermitian matrix, as is
c_q^{AP}.  The kernels contract the full anticommutator {T_A, T_B}, so
off-diagonal spurion entries are picked up whenever the anticommutator has
off-diagonal entries; the kaon sources below are off-diagonal single matrix
elements, the others are diagonal.

Units.  Purely dimensionless flavor algebra; no masses enter except through
B0mq.  Deterministic, no random numbers.

Invocation:
    .venv/bin/python scripts/derive_radiative_flavor_algebra.py
"""

from __future__ import annotations

import sympy as sp

# --- meson sources (u,d,s) -----------------------------------------------
# Normalization follows the convention recorded in the memory, i.e. the
# on-shell radiative extract of Ruiz-Femenia:2003jdx: octet fields enter as
# Phi = lambda^a phi^a / sqrt(2), the singlet is added as omega_1/sqrt(3), and
# every nonet state carries unit Tr(V^2).  The eta and eta' textures are the
# ones of Aloni:2018vki (sin(theta_eta_eta') ~ -1/3, cos ~ 2 sqrt(2)/3),
# rescaled to unit Tr(V^2).
SQ2, SQ3, SQ6 = sp.sqrt(2), sp.sqrt(3), sp.sqrt(6)

PI0 = sp.diag(1, -1, 0) / SQ2
ETA = sp.diag(1, 1, -1) / SQ3
ETAPRIME = sp.diag(1, 1, 2) / SQ6

RHO = sp.diag(1, -1, 0) / SQ2
OMEGA = sp.diag(1, 1, 0) / SQ2
PHI = sp.diag(0, 0, -1)

# Kaon states use the single off-diagonal element texture of Aloni:2018vki:
# each kaon is one distinct field, given as one non-zero matrix element.  Only
# the two electric-charge conserving channels are used, K*0 -> K0 gamma and
# K*+ -> K+ gamma.  The charge operator is Q = diag(2/3, -1/3, -1/3), so an
# element T_ij carries charge Q_i - Q_j: K0 = E_23 is neutral (Q_2 - Q_3 = 0)
# and K+ = E_13 has charge +1 (Q_1 - Q_3 = 1).  Pairings that change charge, such
# as K*+ -> K- gamma with Delta(Q) = -2, are not channels and must not be
# listed.
#
# A single element squares to zero, so Tr(V^2) = 0 and the unit-trace
# requirement cannot be imposed on these matrices.  The anticommutator is
# {T*, T} = {T^T, T}, not {T, T} = 2 T^2 = 0: the transpose comes from crossing
# symmetry, so the decay and its crossed partner gamma -> K*+ K- give the same
# {T*, T} and the same D.
K0 = sp.Matrix([[0, 0, 0], [0, 0, 1], [0, 0, 0]])
KPLUS = sp.Matrix([[0, 0, 1], [0, 0, 0], [0, 0, 0]])

# Charged pions follow the same single off-diagonal element texture, with the
# negative-charge states given by the transpose: pi+ ~ E_12 and pi- ~ E_21.
# Charge check: Q_1 - Q_2 = 2/3 - (-1/3) = +1, as required.  The overall factor
# is left at unity, matching the kaon convention above; it is common to T_V and
# T_P, so it rescales D but leaves the weight pattern, and hence the redundancy
# grouping, unchanged.
PIPLUS = sp.Matrix([[0, 1, 0], [0, 0, 0], [0, 0, 0]])


def vector_source(pseudoscalar: sp.Matrix) -> sp.Matrix:
    """Flavor source of the vector meson, given the pseudoscalar source.

    The two are related by crossing symmetry, which transposes the flavor
    matrix: T* = T^T.  Substituting the same matrix on both sides instead, i.e.
    {T, T} = 2 T^2, is not a shortcut but a different and wrong object: for a
    single-element matrix it vanishes identically, and it would make a channel
    vanish whenever its crossed partner does not, contradicting crossing
    symmetry.  Diagonal sources are unaffected since (T^T)^2 = T^2 there.

    This applies to the kaon pairs, where the vector is the crossing partner of
    the pseudoscalar.  It does not apply to the nonet channels, where T_V and
    T_P are different states that are not crossing partners of each other, and
    their matrices are carried independently in CHANNELS.
    """
    return pseudoscalar.T


CHANNELS: tuple[tuple[str, sp.Matrix, sp.Matrix], ...] = (
    ("omega -> pi0 gamma", OMEGA, PI0),
    ("rho   -> pi0 gamma", RHO, PI0),
    ("omega -> eta  gamma", OMEGA, ETA),
    ("rho   -> eta  gamma", RHO, ETA),
    ("phi   -> eta  gamma", PHI, ETA),
    ("K*0   -> K0  gamma", vector_source(K0), K0),
    ("K*+   -> K+  gamma", vector_source(KPLUS), KPLUS),
    ("rho+  -> pi+  gamma", vector_source(PIPLUS), PIPLUS),
    ("eta'  -> rho0 gamma", RHO, ETAPRIME),
    ("eta'  -> omega gamma", OMEGA, ETAPRIME),
)

# --- symbolic inputs -------------------------------------------------------
C00, C11, C22 = sp.symbols("C00 C11 C22", real=True)
X01, X02, X12 = sp.symbols("X01 X02 X12", real=True)   # real parts of Ma off-diagonal
Y01, Y02, Y12 = sp.symbols("Y01 Y02 Y12", real=True)   # imaginary parts
Q00, Q11, Q22 = sp.symbols("q00 q11 q22", real=True)
M_U, M_D, M_S = sp.symbols("m_u m_d m_s", positive=True)
E, F_A = sp.symbols("e f_a", nonzero=True)

MA = sp.Matrix([
    [C00, X01 + sp.I * Y01, X02 + sp.I * Y02],
    [X01 - sp.I * Y01, C11, X12 + sp.I * Y12],
    [X02 - sp.I * Y02, X12 - sp.I * Y12, C22],
])
CQ = sp.diag(Q00, Q11, Q22)
B0MQ = sp.diag(M_U, M_D, M_S)



def anticommutator(vector: sp.Matrix, pseudoscalar: sp.Matrix) -> sp.Matrix:
    """Full matrix {T_V, T_P}.

    This must stay a full matrix.  Projecting onto the diagonal first and then
    tracing would silently discard every off-diagonal entry of M_a^{-1}(C) that
    an off-diagonal anticommutator probes, so the kernels contract the full
    product.  For diagonal meson sources the anticommutator is diagonal and the
    two agree, which is why this was not visible before.
    """
    return sp.simplify(vector * pseudoscalar + pseudoscalar * vector)


def kernels(sources: sp.Matrix) -> tuple[sp.Expr, sp.Expr]:
    """Return (D, M) for a full anticommutator {T_V, T_P}."""
    total = CQ + MA
    mass = total * B0MQ + B0MQ * total
    return sp.simplify(sp.trace(sources * MA)), sp.simplify(sp.trace(sources * mass))


def operator_signature(d: sp.Expr, m: sp.Expr) -> tuple[list, list]:
    c = [4 * E * d / F_A, 4 * E * d / F_A, -8 * E * m / F_A, sp.Integer(0),
         4 * E * d / F_A, 4 * E * d / F_A, 4 * E * d / F_A]
    v = [-2 * d / F_A, 4 * m / F_A, -2 * d / F_A, -2 * d / F_A]
    return c, v


def main() -> None:
    print("Meson sources (u,d,s); convention Tr[V^2] = 1 (memory convention)")
    for label, matrix in (("pi0", PI0), ("eta", ETA), ("rho", RHO),
                          ("omega", OMEGA), ("phi", PHI)):
        trace_of_square = sp.simplify(sp.trace(matrix * matrix))
        print(f"  {label:6s} {sp.diag(*matrix.diagonal())}   Tr(V^2) = {trace_of_square}")
        assert sp.simplify(trace_of_square - 1) == 0 or trace_of_square == 0, \
            f"{label} violates the unit Tr[V^2] = 1 convention"

    print("\n" + "=" * 78)
    print("1. Source anticommutators {T_V, T_P} (full matrices)")
    print("=" * 78)
    patterns: dict[str, sp.Matrix] = {}
    for name, vector, pseudoscalar in CHANNELS:
        pattern = anticommutator(vector, pseudoscalar)
        patterns[name] = pattern
        print(f"  {name:22s} {pattern}")

    print("\n" + "=" * 78)
    print("2. Spurion combinations, general Hermitian M_a^{-1}(C) and c_q^AP")
    print("   B0mq = diag(m_u, m_d, m_s) with no isospin assumption")
    print("=" * 78)
    for name, _, _ in CHANNELS:
        d, m = kernels(patterns[name])
        print(f"\n  {name}")
        print(f"    D = {sp.expand(d)}")
        print(f"    M = {sp.expand(m)}")
        c, v = operator_signature(d, m)
        print(f"    (c1..c7) propto {sp.factor(c[0])} ... ; c3 propto {sp.factor(c[2])} ; c4 = 0")
        print(f"    (d1..d4) propto {sp.factor(v[0])} ... ; d2 propto {sp.factor(v[1])}")

    print("\n" + "=" * 78)
    print("3. Redundancy: distinct source patterns and the channels sharing them")
    print("=" * 78)
    groups: dict[str, list[str]] = {}
    for name, pattern in patterns.items():
        weights = [sp.simplify(w) for w in pattern.diagonal()]
        nonzero = [w for w in weights if w != 0]
        scale = max(abs(w) for w in nonzero) if nonzero else sp.Integer(0)
        key = sp.srepr([sp.simplify(w / scale) for w in weights]) if scale != 0 else "ZERO"
        groups.setdefault(key, []).append(name)
    for index, (_, members) in enumerate(sorted(groups.items(), key=lambda kv: kv[1][0]), 1):
        print(f"  pattern {index}: {', '.join(members)}")
    print(f"\n  {len(CHANNELS)} channels share {len(groups)} distinct patterns.")
    for _, members in groups.items():
        if len(members) > 1:
            print(f"    REDUNDANT: {members}")

    print("\n" + "=" * 78)
    print("4. Reach of the diagonal spurion")
    print("=" * 78)
    for name, pattern in patterns.items():
        weights = [sp.simplify(w) for w in pattern.diagonal()]
        reach = [f"C{i}{i}" for i, w in enumerate(weights) if w != 0]
        source = patterns[name]
        offdiagonal = [(i, j) for i in range(3) for j in range(3)
                       if i != j and source[i, j] != 0]
        if offdiagonal:
            reach.append("off-diagonal " + ",".join(f"C{i}{j}" for i, j in offdiagonal))
        nonzero = [w for w in weights if w != 0]
        scale = max(abs(w) for w in nonzero) if nonzero else sp.Integer(0)
        shape = sp.simplify(weights[0] / scale) if scale != 0 else sp.Integer(0)
        print(f"  {name:22s} weights {weights}  shape~({shape},+-{shape},..)  reaches {reach}")

    print("\n" + "=" * 78)
    print("5. Special cases (isospin limits), for orientation only")
    print("=" * 78)
    iso = {M_U: sp.Symbol("m_l"), M_D: sp.Symbol("m_l"), M_S: sp.Symbol("m_s")}
    for name, pattern in patterns.items():
        d, m = kernels(pattern)
        print(f"  {name:22s} D = {sp.expand(d).subs(iso)}")
        print(f"  {'':22s} M = {sp.expand(m).subs(iso)}")


if __name__ == "__main__":
    main()