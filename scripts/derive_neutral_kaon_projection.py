"""Derive flavor-to-mass-eigenstate neutral-kaon amplitude projections.

Invocation:
    .venv/bin/python scripts/derive_neutral_kaon_projection.py

The convention is |K_S> = p|K0> + q|K0bar> and
|K_L> = p|K0> - q|K0bar>.  No CP assumption is made until the explicitly
labelled CP-limit substitution p = q = 1/sqrt(2).
"""

import sympy as sp
try:
    from scripts.symbolic_utils import show
except ModuleNotFoundError:  # direct ``python scripts/...py`` invocation
    from symbolic_utils import show



p, q = sp.symbols("p q")
A_0, A_bar = sp.symbols(r"\mathcal{A}_0 \overline{\mathcal{A}}_0")
A_S, A_L = sp.symbols(r"\mathcal{A}_S \mathcal{A}_L")

# Columns of U are the K_S and K_L kets in the (K0, K0bar) flavor basis.
U = sp.Matrix([[p, p], [q, -q]])
A_flavor = sp.Matrix([A_0, A_bar])
A_mass = sp.simplify(U.conjugate().T * A_flavor)
A_flavor_inverse = sp.simplify(U.conjugate().T.inv() * sp.Matrix([A_S, A_L]))

show("Mass eigenstates as flavor-basis columns", U)
show("Projected amplitudes (A_S, A_L)", A_mass)
show("Inverse amplitudes (A_0, A_0bar)", A_flavor_inverse)

sqrt2 = sp.sqrt(2)
cp_limit = {p: 1 / sqrt2, q: 1 / sqrt2}
A_mass_cp = sp.simplify(A_mass.subs(cp_limit))
show("Projected amplitudes in the CP-conserving mixing limit", A_mass_cp)

# Standard epsilon_K parametrization, with real and imaginary parts kept
# explicit so SymPy knows that the normalization is real and positive.
epsilon_re, epsilon_im = sp.symbols(
    r"\operatorname{Re}\epsilon_K \operatorname{Im}\epsilon_K", real=True
)
epsilon_k = epsilon_re + sp.I * epsilon_im
normalization = sp.sqrt(2 * (1 + epsilon_re**2 + epsilon_im**2))
p_epsilon = (1 + epsilon_k) / normalization
q_epsilon = (1 - epsilon_k) / normalization
A_mass_epsilon = sp.simplify(
    A_mass.subs({p: p_epsilon, q: q_epsilon})
)
state_overlap = sp.simplify(
    sp.conjugate(p_epsilon) * p_epsilon
    - sp.conjugate(q_epsilon) * q_epsilon
)
show("p(epsilon_K)", p_epsilon)
show("q(epsilon_K)", q_epsilon)
show("Projected amplitudes including epsilon_K", A_mass_epsilon)
show("Nonorthogonality <K_S|K_L>", state_overlap)

# A common tree-level special case has one hermitian flavor coupling c and a
# common hadronic/kinematic factor F.  Hermiticity gives conjugate flavor
# coefficients, while A_0bar = A_0* additionally requires F to be real.
c = sp.Symbol("c")
F = sp.Symbol("F")
A_hermitian_pair = sp.Matrix([F * c, F * sp.conjugate(c)])
A_mass_hermitian_cp = sp.simplify(
    (U.conjugate().T * A_hermitian_pair).subs(cp_limit).expand(complex=True)
)
show(
    "CP-limit amplitudes for a common factor F and hermitian coupling pair",
    A_mass_hermitian_cp,
)

F_real = sp.Symbol("F_R", real=True)
A_mass_real_tree = sp.simplify(A_mass_hermitian_cp.subs(F, F_real))
show("Further special case with real common tree-level factor", A_mass_real_tree)

# Diagnostics: the projection follows from bra coefficients, hence p* and q*.
assert A_mass[0] == sp.conjugate(p) * A_0 + sp.conjugate(q) * A_bar
assert A_mass[1] == sp.conjugate(p) * A_0 - sp.conjugate(q) * A_bar
assert sp.simplify(U.conjugate().T * A_flavor_inverse - sp.Matrix([A_S, A_L])) == sp.zeros(2, 1)
assert sp.simplify(A_mass_cp[0] - (A_0 + A_bar) / sqrt2) == 0
assert sp.simplify(A_mass_cp[1] - (A_0 - A_bar) / sqrt2) == 0
assert sp.simplify(
    sp.conjugate(p_epsilon) * p_epsilon
    + sp.conjugate(q_epsilon) * q_epsilon
    - 1
) == 0
assert sp.simplify(
    state_overlap
    - 2 * epsilon_re / (1 + epsilon_re**2 + epsilon_im**2)
) == 0
assert sp.simplify(A_mass_real_tree[0] - sqrt2 * F_real * sp.re(c)) == 0
assert sp.simplify(A_mass_real_tree[1] - sqrt2 * sp.I * F_real * sp.im(c)) == 0
