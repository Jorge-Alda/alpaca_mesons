#!/usr/bin/env python3
"""Audit the diagonal-spurion neutral-vector traces in the a -> 3P classifier.

Invocation:
    .venv/bin/python scripts/audit_neutral_vector_classification.py

The calculation is exact, symbolic, and dimensionless flavor algebra.  It
imports the classifier under audit so that the concrete checks use its actual
pseudoscalar and vector generators.
"""

from contextlib import redirect_stdout
from io import StringIO

import sympy as sp


def show(label: str, expression: sp.Basic) -> None:
    """Print one deterministic audit step using SymPy's LaTeX renderer."""
    print(f"{label}:\n{sp.latex(expression)}\n")


def commutator(left: sp.Matrix, right: sp.Matrix) -> sp.Matrix:
    """Return an explicit matrix commutator."""
    return left * right - right * left


# Importing this script executes its own diagnostics and table generation.
# Suppress that unrelated table while retaining all import-time assertions.
with redirect_stdout(StringIO()):
    import classify_a3p_vector_exchange as classifier


print("Invocation: .venv/bin/python scripts/audit_neutral_vector_classification.py")
print("Assumptions:")
print("- C_A=diag(c_u,c_d,c_s) and D_V=diag(d_u,d_d,d_s) are generic diagonal matrices.")
print("- T_P is a fully generic 3x3 matrix; no Hermiticity, diagonality, or normalization is assumed.")
print("- Matrix multiplication and the ordinary finite-dimensional trace are used.")
print("- The algebra is exact and dimensionless; no numerical or mixing approximation enters the proof.\n")


# ---------------------------------------------------------------------------
# 1. Generic 3x3 proof
# ---------------------------------------------------------------------------
c_u, c_d, c_s = sp.symbols("c_u c_d c_s")
d_u, d_d, d_s = sp.symbols("d_u d_d d_s")
t_entries = sp.symbols("t_uu t_ud t_us t_du t_dd t_ds t_su t_sd t_ss")

C_A = sp.diag(c_u, c_d, c_s)
D_V = sp.diag(d_u, d_d, d_s)
T_P = sp.Matrix(3, 3, t_entries)
commutator_PC = sp.expand(commutator(T_P, C_A))
weighted_commutator = sp.expand(D_V * commutator_PC)
weighted_diagonal = sp.Tuple(*(weighted_commutator[index, index] for index in range(3)))
generic_trace = sp.expand(sp.trace(weighted_commutator))

show("Generic diagonal C_A", C_A)
show("Generic diagonal neutral-vector generator D_V", D_V)
show("Fully generic 3x3 pseudoscalar generator T_P", T_P)
show("Intermediate commutator [T_P,C_A]", commutator_PC)
show("Intermediate product D_V[T_P,C_A]", weighted_commutator)
show("Diagonal entries of D_V[T_P,C_A]", weighted_diagonal)
show("Final generic trace Tr(D_V[T_P,C_A])", generic_trace)

assert commutator_PC.diagonal() == sp.zeros(1, 3)
assert weighted_diagonal == sp.Tuple(0, 0, 0)
assert generic_trace == 0


# ---------------------------------------------------------------------------
# 2. Exact classifier-generator audit for diagonal C_A
# ---------------------------------------------------------------------------
neutral_vectors = ("rho0", "omega", "phi")
pseudoscalar_names = tuple(classifier.PSEUDOSCALARS)
first_trace_table = sp.Matrix(
    [
        [
            classifier.simplified(
                sp.trace(
                    classifier.VECTORS[vector_name].generator
                    * commutator(
                        classifier.PSEUDOSCALARS[pseudoscalar_name].generator,
                        classifier.C_DIAGONAL,
                    )
                )
            )
            for pseudoscalar_name in pseudoscalar_names
        ]
        for vector_name in neutral_vectors
    ]
)

print("Concrete trace-table ordering:")
print(f"- rows: {', '.join(neutral_vectors)}")
print(f"- columns: {', '.join(pseudoscalar_names)}")
show("Classifier first traces Tr(T_V[T_P,C_A])", first_trace_table)

assert classifier.C_DIAGONAL == C_A
assert all(classifier.VECTORS[name].generator.is_diagonal() for name in neutral_vectors)
assert first_trace_table == sp.zeros(len(neutral_vectors), len(pseudoscalar_names))

# Also inspect the classifier's retained products, not only the first traces.
retained_neutral = tuple(
    (channel, contribution.vector)
    for (scenario, channel), contributions in classifier.all_results.items()
    if scenario == "diag C_A"
    for contribution in contributions
    if contribution.vector in neutral_vectors
)
show("Number of retained neutral-vector products in the diagonal-C_A scenario", sp.Integer(len(retained_neutral)))
assert retained_neutral == ()


# ---------------------------------------------------------------------------
# 3. Counterexamples to an over-broad all-vector zero
# ---------------------------------------------------------------------------
nonneutral_examples = (
    ("rho-", "pi+", c_d - c_u),
    ("K*-", "K+", c_s - c_u),
    ("K*0bar", "K0", c_s - c_d),
)
for vector_name, pseudoscalar_name, expected in nonneutral_examples:
    example_trace = classifier.simplified(
        sp.trace(
            classifier.VECTORS[vector_name].generator
            * commutator(
                classifier.PSEUDOSCALARS[pseudoscalar_name].generator,
                classifier.C_DIAGONAL,
            )
        )
    )
    show(
        f"Nonneutral example Tr(T_{vector_name}[T_{pseudoscalar_name},C_A])",
        example_trace,
    )
    assert sp.simplify(example_trace - expected) == 0
    assert example_trace != 0


print("Diagnostic assertions: all passed.")
print("Exact conclusion: for generic diagonal C_A and any diagonal D_V, Tr(D_V[T_P,C_A])=0 for every 3x3 T_P.")
print("Classifier conclusion: rho0, omega, and phi first traces vanish for all nine classifier pseudoscalars under diagonal C_A, and no such products are retained.")
print("Scope check: charged rho and strange K* first traces can instead be proportional to flavor differences and need not vanish.")
print("Detected classifier issue: none in the audited neutral-vector classification; ideal omega-phi mixing is stronger than needed for this zero.")
