"""Derive the single-epsilon obstruction for a -> P1 P2 P3.

Run: .venv/bin/python scripts/audit_alp_three_p_odd_parity.py
Natural units; vacuum Lorentz invariance, spin-zero external states,
p_a = p_1 + p_2 + p_3. No on-shell mass or parity assumption is needed.
This is not an all-orders absence theorem for two anomalous insertions.
"""

import sympy as sp
from itertools import product


def show(label, expression):
    print(f"{label}:\n{sp.latex(expression)}\n")


def main():
    # Explicit components are appropriate here: the objects are four-momenta,
    # not abstract flavour matrices. The determinant is SymPy's native
    # Levi-Civita contraction of the four column vectors.
    momenta = [sp.Matrix(sp.symbols(f"p_{i}_0:4")) for i in range(1, 4)]
    pa = sp.Matrix(sp.symbols("p_a_0:4"))
    epsilon_before = sp.Matrix.hstack(pa, *momenta).det()
    show("epsilon(p_a,p_1,p_2,p_3), before conservation", epsilon_before)
    substitutions = dict(zip(pa, sum(momenta, sp.zeros(4, 1))))
    epsilon_after = sp.factor(epsilon_before.subs(substitutions))
    show("epsilon after p_a = p_1 + p_2 + p_3", epsilon_after)

    # Every tree internal momentum is a linear combination of the three
    # independent external momenta. This derives the stronger determinant
    # statement for four arbitrary such momenta, without assuming particular
    # routing or pole masses.
    basis = sp.Matrix.hstack(*momenta)
    coefficients = sp.Matrix(3, 4, sp.symbols("b_0:12"))
    routed = basis * coefficients
    # Avoid expanding a huge determinant: native column operations remove
    # the first three linear combinations once a nonsingular 3x3 minor is
    # chosen. Equivalently rank(A B) <= rank(A) <= 3. An augmented 4x4
    # determinant with the fourth momentum expressed in that basis suffices.
    fourth = basis * sp.Matrix(sp.symbols("x_1:4"))
    routed_determinant = sp.factor(sp.Matrix.hstack(*momenta, fourth).det())
    show("epsilon(p_1,p_2,p_3,x_1 p_1+x_2 p_2+x_3 p_3)", routed_determinant)
    show("general tree momentum routing matrix", routed)

    # Derive tree counting from incidence and connected-tree identities.
    # Bilinear mixings are incorporated in propagators/physical fields and
    # there are no tadpoles; every interaction vertex has valence >= 3.
    E, I, V, total_valence = sp.symbols("E I V N", integer=True)
    incidence = sp.Eq(total_valence, 2 * I + E)
    tree = sp.Eq(I, V - 1)
    show("incidence counting", incidence)
    show("connected tree identity", tree)
    external_count = sp.solve(incidence.subs(I, sp.solve(tree, I)[0]), E)[0]
    show("external count N=sum(n_v)", sp.Eq(E, external_count))
    n3, n4 = sp.symbols("n_3 n_4", integer=True, nonnegative=True)
    # E=4 bounds each valence at four and sum(n_v-2)=2. The native
    # diophantine solver derives all integer solutions, then positivity
    # selects the allowed vertex counts.
    equation = sp.expand(external_count.subs({total_valence: 3*n3+4*n4,
                                             V: n3+n4}) - 4)
    show("four-leg cubic/quartic counting equation", sp.Eq(equation, 0))
    solutions = sp.diophantine(equation, syms=(n3, n4))
    allowed = set()
    for solution in solutions:
        parameter = next(iter(set().union(*(s.free_symbols for s in solution))))
        # n_3,n_4 <= 2 follows directly from n_3+2 n_4=2.
        for value in range(-2, 3):
            evaluated = tuple(s.subs(parameter, value) for s in solution)
            if all(s.is_nonnegative for s in evaluated):
                allowed.add(evaluated)
    show("allowed counts (n_3,n_4)", sp.FiniteSet(*allowed))

    # The two-odd-cubic possibility uses only the VVP/VJP anomalous
    # vertices and their a/P source substitutions in the declared action.
    # Integer incidence vectors (a,P,V,gamma) give external multiplicities
    # after removing the two endpoints of the one internal line.
    alp_vertices = ((1,0,2,0), (1,0,1,1), (1,0,0,2))
    meson_vertices = ((0,1,2,0), (0,1,1,1), (0,1,0,2))
    survivors = set()
    for left, right in product(alp_vertices, meson_vertices):
        for species in (2, 3):
            if left[species] and right[species]:
                external = sp.Matrix(left) + sp.Matrix(right)
                external[species] -= 2
                survivors.add(tuple(external))
    show("two-odd-cubic external counts (a,P,V,gamma)", sp.FiniteSet(*survivors))
    target = (1,3,0,0)
    assert epsilon_after == 0
    assert routed_determinant == 0
    assert allowed == {(2,0), (0,1)}
    assert target not in survivors
    print("Assertions passed: single-epsilon scalar amplitude vanishes; "
          "declared odd cubic pairs cannot produce a -> 3P at tree level.")


if __name__ == "__main__":
    main()
