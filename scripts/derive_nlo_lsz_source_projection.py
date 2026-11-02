"""Derive the local NLO invariant source for general flavor matrices.

Run .venv/bin/python scripts/derive_nlo_lsz_source_projection.py [--test].
Natural units; full U(3), arbitrary Hermitian couplings and quark spurion.
Retain local LECs linearly; omit loops and O(1/f_a**2). Products do not commute.
"""
import sys
import unittest
import sympy as sp


def canonical_trace(trace_expression):
    """Use trace transpose and cyclicity, absent from SymPy simplification."""
    argument = trace_expression.arg
    if argument.has(sp.Transpose):
        transposed_argument = argument.T.doit()
        if not transposed_argument.has(sp.Transpose):
            argument = transposed_argument
    factors = []
    for factor in sp.MatMul.make_args(argument):
        if isinstance(factor, sp.MatPow) and factor.exp.is_Integer and factor.exp > 0:
            factors.extend([factor.base] * int(factor.exp))
        else:
            factors.append(factor)
    rotations = [sp.MatMul(*(factors[index:]+factors[:index])) for index in range(len(factors))]
    return sp.Trace(min(rotations, key=sp.srepr))


def expand_standard(expression):
    """Distribute standard expressions and normalize the monomial traces."""
    expanded = sp.expand(sp.expand(expression).doit()).doit()
    return sp.expand(expanded.replace(lambda item: isinstance(item, sp.Trace), canonical_trace))


def trace_gradient(expression, variable):
    """Differentiate terms separately: this SymPy version loses scalar Adds.

    Each derivative uses SymPy's MatrixSymbol derivative; the transpose
    converts the result to the coefficient in tr(variable * coefficient).
    """
    terms = sp.Add.make_args(sp.expand(expression))
    return sum((sp.diff(term, variable).T for term in terms), sp.ZeroMatrix(3, 3))


def derive():
    f_pi, f_a, m_zero_squared = sp.symbols("f_pi f_a m_0^2", nonzero=True, real=True)
    c_g, c_ga = sp.symbols("c_G c_GA", real=True)
    phi = sp.MatrixSymbol("Phi", 3, 3)
    derivative_phi = sp.MatrixSymbol("partial_Phi", 3, 3)
    chi = sp.MatrixSymbol("chi_0", 3, 3)
    axial = sp.MatrixSymbol("c_q^A", 3, 3)
    pseudoscalar = sp.MatrixSymbol("c_q^P", 3, 3)
    invariant_ap = sp.MatrixSymbol("c_q^{AP}", 3, 3)
    trace = sp.Trace
    # Arbitrary-matrix operators from memory/nlo/quadratic-lagrangian.md.
    meson_kinetic = {
        "LO": trace(derivative_phi**2)/4,
        "L_5": 2*trace(chi*derivative_phi**2)/f_pi**2,
        "Lambda_1": trace(derivative_phi)**2/12,
        "L_4": 2*trace(derivative_phi**2)*trace(chi)/f_pi**2,
        "L_18": 2*trace(chi*derivative_phi)*trace(derivative_phi)/f_pi**2,
    }
    mixed_kinetic = {
        "LO": -f_pi*trace(axial*derivative_phi)/f_a,
        "L_5": -4*(trace(axial*chi*derivative_phi)+trace(axial*derivative_phi*chi))/(f_a*f_pi),
        "Lambda_1": -f_pi*trace(axial)*trace(derivative_phi)/(3*f_a),
        "L_4": -8*trace(axial*derivative_phi)*trace(chi)/(f_a*f_pi),
        "L_18": -4*(trace(axial*chi)*trace(derivative_phi)+trace(chi*derivative_phi)*trace(axial))/(f_a*f_pi),
    }
    meson_potential = {
        "LO": -trace(phi**2*chi)/4-m_zero_squared*trace(phi)**2/12,
        "L_8": -2*(trace(phi**2*chi**2)+trace(phi*chi*phi*chi))/f_pi**2,
        "Lambda_2": -trace(phi*chi)*trace(phi)/6,
        "L_6": -4*trace(phi**2*chi)*trace(chi)/f_pi**2,
        "L_7": -4*trace(phi*chi)**2/f_pi**2,
        "L_25": 4*trace(phi*chi**2)*trace(phi)/f_pi**2,
    }
    mixed_potential = {
        "LO": f_pi*trace(phi*(pseudoscalar*chi+chi*pseudoscalar))/(2*f_a)-f_pi*m_zero_squared*c_g*trace(phi)/(3*f_a),
        "L_8": 4*(trace(phi*chi**2*pseudoscalar)+trace(phi*pseudoscalar*chi**2)+2*trace(phi*chi*pseudoscalar*chi))/(f_a*f_pi),
        "Lambda_2": f_pi*(-c_g*trace(phi*chi)+trace(pseudoscalar*chi)*trace(phi))/(3*f_a),
        "L_6": 8*(trace(phi*pseudoscalar*chi)+trace(phi*chi*pseudoscalar))*trace(chi)/(f_a*f_pi),
        "L_7": 16*trace(phi*chi)*trace(pseudoscalar*chi)/(f_a*f_pi),
        "L_25": 8*(c_g*trace(phi*chi**2)-trace(pseudoscalar*chi**2)*trace(phi))/(f_a*f_pi),
    }
    # For tr(Phi G), diff(..., Phi).T is G. U(3) completeness gives
    # Phi(r_A)=2*f_pi*cA/f_a. The quadratic gradient evaluated there is
    # its Hessian applied to r_A, so this derives -M^2*r_A-mu.
    axial_flavor_shift = 2*f_pi*axial/f_a
    kinetic_sums = {}
    for name, expression in meson_kinetic.items():
        kinetic_sums[name] = expand_standard(trace_gradient(expression, derivative_phi).subs(derivative_phi, axial_flavor_shift)+trace_gradient(mixed_kinetic[name], derivative_phi))
    source_substitution = {pseudoscalar: axial-invariant_ap, c_g: c_ga-trace(axial)}
    potential_gradients, shifted_potential_gradients, invariant_flavor_sources = {}, {}, {}
    for name, expression in meson_potential.items():
        shifted_potential_gradients[name] = expand_standard(trace_gradient(expression, phi).subs(phi, axial_flavor_shift))
        potential_gradients[name] = expand_standard(trace_gradient(mixed_potential[name], phi))
        invariant_flavor_sources[name] = expand_standard((f_a/f_pi)*(shifted_potential_gradients[name]+potential_gradients[name]).subs(source_substitution))
    return locals()


class SourceProjectionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = derive()

    def test_kinetic_terms_obey_k_equals_minus_K_r_A(self):
        for expression in self.result["kinetic_sums"].values():
            self.assertEqual(expression, sp.ZeroMatrix(3, 3))

    def test_no_axial_matrix_survives_in_invariant_sources(self):
        for expression in self.result["invariant_flavor_sources"].values():
            self.assertNotIn(self.result["axial"], expression.free_symbols)

    def test_ordered_general_formulas(self):
        r = self.result
        chi, cap, f_pi, c_ga = r["chi"], r["invariant_ap"], r["f_pi"], r["c_ga"]
        expected = {
            "LO": -(chi*cap+cap*chi)/2-r["m_zero_squared"]*c_ga*sp.Identity(3)/3,
            "L_8": -4*(chi**2*cap+2*chi*cap*chi+cap*chi**2)/f_pi**2,
            "Lambda_2": -(c_ga*chi+sp.Trace(cap*chi)*sp.Identity(3))/3,
            "L_6": -8*sp.Trace(chi)*(chi*cap+cap*chi)/f_pi**2,
            "L_7": -16*sp.Trace(cap*chi)*chi/f_pi**2,
            "L_25": 8*(c_ga*chi**2+sp.Trace(cap*chi**2)*sp.Identity(3))/f_pi**2,
        }
        for name, expression in expected.items():
            self.assertEqual(expand_standard(r["invariant_flavor_sources"][name]-expression), sp.ZeroMatrix(3, 3))

    def test_nondiagonal_hermitian_example_rejects_old_shortcuts(self):
        r = self.result
        chi = sp.diag(2, 3, 5)
        cap = sp.Matrix([[1, 2+sp.I, 3], [2-sp.I, 4, 1-sp.I], [3, 1+sp.I, 2]])
        substitution = {r["chi"]: chi, r["invariant_ap"]: cap, r["f_pi"]: 1, r["c_ga"]: 2, r["m_zero_squared"]: 7}
        sources = {name: expression.subs(substitution).doit().as_explicit() for name, expression in r["invariant_flavor_sources"].items()}
        for expression in sources.values():
            self.assertEqual(sp.simplify(expression-expression.H), sp.zeros(3))
        self.assertNotEqual(sources["LO"], -chi*cap-sp.Rational(14, 3)*sp.eye(3))
        self.assertNotEqual(sources["L_8"], -16*chi**2*cap)
        self.assertNotEqual(sources["L_6"], -16*sp.trace(chi)*chi*cap)


def main():
    result = derive()
    print("Natural units; general noncommuting U(3) flavor matrices; local LECs linearly; no loops.")
    for name, expression in result["kinetic_sums"].items():
        print(f"Mixed plus shifted meson kinetic gradient, {name}:\n{sp.latex(expression)}\n")
    for name, expression in result["invariant_flavor_sources"].items():
        print(f"Mixed potential trace gradient, {name}:\n{sp.latex(result['potential_gradients'][name])}\n")
        print(f"Shifted meson potential trace gradient, {name}:\n{sp.latex(result['shifted_potential_gradients'][name])}\n")
        print(f"Invariant flavor source, {name}:\n{sp.latex(expression)}\n")


if __name__ == "__main__":
    if "--test" in sys.argv:
        unittest.main(argv=[sys.argv[0]])
    else:
        main()
