"""Derive finite chiral transformations of higher ALP quark sources.

Run: .venv/bin/python scripts/derive_alp_higher_quark_bilinears.py
Natural units: a and f_a have mass dimension one. The differentiation
variable is a/f_a; divide an nth derivative by f_a**n for d**n/da**n.
All flavour matrices are noncommuting. The source convention is
chi' = g_R chi g_L^dagger and r' = g_R r g_R^dagger
+ i g_R partial_mu g_R^dagger, with the analogous left source.
The displayed Taylor polynomials retain degree four in a/f_a.
This supplies coefficients at a=0 through that degree without assuming
that the untransformed source terminates at first order.
"""

import sympy as sp


def show(label, expression):
    print(f"{label}:\n{sp.latex(expression)}\n")


def exponential_polynomial(generator, variable, order):
    """Standard exponential Taylor polynomial with noncommuting generator."""
    return sp.Add(*(variable**power * generator**power / sp.factorial(power)
                    for power in range(order + 1)))


def derivative_coefficients(expression, variable, order):
    """Differentiate the product before evaluating at the origin."""
    return tuple(sp.expand(sp.diff(expression, variable, degree).subs(variable, 0))
                 for degree in range(order + 1))


def derive(order=4):
    variable = sp.Symbol(r"\left(\frac{a}{f_a}\right)", real=True)
    right_generator = sp.Symbol(r"(\delta+\kappa)", commutative=False)
    left_generator = sp.Symbol(r"(\delta-\kappa)", commutative=False)
    mass_derivatives = sp.symbols(f"chi_0:{order + 1}", commutative=False)
    right_current_derivatives = sp.symbols(f"c_R_0:{order}", commutative=False)
    left_current_derivatives = sp.symbols(f"c_L_0:{order}", commutative=False)
    mass_source = sp.Add(*(mass_derivatives[degree] * variable**degree / sp.factorial(degree)
                          for degree in range(order + 1)))
    right_current = sp.Add(*(right_current_derivatives[degree] * variable**degree / sp.factorial(degree)
                            for degree in range(order)))
    left_current = sp.Add(*(left_current_derivatives[degree] * variable**degree / sp.factorial(degree)
                           for degree in range(order)))
    right_group = exponential_polynomial(sp.I * right_generator, variable, order)
    left_group = exponential_polynomial(sp.I * left_generator, variable, order)
    right_group_adjoint = exponential_polynomial(-sp.I * right_generator, variable, order)
    left_group_adjoint = exponential_polynomial(-sp.I * left_generator, variable, order)
    transformed_mass = right_group * mass_source * left_group_adjoint
    # The current sources are r_mu=(partial_mu a/f_a)c_R(a/f_a),
    # l_mu=(partial_mu a/f_a)c_L(a/f_a); take out that common factor.
    transformed_right_current = (right_group * right_current * right_group_adjoint
                                + sp.I * right_group * sp.diff(right_group_adjoint, variable))
    transformed_left_current = (left_group * left_current * left_group_adjoint
                               + sp.I * left_group * sp.diff(left_group_adjoint, variable))
    return {
        "variable": variable,
        "mass_source": mass_source,
        "right_current": right_current,
        "left_current": left_current,
        "mass_derivatives": derivative_coefficients(transformed_mass, variable, order),
        # Differentiating g^dagger consumes one polynomial degree.
        "right_current_derivatives": derivative_coefficients(transformed_right_current, variable, order - 1),
        "left_current_derivatives": derivative_coefficients(transformed_left_current, variable, order - 1),
    }


def main():
    result = derive()
    print("Noncommuting flavour matrices; g_R=exp[i(a/f_a)(delta+kappa)], "
          "g_L=exp[i(a/f_a)(delta-kappa)].")
    print("chi_n denotes d^n chi/d(a/f_a)^n at zero; c_R_n and c_L_n "
          "denote analogous derivatives of derivative-current coefficients.")
    print("Each current expression multiplies partial_mu a/f_a. "
          "The displayed derivative coefficients carry no Taylor factorial.")
    show("Input mass spurion through fourth degree", result["mass_source"])
    show("Input right derivative-current coefficient", result["right_current"])
    show("Input left derivative-current coefficient", result["left_current"])
    for name, label in (("mass_derivatives", "Transformed mass spurion"),
                        ("right_current_derivatives", "Transformed right derivative current"),
                        ("left_current_derivatives", "Transformed left derivative current")):
        for degree, expression in enumerate(result[name]):
            show(f"{label}, derivative {degree} at zero", expression)


if __name__ == "__main__":
    main()
