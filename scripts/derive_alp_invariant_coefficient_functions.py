"""Derive partonic ALP coupling invariants from coefficient functions.

Run: .venv/bin/python scripts/derive_alp_invariant_coefficient_functions.py --order 2

Natural units; x=a/f_a is dimensionless. C_R(x), C_L(x) are Hermitian
flavor matrices multiplying the right/left derivative currents. Their
transformation is C'=g C g^{-1}+i g d(g^{-1})/dx. M(x) multiplies
bar(q_R) M q_L and transforms as g_R M g_L^{-1}. All flavor products
are noncommuting. The field redefinitions obey g_R(0)=g_L(0)=1.

W_eta solves dW_eta/dx=i C_eta W_eta with W_eta(0)=1. Hermiticity of
C_eta makes W_eta unitary. Hence W_R^{-1} M W_L equals W_R^dagger M W_L.
Taylor coefficients printed below are derivatives at x=0, with the Taylor
factorial removed by differentiation. No expansion in flavor matrices is
made; the displayed local Taylor series stops at the requested order in x.
The general constant chiral redefinition uses G_R=delta+kappa and
G_L=delta-kappa. The current coefficients include the inhomogeneous kinetic
connection, while scalar/pseudoscalar coefficients transform homogeneously.
"""

import argparse

import sympy as sp


def show(label, expression):
    """Print deterministic standard LaTeX."""
    print(f"{label}:\n{sp.latex(expression)}\n")


def connection_transport_covariance():
    """Use the derivative-of-inverse identity in the transformed ODE."""
    group, group_derivative, connection, transport = sp.symbols(
        "g dg C W", commutative=False)
    inverse_derivative = -group**-1 * group_derivative * group**-1
    transformed_connection = (group * connection * group**-1
                              + sp.I * group * inverse_derivative)
    derivative_of_transformed_transport = (group_derivative * transport
                                           + group * sp.I * connection * transport)
    transformed_ode_rhs = sp.I * transformed_connection * group * transport
    return (transformed_connection,
            sp.expand(derivative_of_transformed_transport - transformed_ode_rhs))


def transport_series(connection_coefficients, variable, order):
    """Integrate dW/dx=i C W by coefficient extraction, preserving order.

    SymPy's ordinary polynomial operations suffice: the coefficient of
    x^n in dW/dx is (n+1) times the next coefficient of W.
    """
    connection = sp.Add(*(coefficient * variable**degree / sp.factorial(degree)
                          for degree, coefficient in enumerate(connection_coefficients)))
    transport = sp.Integer(1)
    for degree in range(order):
        next_coefficient = sp.expand(sp.I * connection * transport).coeff(variable, degree)
        transport += next_coefficient * variable**(degree + 1) / (degree + 1)
    return sp.expand(transport)


def invariant_mass_coefficients(variable, order=2):
    right_coefficients = sp.symbols(f"C_R0:{order}", commutative=False, hermitian=True)
    left_coefficients = sp.symbols(f"C_L0:{order}", commutative=False, hermitian=True)
    mass_coefficients = sp.symbols(f"M_0:{order + 1}", commutative=False)
    right_transport = transport_series(right_coefficients, variable, order)
    left_transport = transport_series(left_coefficients, variable, order)
    mass = sp.Add(*(coefficient * variable**degree / sp.factorial(degree)
                   for degree, coefficient in enumerate(mass_coefficients)))
    invariant_mass = sp.adjoint(right_transport) * mass * left_transport
    derivatives = tuple(sp.expand(sp.diff(invariant_mass, variable, degree).subs(variable, 0))
                        for degree in range(order + 1))
    return right_transport, left_transport, derivatives


def axial_second_order_shift(variable):
    axial_generator, vector_zero, axial_zero, vector_one, axial_one = sp.symbols(
        r"\kappa c_V0 c_A0 c_V1 c_A1", commutative=False)
    positive_group = sp.Add(*(sp.I**power * axial_generator**power * variable**power
                             / sp.factorial(power) for power in range(3)))
    negative_group = sp.Add(*((-sp.I)**power * axial_generator**power * variable**power
                             / sp.factorial(power) for power in range(3)))
    right_current = vector_zero + axial_zero + variable * (vector_one + axial_one)
    left_current = vector_zero - axial_zero + variable * (vector_one - axial_one)
    transformed_right = (positive_group * right_current * negative_group
                         + sp.I * positive_group * sp.diff(negative_group, variable))
    transformed_left = (negative_group * left_current * positive_group
                        + sp.I * negative_group * sp.diff(positive_group, variable))
    transformed_vector = (transformed_right + transformed_left) / 2
    vector_one_transformed = sp.expand(sp.diff(transformed_vector, variable).subs(variable, 0))
    second_vector_coefficient = sp.Symbol("c_q^{V,2}", commutative=False)
    user_normalization = sp.expand(vector_one_transformed.subs(
        vector_one, sp.I * second_vector_coefficient) / sp.I)
    return vector_one_transformed, user_normalization


def exponential_series(generator, variable, order):
    """Standard exponential polynomial for a noncommuting flavor generator."""
    return sp.Add(*(generator**degree * variable**degree / sp.factorial(degree)
                    for degree in range(order + 1)))


def taylor_series(coefficients, variable):
    """Build a Taylor polynomial whose entries are derivatives at zero."""
    return sp.Add(*(coefficient * variable**degree / sp.factorial(degree)
                    for degree, coefficient in enumerate(coefficients)))


def derivatives_at_zero(expression, variable, order):
    return tuple(sp.expand(sp.diff(expression, variable, degree).subs(variable, 0))
                 for degree in range(order + 1))


def general_chiral_coefficients(variable, order):
    """Differentiate the actual source products and inhomogeneous connections.

    One extra exponential degree is kept because differentiating the inverse
    group matrix consumes a Taylor degree. G_R and G_L are left as the
    physical generator matrices until the explicit delta/kappa substitution.
    """
    right_generator, left_generator = sp.symbols(
        "G_R G_L", commutative=False, hermitian=True)
    vector_generator, axial_generator = sp.symbols(
        r"\delta \kappa", commutative=False, hermitian=True)
    vector_coefficients = sp.symbols(f"c_V0:{order + 1}", commutative=False, hermitian=True)
    axial_coefficients = sp.symbols(f"c_A0:{order + 1}", commutative=False, hermitian=True)
    scalar_coefficients = sp.symbols(f"S_0:{order + 1}", commutative=False, hermitian=True)
    pseudoscalar_coefficients = sp.symbols(f"P_0:{order + 1}", commutative=False, hermitian=True)
    vector_source = taylor_series(vector_coefficients, variable)
    axial_source = taylor_series(axial_coefficients, variable)
    scalar_source = taylor_series(scalar_coefficients, variable)
    pseudoscalar_source = taylor_series(pseudoscalar_coefficients, variable)
    right_source = vector_source + axial_source
    left_source = vector_source - axial_source
    mass_source = scalar_source + sp.I * pseudoscalar_source
    right_group = exponential_series(sp.I * right_generator, variable, order + 1)
    right_inverse = exponential_series(-sp.I * right_generator, variable, order + 1)
    left_group = exponential_series(sp.I * left_generator, variable, order + 1)
    left_inverse = exponential_series(-sp.I * left_generator, variable, order + 1)
    transformed_right = (right_group * right_source * right_inverse
                         + sp.I * right_group * sp.diff(right_inverse, variable))
    transformed_left = (left_group * left_source * left_inverse
                        + sp.I * left_group * sp.diff(left_inverse, variable))
    transformed_vector = (transformed_right + transformed_left) / 2
    transformed_axial = (transformed_right - transformed_left) / 2
    transformed_mass = right_group * mass_source * left_inverse
    transformed_mass_adjoint = sp.adjoint(transformed_mass)
    transformed_scalar = (transformed_mass + transformed_mass_adjoint) / 2
    transformed_pseudoscalar = (transformed_mass - transformed_mass_adjoint) / (2 * sp.I)
    sources = {
        "C_R": transformed_right,
        "C_L": transformed_left,
        "c_V": transformed_vector,
        "c_A": transformed_axial,
        "M=S+iP": transformed_mass,
        "S": transformed_scalar,
        "P": transformed_pseudoscalar,
    }
    coefficients = {name: derivatives_at_zero(source, variable, order)
                    for name, source in sources.items()}
    generator_substitution = {
        right_generator: vector_generator + axial_generator,
        left_generator: vector_generator - axial_generator,
    }
    return coefficients, generator_substitution


def general_leibniz_laws():
    """Represent the general derivative product rule with standard Sum objects.

    The finite-order expressions below are calculated by differentiating the
    exponential/source products themselves. These sums display the same
    product rule for a symbolic nonnegative derivative number n.
    """
    derivative_number, left_power, right_power = sp.symbols("n r s", integer=True, nonnegative=True)
    right_generator, left_generator = sp.symbols("G_R G_L", commutative=False)
    right_coefficient = sp.Function("C_R", commutative=False)
    left_coefficient = sp.Function("C_L", commutative=False)
    scalar_coefficient = sp.Function("S", commutative=False)
    pseudoscalar_coefficient = sp.Function("P", commutative=False)
    source_degree = derivative_number - left_power - right_power
    multinomial = (sp.factorial(derivative_number) / sp.factorial(left_power)
                   / sp.factorial(right_power) / sp.factorial(source_degree))

    def product_derivative(left_generator_value, source_value, right_generator_value):
        summand = (multinomial * (sp.I * left_generator_value)**left_power
                   * source_value * (-sp.I * right_generator_value)**right_power)
        return sp.Sum(summand, (right_power, 0, derivative_number - left_power),
                      (left_power, 0, derivative_number))

    right_homogeneous = product_derivative(right_generator, right_coefficient(source_degree),
                                           right_generator)
    left_homogeneous = product_derivative(left_generator, left_coefficient(source_degree),
                                          left_generator)
    right_full = right_homogeneous + sp.KroneckerDelta(derivative_number, 0) * right_generator
    left_full = left_homogeneous + sp.KroneckerDelta(derivative_number, 0) * left_generator
    mass = product_derivative(right_generator,
                              scalar_coefficient(source_degree) + sp.I * pseudoscalar_coefficient(source_degree),
                              left_generator)
    mass_adjoint = product_derivative(left_generator,
                                      scalar_coefficient(source_degree) - sp.I * pseudoscalar_coefficient(source_degree),
                                      right_generator)
    return {"C_R": right_full, "C_L": left_full,
            "c_V": (right_full + left_full) / 2,
            "c_A": (right_full - left_full) / 2,
            "M=S+iP": mass, "S": (mass + mass_adjoint) / 2,
            "P": (mass - mass_adjoint) / (2 * sp.I)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--order", type=int, default=2,
                        help="largest invariant Taylor derivative to print (default: 2)")
    arguments = parser.parse_args()
    if arguments.order < 1:
        parser.error("--order must be positive")
    variable = sp.Symbol("x", real=True)
    print("Partonic coefficient functions; x=a/f_a; natural units.")
    print("C_R,C_L are Hermitian connections; M is a general chirality-changing source.")
    print("g_eta(0)=1; all products retain noncommuting flavor order.")
    transformed_connection, transport_residual = connection_transport_covariance()
    show("Connection after substituting d(g^{-1})/dx=-g^{-1}(dg/dx)g^{-1}",
         sp.expand(transformed_connection))
    show("Product-rule residual proving W'=g W solves the transformed transport ODE",
         transport_residual)
    right_group, left_group, right_transport, left_transport, mass = sp.symbols(
        "g_R g_L W_R W_L M", commutative=False)
    transformed_invariant = ((right_group * right_transport)**-1
                             * (right_group * mass * left_group**-1)
                             * left_group * left_transport)
    show("Invariant W_R^{-1} M W_L after the field redefinition",
         sp.expand(transformed_invariant))
    right_series, left_series, invariant_derivatives = invariant_mass_coefficients(
        variable, arguments.order)
    show(f"Right transport W_R through order {arguments.order}", right_series)
    show(f"Left transport W_L through order {arguments.order}", left_series)
    for degree, derivative in enumerate(invariant_derivatives):
        show(f"Derivative {degree} at zero of W_R^dagger M W_L", derivative)
    vector_shift, user_normalization = axial_second_order_shift(variable)
    show("Pure axial rotation: transformed derivative c_V1", vector_shift)
    show("Same law when c_V1=i c_q^{V,2}", user_normalization)
    print("With the explicit prefactor i in the latter convention, c_q^{V,2}")
    print("is anti-Hermitian when the Lagrangian is Hermitian; i c_q^{V,2} is Hermitian.")
    print("General constant chiral redefinition: G_R=delta+kappa, G_L=delta-kappa.")
    print("All following coefficients are derivatives at zero and include no Taylor factorial.")
    for source_name, source_law in general_leibniz_laws().items():
        show(f"General nth-derivative Leibniz law for {source_name}", source_law)
    general_coefficients, generator_substitution = general_chiral_coefficients(
        variable, arguments.order)
    for source_name, coefficients in general_coefficients.items():
        for degree, coefficient in enumerate(coefficients):
            show(f"General chiral transformation of {source_name}, derivative {degree}", coefficient)
    for source_name in ("c_V", "c_A", "S", "P"):
        for degree in range(min(arguments.order, 1) + 1):
            show(f"Explicit delta/kappa law for {source_name}, derivative {degree}",
                 sp.expand(general_coefficients[source_name][degree].subs(generator_substitution)))


if __name__ == "__main__":
    main()
