"""Derive source transformations and check their composition.

Run: .venv/bin/python scripts/derive_alp_chiral_source_transformations.py
Natural units; x=a/f_a is dimensionless. Matrices are U(2) in the series
checks, and all products keep their order. Results are checked through x^3.
The anomaly law and extension of X_TC are inputs from the linked memory record.
"""
import sympy as sp

x = sp.Symbol('x', real=True)
I = sp.I
ORDER = 4  # One extra degree is needed when differentiating.


def truncate(expression, order=ORDER):
    """Use SymPy's polynomial coefficients for an exact Taylor truncation."""
    return sp.Add(*(sp.expand(expression).coeff(x, degree) * x**degree
                    for degree in range(order + 1)))


def matrix_truncate(matrix, order=ORDER):
    return matrix.applyfunc(lambda entry: truncate(entry, order))


def multiply(*matrices, order=ORDER):
    result = matrices[0]
    for matrix in matrices[1:]:
        result = matrix_truncate(result * matrix, order)
    return result


def exponential(generator):
    """Taylor-expand an exponential with generator vanishing at x=0."""
    result = sp.eye(2)
    power = sp.eye(2)
    for degree in range(1, ORDER + 1):
        power = multiply(power, generator)
        result += power / sp.factorial(degree)
    return matrix_truncate(result)


def transport(connection):
    """Solve dW/dx=i C W recursively by matching Taylor coefficients."""
    result = sp.eye(2)
    for degree in range(ORDER):
        rhs = multiply(I * connection, result)
        coefficient = rhs.applyfunc(lambda entry: sp.expand(entry).coeff(x, degree))
        result += coefficient * x**(degree + 1) / (degree + 1)
    return result


def transform_connection(group, connection):
    return matrix_truncate(multiply(group, connection, group.H)
                           + I * multiply(group, sp.diff(group.H, x)), ORDER - 1)


def check_series():
    """Exact symbolic x-series for independent noncommuting U(2) generators."""
    sigma1 = sp.Matrix([[0, 1], [1, 0]])
    sigma2 = sp.Matrix([[0, -I], [I, 0]])
    sigma3 = sp.diag(1, -1)
    identity = sp.eye(2)
    right = sigma1 + x * (sigma2 + identity) / 3 + x**2 * sigma3 / 5
    left = sigma2 + x * (sigma3 - identity) / 2 + x**2 * sigma1 / 7
    right_first = exponential(I * (x * sigma3 + x**2 * sigma1 / 2))
    left_first = exponential(I * (x * sigma1 + x**2 * sigma2 / 3))
    right_second = exponential(I * (x * (sigma2 + identity) / 2 + x**3 * sigma1 / 3))
    left_second = exponential(I * (x * (sigma3 - identity) / 3 + x**2 * sigma1 / 4))
    goldstone = exponential(I * (x * sigma1 + x**2 * sigma3 / 2))
    mass = identity + x * sigma3 + I * x**2 * sigma2
    right_transport, left_transport = transport(right), transport(left)
    dressed_goldstone = multiply(right_transport.H, goldstone, left_transport)
    dressed_mass = multiply(right_transport.H, mass, left_transport)
    covariant_derivative = (sp.diff(goldstone, x) - I * multiply(right, goldstone)
                            + I * multiply(goldstone, left))
    derivative_residual = matrix_truncate(sp.diff(dressed_goldstone, x)
        - multiply(right_transport.H, covariant_derivative, left_transport), ORDER - 1)
    transformed_right = transform_connection(right_first, right)
    transformed_left = transform_connection(left_first, left)
    final_right = transform_connection(right_second, transformed_right)
    final_left = transform_connection(left_second, transformed_left)
    combined_right = multiply(right_second, right_first)
    combined_left = multiply(left_second, left_first)
    direct_right = transform_connection(combined_right, right)
    direct_left = transform_connection(combined_left, left)
    final_goldstone = multiply(combined_right, goldstone, combined_left.H)
    final_mass = multiply(combined_right, mass, combined_left.H)
    final_right_transport, final_left_transport = transport(final_right), transport(final_left)
    results = {
        'right connection composition': final_right - direct_right,
        'left connection composition': final_left - direct_left,
        'right transport covariance': final_right_transport - multiply(combined_right, right_transport),
        'left transport covariance': final_left_transport - multiply(combined_left, left_transport),
        'dressed U invariance': multiply(final_right_transport.H, final_goldstone, final_left_transport) - dressed_goldstone,
        'dressed mass invariance': multiply(final_right_transport.H, final_mass, final_left_transport) - dressed_mass,
        'covariant derivative': derivative_residual,
    }
    return {name: matrix_truncate(residual, 3) for name, residual in results.items()}


def derive_anomaly():
    """Use the stated Jacobian and the determinant of a singlet rotation."""
    flavor_number = sp.Symbol('N_f', integer=True, positive=True)
    theta = sp.Function('theta')(x)
    Theta = sp.Function('Theta')(x)
    trace_log_u = sp.Symbol(r'\operatorname{tr}\log u')
    alpha_s = sp.Symbol(r'\alpha_s', positive=True)
    Theta_raw = sp.Symbol(r'\Theta_{raw}', real=True)
    # det exp(+-i theta) 1 = exp(+-i N_f theta).
    determinant_ratio = sp.exp(2 * I * flavor_number * theta)
    logarithmic_derivative = sp.simplify(sp.diff(determinant_ratio, x) / determinant_ratio)
    # Integrating the logarithmic derivative fixes the branch at theta(0)=0.
    log_determinant_ratio = sp.integrate(logarithmic_derivative, x)
    jacobian_shift = sp.simplify(-log_determinant_ratio / (2 * I))
    trace_log_u_shift = log_determinant_ratio / 2
    original = 2 * I * (trace_log_u + I * Theta)
    transformed = 2 * I * (trace_log_u + trace_log_u_shift + I * (Theta + jacobian_shift))
    normalized = 4 * sp.pi * Theta_raw / alpha_s
    raw_shift = sp.simplify(alpha_s * jacobian_shift / (4 * sp.pi))
    return {'singlet determinant ratio': determinant_ratio,
            'derivative of its logarithm': logarithmic_derivative,
            'Jacobian shift of Theta': jacobian_shift,
            'shift of tr log u': trace_log_u_shift,
            'X_TC residual': sp.simplify(transformed - original),
            'raw shift of Theta': raw_shift,
            'normalization residual': sp.simplify(4 * sp.pi * (Theta_raw + raw_shift) / alpha_s - normalized - jacobian_shift)}


def derive_source_laws():
    """Apply the partonic connection law and the transport product rule."""
    group_r, group_l, derivative_group_r, derivative_group_l = sp.symbols(
        'g_R g_L dg_R dg_L', commutative=False)
    current_r, current_l, spurion, goldstone = sp.symbols(
        'C_R C_L chi U', commutative=False)
    transport_r, transport_l, derivative_goldstone = sp.symbols(
        'W_R W_L dU', commutative=False)
    source_r, source_l = sp.symbols('r_mu ell_mu', commutative=False)
    derivative_x = sp.Symbol('d_mu_x', real=True)
    # d(g^-1) = -g^-1 (dg) g^-1, and r_mu = (d_mu x) C_R.
    transformed_r = sp.expand(derivative_x * (
        group_r * current_r * group_r**-1 - I * derivative_group_r * group_r**-1))
    transformed_l = sp.expand(derivative_x * (
        group_l * current_l * group_l**-1 - I * derivative_group_l * group_l**-1))
    transformed_spurion = group_r * spurion * group_l**-1
    dressed_derivative = sp.expand(
        -I * transport_r**-1 * source_r * goldstone * transport_l
        + transport_r**-1 * derivative_goldstone * transport_l
        + I * transport_r**-1 * goldstone * source_l * transport_l)
    covariant_derivative = (derivative_goldstone - I * source_r * goldstone
                            + I * goldstone * source_l)
    residual = sp.expand(dressed_derivative
                         - transport_r**-1 * covariant_derivative * transport_l)
    return {'r_mu transformed': transformed_r, 'ell_mu transformed': transformed_l,
            'chi transformed': transformed_spurion,
            'derivative of dressed U': dressed_derivative,
            'product rule residual': residual}


def main():
    print('Natural units, x=a/f_a; U(2) series through x^3; continuous logarithm branch.')
    print('The Jacobian is an input; X_TC with general Theta is an assumed extension.')
    sources = derive_source_laws()
    for name, expression in sources.items():
        print(f'{name}:\n{sp.latex(expression)}\n')
    assert sources['product rule residual'] == 0
    results = check_series()
    for name, residual in results.items():
        print(f'{name} residual:\n{sp.latex(residual)}\n')
        assert residual == sp.zeros(2)
    anomaly = derive_anomaly()
    for name, expression in anomaly.items():
        print(f'{name}:\n{sp.latex(expression)}\n')
    assert anomaly['X_TC residual'] == 0
    assert anomaly['normalization residual'] == 0


if __name__ == '__main__':
    main()
