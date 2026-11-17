"""Does the double-trace source term enter the vector masses and the KK pi current?

Operator (declared in memory/rcht/foundations/rcht-photon-vector-diagonalization.md):

    L_y = y_V / (2 sqrt 2) <V_{mu nu}> <chi_+ f_+^{mu nu}>,

which, for f_+ = -2 e Q F in the Goldstone vacuum, adds y_V <chi_+ Q> 1 to the
photon source matrix S_gamma = F_V Q + z_V {chi_+, Q} + y_V <chi_+ Q> 1.

Conventions (Alda:2026nrd): Phi = sum_P lambda_P P with tr(lambda_i lambda_j) = 2 delta_ij,
including the singlet direction (eta' component), u = exp(i Phi / (2 f)),
chi = 2 B_0 diag(m_u, m_d, m_s) real, chi_+ = u^dag chi u^dag + u chi^dag u,
f_+ = u F_L u^dag + u^dag F_R u with F_L = F_R = -e Q F, Q = diag(2/3, -1/3, -1/3).
The vector nonet V = sum_i T_i V^i uses <T_i T_j> = delta_ij. The single Lorentz
structure V_{mu nu} F^{mu nu} is represented by commuting symbols, which is
sufficient to count fields and flavor traces. Units: x_q = (chi_+)_qq at the
vacuum in GeV^2, F_V in GeV, z_V and y_V in GeV^-1, e dimensionless.

Results checked by assertions:
  (a) L_y is linear in V and proportional to e F: it gives no quadratic vector mass
      term. Through photon-vector mixing it changes only the O(e^2) shifts
      delta m_i^2 = 2 e^2 R_i^2, which the fit omits for every source term.
  (b) Only <V> appears, so rho^0, rho^+-, K^* never couple; the Goldstone
      expansion of <chi_+ f_+> has no odd powers of Phi; a V->P bilinear vanishes.
      Hence the tree-level KK pi current receives y_V only through the photon
      residues of the ideal omega_I and phi_I states:
      F_l -> F_l + 2 y_V (x_l - x_s) in the omega_I source,
      F_s -> F_s - y_V (x_l - x_s) in the phi_I source, rho unchanged.

Run: .venv/bin/python scripts/derive_double_trace_vector_mass_and_kkpi.py
Tests: .venv/bin/python -m unittest scripts.test_derive_double_trace_vector_mass_and_kkpi
"""
import sympy as sp

EXPANSION_ORDER = 3

e, field_strength, decay_constant = sp.symbols('e F f', positive=True)
y_V, z_V, F_V = sp.symbols('y_V z_V F_V', real=True)
chi_u, chi_d, chi_s = sp.symbols('chi_u chi_d chi_s', real=True)
x_l, x_s = sp.symbols('x_l x_s', real=True)
CHARGE_MATRIX = sp.diag(sp.Rational(2, 3), -sp.Rational(1, 3), -sp.Rational(1, 3))


def hermitian_matrix(prefix):
    """General Hermitian 3x3 matrix with nine real components named by prefix."""
    diagonal = sp.symbols(f'{prefix}_11 {prefix}_22 {prefix}_33', real=True)
    real_parts = sp.symbols(f'{prefix}_12_re {prefix}_13_re {prefix}_23_re', real=True)
    imaginary_parts = sp.symbols(f'{prefix}_12_im {prefix}_13_im {prefix}_23_im', real=True)
    matrix = sp.diag(*diagonal)
    for (row, column), real, imaginary in zip(((0, 1), (0, 2), (1, 2)), real_parts, imaginary_parts):
        matrix[row, column] = real + sp.I*imaginary
        matrix[column, row] = real - sp.I*imaginary
    return matrix, diagonal + real_parts + imaginary_parts


def goldstone_matrix():
    """Phi with nine real fields: octet plus singlet (eta' direction)."""
    return hermitian_matrix('phi')


def vector_nonet():
    """Vector nonet matrix with Lorentz indices suppressed."""
    return hermitian_matrix('V')


def exponential_terms(matrix, order=EXPANSION_ORDER):
    """Terms matrix^n / n! of exp(matrix), for n = 0, ..., order."""
    terms = [sp.eye(3)]
    for n in range(1, order + 1):
        terms.append((terms[-1]*matrix/n).applyfunc(sp.expand))
    return terms


def series_product(*factors, order=EXPANSION_ORDER):
    """Product of matrix power series given as lists of order-n coefficients."""
    result = factors[0]
    for factor in factors[1:]:
        result = [sum((result[m]*factor[n - m] for m in range(n + 1)), sp.zeros(3)).applyfunc(sp.expand)
                  for n in range(order + 1)]
    return result


def constant_series(matrix, order=EXPANSION_ORDER):
    """A field-independent matrix written as a power series in Phi."""
    return [matrix] + [sp.zeros(3)]*order


def source_trace_expansion():
    """Coefficients of <chi_+ f_+> in powers of Phi, up to EXPANSION_ORDER."""
    phi, fields = goldstone_matrix()
    argument = sp.I*phi/(2*decay_constant)
    u = exponential_terms(argument)
    u_dagger = exponential_terms(-argument)
    chi = constant_series(sp.diag(chi_u, chi_d, chi_s))
    external_field_strength = constant_series(-e*CHARGE_MATRIX*field_strength)
    chi_plus = [a + b for a, b in zip(series_product(u_dagger, chi, u_dagger), series_product(u, chi, u))]
    f_plus = [a + b for a, b in zip(series_product(u, external_field_strength, u_dagger),
                                    series_product(u_dagger, external_field_strength, u))]
    coefficients = [sp.expand(term.trace()) for term in series_product(chi_plus, f_plus)]
    return dict(coefficients=coefficients, fields=fields, chi_plus=chi_plus, f_plus=f_plus)


def double_trace_lagrangian(source_trace):
    """L_y = y_V/(2 sqrt 2) <V> <chi_+ f_+> for a general vector nonet."""
    vector, vector_fields = vector_nonet()
    return y_V/(2*sp.sqrt(2))*vector.trace()*source_trace, vector, vector_fields


def ideal_generators():
    """Normalized ideal generators of rho^0, omega_I, phi_I (positive strange phase)."""
    return dict(rho=sp.diag(1, -1, 0)/sp.sqrt(2), omega_I=sp.diag(1, 1, 0)/sp.sqrt(2),
                phi_I=sp.diag(0, 0, 1))


def photon_residues(include_double_trace=True):
    """R_i = <T_i S_gamma> in the isospin limit chi_+ = diag(x_l, x_l, x_s)."""
    chi_plus = sp.diag(x_l, x_l, x_s)
    source = F_V*CHARGE_MATRIX + z_V*(chi_plus*CHARGE_MATRIX + CHARGE_MATRIX*chi_plus)
    if include_double_trace:
        source += y_V*(chi_plus*CHARGE_MATRIX).trace()*sp.eye(3)
    return {name: sp.simplify((generator*source).trace())
            for name, generator in ideal_generators().items()}


def kkpi_source_factors():
    """Photon insertions on the neutral vector lines used by the BaBar current.

    With the ideal-basis neutral propagator matrix (Pi_ll, Pi_ls, Pi_ss) and the
    rho propagator Pi_rho, the photon creates rho^0, omega_I, phi_I with weights
    Pi R. The fit code writes them as source_rho/sqrt 2, source_light/(3 sqrt 2),
    -source_strange/3 (scripts/fit_babar_kkpi_global_vector_sector.py).
    """
    pi_rho, pi_ll, pi_ls, pi_ss = sp.symbols('Pi_rho Pi_ll Pi_ls Pi_ss')
    F_l, F_s = sp.symbols('F_l F_s', real=True)
    residues = photon_residues()
    weights = dict(rho=pi_rho*residues['rho'],
                   omega_I=pi_ll*residues['omega_I'] + pi_ls*residues['phi_I'],
                   phi_I=pi_ls*residues['omega_I'] + pi_ss*residues['phi_I'])

    def code_form(light, strange):
        return dict(source_rho=light*pi_rho,
                    source_light=light*pi_ll - sp.sqrt(2)*strange*pi_ls,
                    source_strange=strange*pi_ss - light*pi_ls/sp.sqrt(2))
    single_light = F_V + 2*z_V*x_l
    single_strange = F_V + 2*z_V*x_s
    shifted_light = single_light + 2*y_V*(x_l - x_s)
    shifted_strange = single_strange - y_V*(x_l - x_s)
    rho_form = code_form(single_light, single_strange)
    isoscalar_form = code_form(shifted_light, shifted_strange)
    code_weights = dict(rho=rho_form['source_rho']/sp.sqrt(2),
                        omega_I=isoscalar_form['source_light']/(3*sp.sqrt(2)),
                        phi_I=-isoscalar_form['source_strange']/3)
    differences = {name: sp.simplify(weights[name] - code_weights[name]) for name in weights}
    return dict(weights=weights, code_weights=code_weights, differences=differences,
                shifted_light=shifted_light, shifted_strange=shifted_strange)


def leading_electromagnetic_mass_shift():
    """O(e^2) pole shift from tensor photon-vector mixing, two-vector check.

    Uses det(s K_V - D) = 0 with K_V = 1 - eps eps^T / (1 + eps^T eps),
    eps_i = sqrt 2 e R_i / M_i (rcht-photon-vector-diagonalization.md).
    """
    mass_1, mass_2 = sp.symbols('M_1 M_2', positive=True)
    residue_1, residue_2 = sp.symbols('R_1 R_2', real=True)
    s = sp.Symbol('s')
    epsilon_vector = sp.Matrix([sp.sqrt(2)*e*residue_1/mass_1, sp.sqrt(2)*e*residue_2/mass_2])
    normalization = 1 + (epsilon_vector.T*epsilon_vector)[0]
    kinetic = sp.eye(2) - epsilon_vector*epsilon_vector.T/normalization
    determinant = sp.together((s*kinetic - sp.diag(mass_1**2, mass_2**2)).det())
    numerator = sp.numer(determinant)
    # Pole near M_1^2: s = M_1^2 + shift, shift = c e^2 + O(e^4).
    c = sp.Symbol('c')
    series = sp.series(numerator.subs(s, mass_1**2 + c*e**2), e, 0, 4).removeO()
    leading_coefficient = sp.solve(sp.expand(series).coeff(e, 2), c)
    return dict(leading_shift=[sp.simplify(value*e**2) for value in leading_coefficient],
                residue=residue_1)


def vector_to_pseudoscalar_bilinear():
    """Contraction of an antisymmetric tensor with p^mu p^nu, the only V-P structure."""
    momentum = sp.symbols('p0:4')
    tensor = sp.Matrix(4, 4, lambda mu, nu: sp.Symbol(f'V_{min(mu, nu)}{max(mu, nu)}')
                       * (1 if mu < nu else (-1 if mu > nu else 0)))
    return sp.expand(sum(tensor[mu, nu]*momentum[mu]*momentum[nu] for mu in range(4) for nu in range(4)))


def derive():
    """Run every step and return the intermediate expressions."""
    expansion = source_trace_expansion()
    coefficients = expansion['coefficients']
    lagrangian, vector, vector_fields = double_trace_lagrangian(sp.Add(*coefficients))
    lagrangian = sp.expand(lagrangian)
    vector_degree = sp.Poly(lagrangian, *vector_fields).total_degree()
    hessian = sp.hessian(lagrangian, vector_fields)
    off_diagonal_vector_coefficients = [sp.diff(lagrangian, field) for field in vector_fields[3:]]
    traceless_derivatives = dict(
        rho0=sp.expand(sp.diff(lagrangian, vector_fields[0]) - sp.diff(lagrangian, vector_fields[1])),
        omega8=sp.expand(sp.diff(lagrangian, vector_fields[0]) + sp.diff(lagrangian, vector_fields[1])
                           - 2*sp.diff(lagrangian, vector_fields[2])))
    vacuum_check = sp.simplify(coefficients[0] + 4*e*field_strength*(sp.diag(chi_u, chi_d, chi_s)*CHARGE_MATRIX).trace())
    residues_with = photon_residues(True)
    residues_without = photon_residues(False)
    residue_shift = {name: sp.simplify(residues_with[name] - residues_without[name]) for name in residues_with}
    return dict(
        expansion_coefficients=coefficients,
        lagrangian_vanishes_without_photon=sp.expand(lagrangian.subs(e, 0)) == 0,
        vector_degree=vector_degree,
        hessian_is_zero=hessian.is_zero_matrix,
        off_diagonal_vector_coefficients=off_diagonal_vector_coefficients,
        traceless_derivatives=traceless_derivatives,
        vacuum_check=vacuum_check,
        residues_with=residues_with,
        residue_shift=residue_shift,
        kkpi=kkpi_source_factors(),
        electromagnetic_mass_shift=leading_electromagnetic_mass_shift(),
        vector_to_pseudoscalar=vector_to_pseudoscalar_bilinear(),
        weak_current_trace=(sp.diag(x_l, x_l, x_s)*sp.Matrix(3, 3, lambda i, j: 1 if (i, j) == (0, 1) else 0)).trace(),
    )


def main():
    result = derive()
    coefficients = result['expansion_coefficients']
    assert result['vacuum_check'] == 0
    assert coefficients[1] == 0 and coefficients[3] == 0
    assert coefficients[2] != 0
    assert result['lagrangian_vanishes_without_photon']
    assert result['vector_degree'] == 1 and result['hessian_is_zero']
    assert all(sp.expand(c) == 0 for c in result['off_diagonal_vector_coefficients'])
    assert all(value == 0 for value in result['traceless_derivatives'].values())
    assert result['residue_shift']['rho'] == 0
    assert sp.simplify(result['residue_shift']['omega_I'] - sp.sqrt(2)*y_V*(x_l - x_s)/3) == 0
    assert sp.simplify(result['residue_shift']['phi_I'] - y_V*(x_l - x_s)/3) == 0
    assert all(value == 0 for value in result['kkpi']['differences'].values())
    residue = result['electromagnetic_mass_shift']['residue']
    assert result['electromagnetic_mass_shift']['leading_shift'] == [2*e**2*residue**2]
    assert result['vector_to_pseudoscalar'] == 0
    assert result['weak_current_trace'] == 0

    print('<chi_+ f_+> at Phi^0:', sp.factor(coefficients[0]))
    print('<chi_+ f_+> at Phi^1 and Phi^3:', coefficients[1], coefficients[3])
    print('<chi_+ f_+> at Phi^2 is nonzero; number of terms:', len(sp.Add.make_args(coefficients[2])))
    print('L_y vanishes at e = 0:', result['lagrangian_vanishes_without_photon'])
    print('degree of L_y in V:', result['vector_degree'], '; Hessian in V is zero:', result['hessian_is_zero'])
    print('off-diagonal (rho^+-, K^*) and traceless neutral (rho^0, omega_8) couplings vanish: True')
    print('photon residues with y_V:', {k: sp.latex(v) for k, v in result['residues_with'].items()})
    print('residue shifts from y_V:', {k: sp.latex(v) for k, v in result['residue_shift'].items()})
    print('KK pi current: omega_I source uses F_l ->', sp.latex(result['kkpi']['shifted_light']),
          '; phi_I source uses F_s ->', sp.latex(result['kkpi']['shifted_strange']), '; rho unchanged')
    print('leading O(e^2) photon-mixing mass shift: delta m^2 =',
          sp.latex(result['electromagnetic_mass_shift']['leading_shift'][0]))
    print('antisymmetric V contracted with p p:', result['vector_to_pseudoscalar'])
    print('<chi_+ T> for the charged weak current generator:', result['weak_current_trace'])
    print('All assertions passed.')


if __name__ == '__main__':
    main()
