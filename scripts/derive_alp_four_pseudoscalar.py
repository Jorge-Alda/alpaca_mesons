"""Derive a -> 4P vector amplitudes and their invariant phase-space integral.

Run: .venv/bin/python scripts/derive_alp_four_pseudoscalar.py
Natural units (GeV), metric (+---), epsilon^{0123}=+1. Pseudoscalar
normalization is Phi=sum lambda_i P_i, Tr(lambda_i lambda_j)=2 delta_ij;
vector generators have unit trace normalization. The ALP sources are already
in the no-mixing basis. This calculation retains the original preliminary-fit
Lagrangian at first order in ALP sources and zero electromagnetic coupling.
It does not add independent VPPP operators. The generic tensor contractions
printed below specify topology and combinatorics; they are not numerical
channel amplitudes.
Vector exchanges use the same constant PDG total widths as a -> 3P:
1/(s-M_R^2+i M_R Gamma_R), with fitted neutral vector masses and the
charged-vector masses used in alpaca_chiral/three_pseudoscalar.py.
"""
from dataclasses import dataclass
from itertools import combinations, combinations_with_replacement, product
from functools import lru_cache

import sympy as sp

from derive_alp_vector_two_pseudoscalar import derive as derive_vvp
from derive_radiative_tff_amplitude import m, s1, s2


@dataclass(frozen=True)
class Vertex:
    """Field multiplicities after selecting exactly one ALP insertion."""

    name: str
    vectors: int
    pseudoscalars: int


# Both odd and even aVV vertices share field multiplicities. Their distinct
# flavour and Lorentz tensors must subsequently be retained separately.
SOURCE_VERTICES = (
    Vertex('aPP', 0, 2), Vertex('aPPP', 0, 3), Vertex('aPPPP', 0, 4),
    Vertex('aVP', 1, 1), Vertex('aVPP', 1, 2), Vertex('aVPPP', 1, 3),
    Vertex('aVV', 2, 0), Vertex('aVVP', 2, 1), Vertex('aVVPP', 2, 2),
)
STRONG_VERTICES = (
    Vertex('PPPP', 0, 4), Vertex('VPP', 1, 2),
    Vertex('VPPPP', 1, 4), Vertex('VVP', 2, 1),
    Vertex('VVPP', 2, 2), Vertex('VVPPP', 2, 3),
)


def enumerate_vector_trees():
    """Obtain every tree with one source, four external P and internal V.

    Every interaction has at least three total fields, including the ALP.
    For a tree with five external fields this bounds the number of vertices
    by three. A source aP two-point insertion is absent in the no-mixing
    basis. Quadratic masses are included in propagators. The small explicit
    graph enumeration is used because SymPy has no typed Feynman-tree
    enumerator; all analytical vertex polynomials below use SymPy.
    """
    results = {}
    for source in SOURCE_VERTICES:
        for strong_count in (1, 2):
            for strong in combinations_with_replacement(STRONG_VERTICES, strong_count):
                vertices = (source,) + strong
                vertex_count = len(vertices)
                possible_edges = tuple(combinations(range(vertex_count), 2))
                for edges in combinations(possible_edges, vertex_count - 1):
                    reached = {0}
                    for _ in range(vertex_count):
                        reached.update(right for left, right in edges if left in reached)
                        reached.update(left for left, right in edges if right in reached)
                    if len(reached) != vertex_count:
                        continue
                    for species in product(('P', 'V'), repeat=len(edges)):
                        if 'V' not in species:
                            continue
                        remaining_v = [vertex.vectors for vertex in vertices]
                        remaining_p = [vertex.pseudoscalars for vertex in vertices]
                        for (left, right), particle in zip(edges, species):
                            remaining = remaining_v if particle == 'V' else remaining_p
                            remaining[left] -= 1
                            remaining[right] -= 1
                        if min(remaining_v + remaining_p) < 0:
                            continue
                        if any(remaining_v) or sum(remaining_p) != 4:
                            continue
                        key = tuple(vertex.name for vertex in vertices)
                        # Equivalent internal labels are counted only once here.
                        # External assignments are summed separately, without a
                        # symmetry factor in the amplitude.
                        results.setdefault(key, (vertices, edges, species, tuple(remaining_p)))
    return tuple(results[key] for key in sorted(results))


def derive_odd_numerator():
    """Derive the two off-shell vector numerator from the operator contraction."""
    calculation = derive_vvp()
    source_mass_squared, first_virtuality, second_virtuality = sp.symbols(
        'm_X_squared s_R s_S', real=True)
    numerator = sp.expand(calculation['strong'].subs({
        m**2: source_mass_squared, s1: first_virtuality, s2: second_virtuality}))
    coefficient_row = sp.Matrix([[sp.diff(numerator, coefficient)
                                 for coefficient in sp.symbols('c_1:8 d_1:5', real=True)]])
    return numerator, sp.simplify(coefficient_row * calculation['nulls'])


def derive_tensor_response():
    """Derive the propagator on a PP current using the existing tensor Hessian."""
    from derive_alp_vector_two_pseudoscalar import derive_tensor_exchange
    hessian, response, expected = derive_tensor_exchange()
    return hessian, response, sp.simplify(response - expected)


def derive_pairing_signs():
    """Reduce each two-pair epsilon contraction to epsilon(p1,p2,p3,p4)."""
    # A symbolic 4x4 matrix supplies four independent external momenta;
    # determinant multilinearity gives the coefficient without imposing
    # masses or kinematic restrictions.
    momenta = sp.eye(4)
    results = []
    for second in (1, 2, 3):
        first_pair = (0, second)
        second_pair = tuple(index for index in range(4) if index not in first_pair)
        ordered = first_pair + second_pair
        coefficient = sp.det(sp.Matrix.vstack(*(momenta[index, :] for index in ordered)))
        results.append((tuple(index + 1 for index in ordered), coefficient))
    return results


def derive_scalar_mass_vertices():
    """Extract the fitted vector mass terms after removing the vector source.

    The change of flavour variables moves a pure-gradient c_q^VS source
    into delta chi_plus=-2 i a [c_q^VS,chi]/f_a. This is a representation
    change at first order in a, not an additional phenomenological operator.
    The invariant M_V and Delta_V masses have zero variation. Both retained
    quark-mass terms, including kappa_V, are varied explicitly.
    """
    source, mass, first, second = [sp.MatrixSymbol(name, 3, 3) for name in
                                 ('c_q^VS', 'chi', 'T_R', 'T_S')]
    alp_field, first_field, second_field = sp.symbols('a V_R V_S')
    alp_scale, mass_coefficient, double_trace_coefficient = sp.symbols(
        'f_a e_m^V kappa_V', real=True)
    vector = first_field * first + second_field * second
    delta_chi_plus = -2 * sp.I * alp_field * (source * mass - mass * source) / alp_scale
    mass_lagrangian = (-mass_coefficient * sp.Trace(vector * vector * delta_chi_plus)
                       + double_trace_coefficient * sp.Trace(vector)
                       * sp.Trace(vector * delta_chi_plus) / 4)
    # Polarization extracts a V_R V_S without differentiating matrix traces;
    # SymPy Trace's scalar derivative is unsupported for a product of traces.
    vertex = sum((-1)**(3 - sum(values)) * mass_lagrangian.subs(
        dict(zip((alp_field, first_field, second_field), values))).doit()
        for values in product((0, 1), repeat=3))
    return expand_traces(vertex)


def expand_traces(expression):
    """Apply SymPy's matrix expansion and trace linearity to display results."""
    return sp.expand(expression.replace(
        lambda term: isinstance(term, sp.Trace),
        lambda term: sp.Trace(sp.expand(term.arg)).doit()))


def derive_axial_vp_source():
    """Polarize i G_V Tr(V u_mu u_nu)/sqrt(2).

    Derivatives are kept as independent commuting field labels here; hence
    these are Lagrangian coefficients before Fourier or vertex factors.
    No second ALP rotation is applied to c_q^A.
    """
    vector, generator, axial = [sp.MatrixSymbol(name, 3, 3) for name in
                                ('T_R', 'lambda_i_dagger', 'c_q^A')]
    decay_constant, alp_scale, coupling = sp.symbols('f_pi f_a G_V', positive=True)
    alp_mu, alp_nu, meson_mu, meson_nu = sp.symbols(
        'partial_mu_a partial_nu_a partial_mu_P_i partial_nu_P_i')
    vielbein_mu = -generator*meson_mu/decay_constant + 2*axial*alp_mu/alp_scale
    vielbein_nu = -generator*meson_nu/decay_constant + 2*axial*alp_nu/alp_scale
    lagrangian = sp.I*coupling*sp.Trace(vector*(
        vielbein_mu*vielbein_nu))/sp.sqrt(2)
    expanded = expand_traces(lagrangian)
    first_coefficient = sp.expand(expanded).coeff(alp_mu).coeff(meson_nu)
    second_coefficient = sp.expand(expanded).coeff(alp_nu).coeff(meson_mu)
    # V^{mu nu}=-V^{nu mu}: relabel the second term's summed indices.
    coefficient = sp.expand(first_coefficient-second_coefficient)
    expected = -sp.sqrt(2)*sp.I*coupling*sp.Trace(
        vector*(axial*generator-generator*axial))/(decay_constant*alp_scale)
    return expanded, coefficient, sp.simplify(coefficient-expand_traces(expected))


def derive_axial_pseudoscalar_vv():
    """Insert the supplied no-mixing sources into the derived odd numerator."""
    first, second, axial, pseudoscalar, quark_mass = [sp.MatrixSymbol(name, 3, 3)
        for name in ('T_R', 'T_S', 'c_q^A', 'c_q^P', 'chi')]
    alp_scale = sp.Symbol('f_a', positive=True)
    numerator, _ = derive_odd_numerator()
    symmetric = first*second+second*first
    derivative_trace = sp.Trace(symmetric*axial)/alp_scale
    mass_trace = -sp.Trace(symmetric*(pseudoscalar*quark_mass
                                     +quark_mass*pseudoscalar))/alp_scale
    return expand_traces(numerator.subs({sp.Symbol('D', real=True): derivative_trace,
                                        sp.Symbol('P', real=True): mass_trace}))


def derive_absent_axial_pseudoscalar_pp():
    """Extract a Phi^2 from the LO cA,cP sources, keeping matrix order.

    The polynomial degree is sufficient: LO kinetic terms require the
    mesonic vielbein through Phi^3 and the axial one through Phi^2, while
    chi_plus contains the sine of its anticommutator map. The anomaly can
    alter bilinear mixing but cannot supply an a Phi^2 vertex.
    """
    field_order, alp_order = sp.symbols('t z')
    field, derivative, axial, mass_source = [sp.MatrixSymbol(name, 3, 3)
        for name in ('phi', 'partial_mu_phi', 'c_q^A', 'chi_P_over_a')]
    alp_derivative, alp_field, alp_scale, decay_constant = sp.symbols(
        'partial_mu_a a f_a f_pi', real=True)
    def commutator(first, second):
        return first*second-second*first
    def anticommutator(first, second):
        return first*second+second*first
    ordinary = -2*field_order*derivative + field_order**3*commutator(
        field, commutator(field, derivative))/3
    source = 2*alp_order*alp_derivative/alp_scale*(axial-field_order**2*
        commutator(field, commutator(field, axial))/2)
    kinetic = decay_constant**2*sp.Trace(ordinary*source+source*ordinary)/4
    mass = decay_constant**2*alp_order*alp_field*sp.Trace(
        2*field_order*anticommutator(field, mass_source)
        -field_order**3*anticommutator(field, anticommutator(
            field, anticommutator(field, mass_source)))/3)/4
    kinetic, mass = expand_traces(kinetic), expand_traces(mass)
    return kinetic, mass, (sp.expand(kinetic+mass).coeff(alp_order)
                           .coeff(field_order, 2))


def derive_axial_routing_signs():
    """Reduce epsilon(q_R,p_i,p_j,p_k) for q_R=p_a-p_i.

    P_i leaves the axial source, P_j leaves VVP, and (P_k,P_l)
    leave VPP. The sign multiplies epsilon(p1,p2,p3,p4).
    Determinants use momentum conservation before any mass-shell choice.
    """
    from itertools import permutations
    momenta = sp.eye(4)
    alp_momentum = sum((momenta[index, :] for index in range(4)), sp.zeros(1, 4))
    result = []
    for first, second in permutations(range(4), 2):
        third, fourth = [index for index in range(4) if index not in (first, second)]
        routed = alp_momentum-momenta[first, :]
        coefficient = sp.det(sp.Matrix.vstack(routed, momenta[first, :],
                                               momenta[second, :], momenta[third, :]))
        result.append(((first+1, second+1, third+1, fourth+1), coefficient))
    return result


def derive_contracted_amplitudes():
    """Sew the canonical odd vertex with classical full-tensor responses.

    Fourier convention: the ALP has exp(-i p_a x), outgoing mesons
    exp(+i p_i x). The odd-vertex normalization is inherited from the
    executable VVP operator contraction and its canonical VV amplitude:
    N epsilon(e_R,e_S,q_R,q_S)/(M_R M_S). Each tensor wave is
    i(q wedge e)/M. Matching that wave to a response includes its Fourier
    i; no additional propagator i or vertex i is inserted afterwards.

    Traces are scalar symbols here solely to display the contraction.
    The printed trace definitions specify the physical matrices.
    """
    coupling, decay_constant, alp_scale = sp.symbols('G_V f_pi f_a', positive=True)
    first_virtuality, second_virtuality = sp.symbols('s_R s_S', real=True)
    first_mass, second_mass, first_width, second_width = sp.symbols(
        'M_R M_S Gamma_R Gamma_S', positive=True)
    first_denominator = first_virtuality-first_mass**2+sp.I*first_mass*first_width
    second_denominator = second_virtuality-second_mass**2+sp.I*second_mass*second_width
    first_pair_trace, second_pair_trace, axial_trace = sp.symbols(
        'Tr_T_R_comm_lambda_i_lambda_j Tr_T_S_comm_lambda_k_lambda_l Tr_T_R_comm_cA_lambda_i')
    alp_numerator, meson_numerator = sp.symbols('N_a_RS N_j_RS')
    # iG/sqrt2 Tr(V u_mu u_nu), two outgoing derivatives give i*i=-1.
    pair_current = -sp.I*coupling*second_pair_trace/(2*sp.sqrt(2)*decay_constant**2)
    pair_response = 2*pair_current/second_denominator
    # wedge(pk+pl,pk)=-wedge(pk,pl), hence divide by -i.
    pair_polarization = sp.simplify(pair_response/(-sp.I))
    # The incoming ALP derivative (-i) and outgoing meson derivative (+i)
    # have product +1. Convert V pa pi to half its wedge before the Hessian.
    axial_current = -sp.sqrt(2)*sp.I*coupling*axial_trace/(2*decay_constant*alp_scale)
    axial_response = 2*axial_current/first_denominator
    # At the odd vertex the crossed vector momentum is -q_R and
    # wedge(-q_R,pi)=-wedge(pa,pi).
    axial_polarization = sp.simplify(axial_response/(-sp.I))
    first_pair_polarization = pair_polarization.subs({
        second_pair_trace:first_pair_trace, second_denominator:first_denominator})
    basis = sp.eye(4)
    first_vector = basis[0, :]+basis[1, :]
    second_vector = basis[2, :]+basis[3, :]
    pair_epsilon = sp.det(sp.Matrix.vstack(basis[0, :],basis[2, :],
                                           first_vector,second_vector))
    crossed_vector = -(basis[1, :]+basis[2, :]+basis[3, :])
    axial_epsilon = sp.det(sp.Matrix.vstack(basis[0, :],basis[2, :],
                                            crossed_vector,second_vector))
    pair_amplitude = sp.simplify(alp_numerator*first_pair_polarization*
                                 pair_polarization*pair_epsilon)
    axial_amplitude = sp.simplify(meson_numerator*axial_polarization*
                                  pair_polarization*axial_epsilon)
    return {'pair_current':pair_current, 'pair_response':pair_response,
            'axial_current':axial_current, 'axial_response':axial_response,
            'pair_epsilon':pair_epsilon, 'axial_epsilon':axial_epsilon,
            'pair_amplitude':pair_amplitude, 'axial_amplitude':axial_amplitude}


def show(label, expression):
    """Print deterministic symbolic expressions in the standard LaTeX format."""
    print(f'{label}:\n{sp.latex(expression)}\n')


@lru_cache(maxsize=1)
def derive_invariant_phase_space():
    """Derive the five-invariant measure and its one-variable integration.

    Let s=s12, t=s34, z1=s13+s14, z2=s13+s23 and u=s13.
    These are linear combinations of Mandelstam invariants, with unit
    Jacobian relative to (s12,s34,s13,s14,s23). The physical interval in u
    is fixed by -det(p_i.p_j)>=0. Both orientations of the final momenta
    are included; a parity-conserving pseudoscalar amplitude squared is
    identical for them. There is no identical-particle factor in dPhi4.

    The last integration is analytical whenever the scalar amplitude is
    independent of u. Polynomial moments are also provided exactly.
    For the full crossed vector-exchange sum this independence is generally
    absent, so the measure alone does not reduce five integrations to four.
    """
    alp_mass_squared = sp.Symbol('m_a_squared', positive=True)
    mass_squared = sp.symbols('m1_squared m2_squared m3_squared m4_squared', nonnegative=True)
    first_pair, second_pair, first_sum, second_sum, crossed_pair = sp.symbols(
        's12 s34 s13_plus_s14 s13_plus_s23 s13', real=True)
    pair_invariants = {(0, 1): first_pair, (2, 3): second_pair,
                       (0, 2): crossed_pair, (0, 3): first_sum-crossed_pair,
                       (1, 2): second_sum-crossed_pair}
    pair_invariants[1, 3] = (alp_mass_squared+2*sum(mass_squared)
                             -sum(pair_invariants.values()))
    gram = sp.diag(*mass_squared)
    for (first, second), invariant in pair_invariants.items():
        gram[first, second] = gram[second, first] = (
            invariant-mass_squared[first]-mass_squared[second])/2
    epsilon_squared = sp.Poly(-gram.det(method='domain-ge'), crossed_pair)
    quadratic_coefficient = sp.factor(epsilon_squared.coeff_monomial(crossed_pair**2))
    linear_coefficient = epsilon_squared.coeff_monomial(crossed_pair)
    constant_coefficient = epsilon_squared.coeff_monomial(1)
    kallen = sp.expand(alp_mass_squared**2+first_pair**2+second_pair**2
                       -2*alp_mass_squared*first_pair-2*alp_mass_squared*second_pair
                       -2*first_pair*second_pair)
    def pair_kallen(invariant, first_mass_squared, second_mass_squared):
        return invariant**2+first_mass_squared**2+second_mass_squared**2-2*(
            invariant*first_mass_squared+invariant*second_mass_squared
            +first_mass_squared*second_mass_squared)
    first_sum_centre = (2*mass_squared[0]+mass_squared[2]+mass_squared[3]
        +(first_pair+mass_squared[0]-mass_squared[1])
        *(alp_mass_squared-first_pair-second_pair)/(2*first_pair))
    second_sum_centre = (mass_squared[0]+mass_squared[1]+2*mass_squared[2]
        +(second_pair+mass_squared[2]-mass_squared[3])
        *(alp_mass_squared-first_pair-second_pair)/(2*second_pair))
    first_sum_half_range = sp.sqrt(kallen*pair_kallen(
        first_pair, mass_squared[0], mass_squared[1]))/(2*first_pair)
    second_sum_half_range = sp.sqrt(kallen*pair_kallen(
        second_pair, mass_squared[2], mass_squared[3]))/(2*second_pair)
    centre = -linear_coefficient/(2*quadratic_coefficient)
    radius_squared = (linear_coefficient**2-4*quadratic_coefficient*constant_coefficient)/(4*quadratic_coefficient**2)
    # Recursive dPhi4 = ds dt/(2pi)^2 times three dPhi2 factors.
    # Integrate the overall orientation (4pi) and common azimuth (2pi).
    angular_coefficient = sp.simplify((4*sp.pi)*(2*sp.pi)/(2*sp.pi)**2/(32*sp.pi**2)**3)
    # The Jacobian to five invariants is sqrt(lambda_a lambda_12 lambda_34)
    # *sqrt(epsilon_squared)/(s*t). The angular measure has that same
    # first factor/(s*t); phi and 2pi-phi give two inverse images.
    invariant_coefficient = sp.simplify(2*angular_coefficient/alp_mass_squared)
    angle = sp.Symbol('theta', real=True)
    angular_moments = tuple(sp.integrate(sp.sin(angle)**2*sp.cos(angle)**order,
                                        (angle, 0, sp.pi)) for order in range(5))
    # epsilon² supplies lambda_a lambda_12 lambda_34/(64*s*t).
    # The full relative azimuth is 0..2pi, so its sin² integral is pi.
    # Include the decay flux 1/(2ma), separately from dPhi4.
    single_pair_width_coefficient = sp.simplify(
        angular_coefficient*sp.Rational(1, 64)*sp.Rational(1, 2)
        *sp.Rational(4, 3)**2*sp.pi)
    return {'gram': gram, 'pair_invariants': pair_invariants,
            'epsilon_squared': epsilon_squared.as_expr(),
            'quadratic_coefficient': quadratic_coefficient,
            'quadratic_coefficient_residual': sp.factor(quadratic_coefficient+kallen/16),
            'centre': centre, 'radius_squared': radius_squared,
            'first_sum_bounds': (first_sum_centre-first_sum_half_range,
                                 first_sum_centre+first_sum_half_range),
            'second_sum_bounds': (second_sum_centre-second_sum_half_range,
                                  second_sum_centre+second_sum_half_range),
            'angular_coefficient': angular_coefficient,
            'invariant_coefficient': invariant_coefficient,
            'angular_moments': angular_moments,
            'single_pair_width_coefficient': single_pair_width_coefficient}


def print_invariant_phase_space():
    """Display Mandelstam-only width and exact conditional reduction."""
    result = derive_invariant_phase_space()
    show('Gram matrix p_i dot p_j; s24 fixed by sum sij=ma²+2 sum mi²', result['gram'])
    show('Coefficient of s13 squared in -det Gram', result['quadratic_coefficient'])
    show('Angular phase-space coefficient before sqrt(lambda_a lambda_12 lambda_34)/(ma² s12 s34)',
         result['angular_coefficient'])
    show('Invariant phase-space coefficient before ds12 ds34 dz1 dz2 ds13/sqrt(-det Gram)',
         result['invariant_coefficient'])
    show('Integral_0^pi sin(theta)^2 cos(theta)^n dtheta, n=0,...,4',
         sp.Tuple(*result['angular_moments']))
    show('Independent bounds for s13+s14 at fixed s12,s34',
         sp.Tuple(*result['first_sum_bounds']))
    show('Independent bounds for s13+s23 at fixed s12,s34',
         sp.Tuple(*result['second_sum_bounds']))
    show('Single-pair angular integration coefficient including decay flux',
         result['single_pair_width_coefficient'])
    print('For M=epsilon(p1,p2,p3,p4) F(sij), |M|²=(-det Gram)|F|². '
          'Gamma=1/(16384 pi^6 ma^3 S) integral ds12 ds34 dz1 dz2 ds13 '
          'sqrt(-det Gram)|F|², where S=product n_identical!. '
          's12 lies between (m1+m2)^2 and (ma-m3-m4)^2; '
          's34 lies between (m3+m4)^2 and (ma-sqrt(s12))^2. '
          'The two sums lie in the displayed independent intervals. '
          'At fixed s12,s34,z1,z2, write -det Gram=lambda(ma²,s12,s34)/16 '
          '*[r²-(s13-c)²], with c=-B/(2A), r²=(B²-4AC)/(4A²), '
          'where A,B,C are its quadratic coefficients. '
          'Then s13=c+r cos(theta), theta in [0,pi]. '
          'If F is independent of s13 the last integral equals '
          'pi sqrt(lambda) r² |F|²/8. '
          'For polynomial |F|² use the printed moments after substitution. '
          'For rational tree amplitudes the same substitution makes this '
          'a rational trigonometric integral and it can in principle be evaluated '
          'analytically term by term. The adopted constant PDG widths place '
          'the vector poles away from the real integration interval. '
          'A single aVV pairing has F(s12,s34): all three angles integrate exactly, '
          'giving dGamma/(ds12 ds34)=lambda_a^(3/2) lambda_12^(3/2) '
          'lambda_34^(3/2)|F|²/(1179648 pi^5 ma^3 s12² s34² S). '
          'Crossed-pair interference and the axial aVP family generally retain '
          'all five invariant dependences. No narrow-width approximation is used.')


def derive_semicircle_integrals():
    """Derive the rational-amplitude inner integral using a semicircle weight.

    With u=c+x, the invariant measure contains sqrt(r²-x²). Polynomial
    division and partial fractions therefore reduce any rational amplitude
    squared with no poles on the interval to polynomial moments and the
    displayed Cauchy integral. Its square root is analytic off [-r,r] and
    tends to z at infinity. Repeated poles follow by differentiating in z.
    This includes the finite complex resonance poles from the adopted
    a -> 3P constant-width prescription; the integral identity itself does
    not determine the physical widths.
    """
    radius, pole_distance = sp.symbols('r z', positive=True)
    tangent, positive_ratio = sp.symbols('t b', positive=True)
    # For real z>r the tangent substitution gives
    # int_0^pi dtheta/(z-r cos theta)
    # =2/(z+r) int_0^infty dt/(t²+(z-r)/(z+r)).
    elementary_integral = sp.integrate(1/(tangent**2+positive_ratio**2),
                                        (tangent, 0, sp.oo))
    inverse_linear = sp.simplify(2/(pole_distance+radius)*elementary_integral.subs(
        positive_ratio, sp.sqrt((pole_distance-radius)/(pole_distance+radius))))
    cosine = sp.Symbol('cos_theta', real=True)
    quotient, remainder = sp.div(radius**2*(1-cosine**2),
                                  radius*cosine-pole_distance, cosine)
    # Polynomial quotient is -r cos(theta)-z; its integral is -pi*z.
    cauchy_integral = sp.pi*(-pole_distance+sp.sqrt(pole_distance**2-radius**2))
    return {'elementary_integral': elementary_integral,
            'inverse_linear_integral': inverse_linear,
            'quotient': quotient, 'remainder': remainder,
            'division_residual': sp.expand(radius**2*(1-cosine**2)
                -(radius*cosine-pole_distance)*quotient-remainder),
            'cauchy_integral': cauchy_integral,
            'first_five_polynomial_moments': tuple(
                0 if order % 2 else sp.pi*radius**(order+2)*sp.catalan(order//2)
                /2**(order+1) for order in range(5))}


def print_semicircle_integrals():
    """Explain the exact one-invariant reduction for rational amplitudes."""
    result = derive_semicircle_integrals()
    for label, expression in result.items():
        show('Semicircle rational integration: '+label, sp.Tuple(*expression)
             if isinstance(expression, tuple) else expression)
    print('For a rational |F(u)|² write R(c+x)=P(x)+sum a_(z,k)/(x-z)^k '
          'by polynomial division and partial fractions. The inner integral is '
          'sqrt(lambda_a)/4 times [sum P_n mu_n + '
          'sum a_(z,k) (d/dz)^(k-1) J(z)/(k-1)!], '
          'mu_(2n)=pi r^(2n+2) Catalan(n)/2^(2n+1), mu_(2n+1)=0, '
          'J(z)=pi[sqrt(z²-r²)-z], sqrt(z²-r²)~z at infinity. '
          'Hence the complete rational tree amplitude permits exactly four '
          'remaining numerical integrations, including crossed interference. '
          'The adopted a -> 3P constant PDG widths leave no vector pole on '
          'the real integration interval and preserve rationality. General running widths '
          'or nonrational form factors need separate analysis. Near coincident '
          'poles use repeated-pole limits to avoid numerical cancellation.')


@lru_cache(maxsize=1)
def derive_flavour_selection_rules():
    """Enumerate nonzero flavour routings before summing Lorentz amplitudes.

    Charged generators at external outgoing legs are conjugated. The two
    ends of a charged propagator carry mutually conjugate generators.
    Neutral mixing entries remain symbolic: a zero means a polynomial
    identity, not cancellation at a chosen angle or fitted coefficient.
    All three neutral pseudoscalars have independent symbolic u,d,s entries.
    This retains their isospin admixtures: the audit is generic support before
    inserting the Alda:2026nrd fitted neutral generators. Independent entries
    are only a test of polynomial support, not a new mixing prescription.
    Orthogonality or fitted-parameter cancellations require substitution.
    Traces alone cannot establish cancellation between distinct routings.
    """
    def elementary(row, column):
        matrix = sp.zeros(3)
        matrix[row, column] = 1
        return matrix

    omega_n, omega_s, phi_n, phi_s = sp.symbols('omega_n omega_s phi_n phi_s', real=True)
    pi0_u, pi0_d, pi0_s = sp.symbols('pi0_u pi0_d pi0_s', real=True)
    eta_u, eta_d, eta_s = sp.symbols('eta_u eta_d eta_s', real=True)
    etap_u, etap_d, etap_s = sp.symbols('etap_u etap_d etap_s', real=True)
    nonstrange = sp.diag(1, 1, 0)/sp.sqrt(2)
    strange = sp.diag(0, 0, 1)
    vectors = {'rho0':sp.diag(1, -1, 0)/sp.sqrt(2),
               'rho+':elementary(0, 1), 'rho-':elementary(1, 0),
               'K*+':elementary(0, 2), 'K*-':elementary(2, 0),
               'K*0':elementary(1, 2), 'K*0bar':elementary(2, 1),
               'omega':omega_n*nonstrange+omega_s*strange,
               'phi':phi_n*nonstrange+phi_s*strange}
    mesons = {'pi0':sp.diag(pi0_u, pi0_d, pi0_s),
              'pi+':sp.sqrt(2)*elementary(1, 0),
              'pi-':sp.sqrt(2)*elementary(0, 1),
              'K+':sp.sqrt(2)*elementary(2, 0),
              'K-':sp.sqrt(2)*elementary(0, 2),
              'K0':sp.sqrt(2)*elementary(2, 1),
              'K0bar':sp.sqrt(2)*elementary(1, 2),
              'eta':sp.diag(eta_u, eta_d, eta_s),
              'etap':sp.diag(etap_u, etap_d, etap_s)}
    sources = {'diagonal':sp.diag(*sp.symbols('c_u c_d c_s')),
               'ds':elementary(1, 2), 'sd':elementary(2, 1)}
    quark_mass = sp.diag(*sp.symbols('chi_u chi_d chi_s', positive=True))
    pair_traces = {}
    for first, second in combinations_with_replacement(mesons, 2):
        commutator = mesons[first]*mesons[second]-mesons[second]*mesons[first]
        pair_traces[first, second] = {name:sp.expand(sp.trace(matrix*commutator))
            for name, matrix in vectors.items() if sp.expand(sp.trace(matrix*commutator)) != 0}
    odd_traces, axial_traces, strong_traces = {}, {}, {}
    for first, first_matrix in vectors.items():
        for second, second_matrix in vectors.items():
            source_product = first_matrix.adjoint()*second_matrix.adjoint()+second_matrix.adjoint()*first_matrix.adjoint()
            strong_product = first_matrix*second_matrix.adjoint()+second_matrix.adjoint()*first_matrix
            for sector, source in sources.items():
                odd_traces[sector, first, second] = (
                    sp.expand(sp.trace(source_product*source)),
                    sp.expand(sp.trace(source_product*(source*quark_mass+quark_mass*source))))
            for name, generator in mesons.items():
                strong_traces[first, second, name] = (
                    sp.expand(sp.trace(strong_product*generator)),
                    sp.expand(sp.trace(strong_product*(generator*quark_mass+quark_mass*generator))))
        for sector, source in sources.items():
            for name, generator in mesons.items():
                axial_traces[sector, first, name] = sp.expand(sp.trace(
                    first_matrix.adjoint()*(source*generator-generator*source)))
    channels = {sector:{} for sector in sources}
    for final_state in combinations_with_replacement(mesons, 4):
        for sector in sources:
            pair_routes, axial_routes = set(), set()
            for first_indices in combinations(range(4), 2):
                second_indices = tuple(index for index in range(4) if index not in first_indices)
                first_pair = tuple(final_state[index] for index in first_indices)
                second_pair = tuple(final_state[index] for index in second_indices)
                for first_vector in pair_traces[first_pair]:
                    for second_vector in pair_traces[second_pair]:
                        derivative, mass = odd_traces[sector, first_vector, second_vector]
                        if derivative != 0 or mass != 0:
                            pair_routes.add((first_pair, second_pair, first_vector, second_vector))
                for first_index in second_indices:
                    second_index = next(index for index in second_indices if index != first_index)
                    for second_vector in pair_traces[first_pair]:
                        for first_vector in vectors:
                            if (axial_traces[sector, first_vector, final_state[first_index]] != 0
                                and any(trace != 0 for trace in strong_traces[
                                    first_vector, second_vector, final_state[second_index]])):
                                axial_routes.add((final_state[first_index], final_state[second_index],
                                                  first_pair, first_vector, second_vector))
            if pair_routes or axial_routes:
                channels[sector][final_state] = {'aVV':sorted(pair_routes), 'aVP':sorted(axial_routes)}
    return {'pair_traces':pair_traces, 'odd_traces':odd_traces,
            'axial_traces':axial_traces, 'strong_traces':strong_traces,
            'channels':channels}


def print_flavour_selection_rules():
    """Print all structural routings as deterministic JSON, without fit cuts."""
    import json
    result = derive_flavour_selection_rules()
    print(json.dumps({sector:[{'final_state':state, **routes}
                             for state, routes in channels.items()]
                      for sector, channels in result['channels'].items()}, indent=2))


def main():
    """Print the field-counting derivation and independently reduced vertices."""
    print(__doc__)
    for vertices, edges, species, external_counts in enumerate_vector_trees():
        print('Tree:', ' + '.join(vertex.name for vertex in vertices))
        print('Internal lines:', tuple((vertices[left].name, vertices[right].name, particle)
                                       for (left, right), particle in zip(edges, species)))
        print('External P at each vertex:', external_counts)
    numerator, projection = derive_odd_numerator()
    show('VVP numerator: D=Tr({T_R,T_S} U_X)/2; P=Tr({T_R,T_S} Y_X)/2', numerator)
    show('Original preliminary-fit coefficient null projection', projection)
    _, response, residual = derive_tensor_response()
    show('Full tensor propagator response on J times p1 wedge p2', response)
    show('Response minus 2 J (p1 wedge p2)/(q squared - M_R squared)', residual)
    for pairing, sign in derive_pairing_signs():
        show(f'Epsilon ordering {pairing} relative to (1,2,3,4)', sign)
    from derive_alp_vector_two_pseudoscalar import derive_general_vector_source_vertices
    show('Rotated scalar-source aPP Lagrangian coefficient',
         expand_traces(derive_general_vector_source_vertices()))
    show('Rotated scalar-source aVV Lagrangian coefficient, contracted Lorentz indices',
         derive_scalar_mass_vertices())
    axial_lagrangian, axial_coefficient, axial_residual = derive_axial_vp_source()
    show('G_V axial-source polarized Lagrangian', axial_lagrangian)
    show('Coefficient of partial_mu_a partial_nu_P_i in the axial source', axial_coefficient)
    show('Axial-source coefficient minus -sqrt(2) i G_V Tr(T_R[c_q^A,lambda_i_dagger])/(f_pi f_a)', axial_residual)
    kinetic, mass, absent_pp = derive_absent_axial_pseudoscalar_pp()
    show('LO axial kinetic source through t cubed, t counts Phi', kinetic)
    show('LO pseudoscalar mass source through t cubed, t counts Phi', mass)
    show('LO aPP coefficient for c_q^A and c_q^P', absent_pp)
    show('Odd aVV numerator with explicit no-mixing axial and pseudoscalar sources',
         derive_axial_pseudoscalar_vv())
    for assignment, sign in derive_axial_routing_signs():
        show(f'epsilon(q_R,p_i,p_j,p_k), (i,j,k,l)={assignment}', sign)
    for label, expression in derive_contracted_amplitudes().items():
        show('Sewn canonical scalar amplitude: '+label, expression)
    print('For the pair family sum the displayed pair_amplitude over '
          '(ij)(kl)=(12)(34),(13)(24),(14)(23), with permutation signs '
          '+,-,+, and over real orthonormal vector generators R,S. '
          'For the axial family sum axial_amplitude over ordered distinct '
          'i,j, put the remaining labels in k<l, multiply by the parity of '
          '(i,j,k,l), and sum R,S. The common Lorentz factor is '
          'epsilon(p1,p2,p3,p4). N_a_RS is the odd numerator with '
          'U=2cA/fa,Y=-2{cP,chi}/fa and qR²=sij,qS²=skl. '
          'N_j_RS uses U=-lambda_j_dagger/fpi,Y={lambda_j_dagger,chi}/fpi, '
          'source mass squared=mj², qR²=(pa-pi)²,qS²=(pk+pl)². '
          'N already includes both vector assignments; do not double it. '
          'For every internal vector use 1/(s-M_R²+i M_R Gamma_R), '
          'with the appropriate two- or three-particle invariant s. This is '
          'minus 1/(M_R²-s-i M_R Gamma_R), the denominator convention used '
          'in a -> 3P. Adopt the same constant PDG total widths in GeV: '
          'rho0 0.1474; rho+ and rho- 0.1491; K*+ and K*- 0.0514; '
          'K*0 and K*0bar 0.0473; omega 0.00868; phi 0.004249. '
          'Neutral rho0, omega and phi masses remain the fitted values. '
          'The other masses in GeV match a -> 3P: rho+ and rho- 0.77511; '
          'K*+ and K*- 0.89167; K*0 and K*0bar 0.89555. '
          'These widths are a phenomenological prescription applied after '
          'the zero-width tensor contraction, with no running-width or '
          'narrow-width approximation.')
    print('For c_q^A,c_q^P only, the aPP+VPP+VPP family is absent. '
          'Its LO source coefficient is zero before introducing an internal '
          'pseudoscalar propagator. The remaining families are odd aVV '
          'and axial aVP+VVP+VPP. The polarized axial source is a '
          'Lagrangian coefficient; the sewn expressions above include '
          'Fourier factors through matching to the canonical odd amplitude.')
    print('aVV has separate odd VVP-source and even vector-source tensors. '
          'aVP is the axial G_V source; aVPP and even aVV arise from the '
          'vector source. aPP followed by VPP and VPP has one internal '
          'pseudoscalar and requires the LO chiral source vertex. '
          'The higher field terms in the same retained interactions cannot '
          'contribute additional vector trees at this external multiplicity. '
          'Pure pseudoscalar trees and electromagnetic exchanges are outside '
          'this vector-exchange subset. Each allowed flavour generator and '
          'each inequivalent assignment of four external labelled legs must '
          'be summed. Identical-particle phase-space factors do not multiply '
          'the amplitude.')
    print_invariant_phase_space()
    print_semicircle_integrals()


if __name__ == '__main__':
    import sys
    import unittest

    class DerivationTests(unittest.TestCase):
        """Focused checks of topology completeness and off-shell contractions."""

        def test_tree_inventory(self):
            expected = {('aPP', 'VPP', 'VPP'), ('aVP', 'VPP', 'VVP'),
                        ('aVPP', 'VPP'), ('aVV', 'VPP', 'VPP')}
            self.assertEqual({tuple(vertex.name for vertex in result[0])
                              for result in enumerate_vector_trees()}, expected)

        def test_vpp_flavour_selection(self):
            result = derive_flavour_selection_rules()
            self.assertEqual(set(result['pair_traces']['pi+', 'pi-']), {'rho0'})
            self.assertEqual(result['pair_traces']['pi0', 'eta'], {})
            self.assertEqual(result['pair_traces']['K0', 'K0'], {})
            self.assertEqual(set(result['pair_traces']['K+', 'K-']), {'rho0', 'omega', 'phi'})

        def test_axial_isoscalar_and_pure_neutral_exclusions(self):
            result = derive_flavour_selection_rules()
            self.assertEqual(result['axial_traces']['diagonal', 'rho0', 'pi0'], 0)
            for state in result['channels']['diagonal']:
                self.assertFalse(all(name in {'pi0', 'eta', 'etap'} for name in state))
            self.assertIn(('pi+', 'pi+', 'pi-', 'pi-'), result['channels']['diagonal'])

        def test_source_quantum_numbers_and_explicit_isospin_limit(self):
            result = derive_flavour_selection_rules()
            charge = {'pi+':1, 'pi-':-1, 'K+':1, 'K-':-1}
            strangeness = {'K+':1, 'K-':-1, 'K0':1, 'K0bar':-1}
            for sector, expected_strangeness in (('diagonal', 0), ('ds', 1), ('sd', -1)):
                for state in result['channels'][sector]:
                    self.assertEqual(sum(charge.get(name, 0) for name in state), 0)
                    self.assertEqual(sum(strangeness.get(name, 0) for name in state), expected_strangeness)
            derivative, mass = result['strong_traces']['rho+', 'rho+', 'pi0']
            self.assertNotEqual(derivative, 0)
            self.assertNotEqual(mass, 0)
            pi0_u, pi0_d, pi0_s = sp.symbols('pi0_u pi0_d pi0_s', real=True)
            isospin_limit = {pi0_u:1, pi0_d:-1, pi0_s:0}
            self.assertEqual(derivative.subs(isospin_limit), 0)
            chi_u, chi_d = sp.symbols('chi_u chi_d', positive=True)
            self.assertNotEqual(mass.subs(isospin_limit), 0)
            self.assertEqual(mass.subs(isospin_limit).subs(chi_u, chi_d), 0)

        def test_pairing_parity(self):
            self.assertEqual([sign for _, sign in derive_pairing_signs()], [1, -1, 1])

        def test_off_shell_null_projection(self):
            self.assertEqual(derive_odd_numerator()[1], sp.zeros(1, 4))

        def test_full_tensor_response(self):
            self.assertEqual(derive_tensor_response()[2], sp.zeros(6, 1))

        def test_conserved_vector_source_mass_vertex(self):
            source = sp.MatrixSymbol('c_q^VS', 3, 3)
            self.assertEqual(sp.simplify(derive_scalar_mass_vertices().subs(
                source, sp.Identity(3)).doit()), 0)

        def test_axial_source_antisymmetry(self):
            self.assertEqual(derive_axial_vp_source()[2], 0)

        def test_no_axial_or_pseudoscalar_a_pp(self):
            self.assertEqual(derive_absent_axial_pseudoscalar_pp()[2], 0)

        def test_axial_routing(self):
            self.assertEqual(len(derive_axial_routing_signs()), 12)
            self.assertTrue(all(abs(sign) == 1 for _, sign in derive_axial_routing_signs()))

        def test_sewn_epsilon_signs(self):
            sewn = derive_contracted_amplitudes()
            self.assertEqual(sewn['pair_epsilon'], -1)
            self.assertEqual(sewn['axial_epsilon'], 1)

        def test_gram_quadratic_and_phase_space_normalization(self):
            result = derive_invariant_phase_space()
            self.assertEqual(result['quadratic_coefficient_residual'], 0)
            self.assertEqual(result['angular_coefficient'], 1/(16384*sp.pi**6))
            self.assertEqual(result['invariant_coefficient'],
                1/(8192*sp.pi**6*sp.Symbol('m_a_squared', positive=True)))
            self.assertEqual(result['angular_moments'],
                             (sp.pi/2, 0, sp.pi/8, 0, sp.pi/16))

        def test_semicircle_partial_fraction_integrals(self):
            result = derive_semicircle_integrals()
            self.assertEqual(result['division_residual'], 0)
            # An independent convergent geometric expansion at large z
            # reproduces the first nonzero semicircle moments.
            inverse_distance = sp.Symbol('w', positive=True)
            radius, distance = sp.symbols('r z', positive=True)
            expanded = sp.series(result['cauchy_integral'].subs(distance, 1/inverse_distance),
                                  inverse_distance, 0, 7).removeO().expand()
            for order, moment in enumerate(result['first_five_polynomial_moments']):
                self.assertEqual(sp.simplify(expanded.coeff(inverse_distance, order+1)+moment), 0)

    if '--selection-rules' in sys.argv:
        print_flavour_selection_rules()
    elif '--phase-space' in sys.argv:
        print_invariant_phase_space()
        print_semicircle_integrals()
    elif '--test' in sys.argv:
        unittest.main(argv=[sys.argv[0]])
    else:
        main()
