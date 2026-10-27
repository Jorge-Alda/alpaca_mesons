"""Derive the Born overlap with the geometric three-pion spin-one current.

Run: .venv/bin/python scripts/derive_alp_three_pion_born_spin_one.py
Inputs: memory/alp/lepton-pair/alp-three-pion-isoscalar-born-projection.md
and alp-three-pion-discontinuity.md. Natural units, metric (+---).
The normalized SO(3) orientation average holds at fixed Dalitz invariants.
This is the component selected by the vacuum current, not a complete
helicity decomposition of the ALP production amplitude.
"""

import sympy as sp


def derive():
    """Calculate azimuthal overlap, polar integral and spin-one reconstruction."""
    azimuth = sp.Symbol('varphi', real=True)
    cosine = sp.Symbol('z', real=True)
    photon_momentum = sp.Symbol('abs_q', positive=True)
    pion_momentum = sp.Symbol('abs_p_r', positive=True)
    normal_magnitude = sp.Symbol('abs_J_3pi', positive=True)
    transverse_difference, longitudinal_difference = sp.symbols('v_perp v_parallel', real=True)
    energy_difference, alp_energy = sp.symbols('E_s_minus_E_t E_a', real=True)
    normal = sp.Matrix([0, normal_magnitude, 0])
    pion = sp.Matrix([0, 0, pion_momentum])
    difference = sp.Matrix([transverse_difference, 0, longitudinal_difference])
    photon = photon_momentum * sp.Matrix([
        sp.sqrt(1-cosine**2)*sp.cos(azimuth),
        sp.sqrt(1-cosine**2)*sp.sin(azimuth), cosine])
    source_difference = alp_energy*energy_difference + photon.dot(difference)
    azimuth_integrand = normal.dot(pion.cross(photon))*source_difference
    azimuth_average = sp.simplify(sp.integrate(azimuth_integrand, (azimuth, 0, 2*sp.pi))/(2*sp.pi))
    geometric_product = normal.dot(pion.cross(difference))
    azimuth_ratio = sp.factor(azimuth_average/geometric_product)

    # Orient the fixed Dalitz triangle in a plane; SO(3) restores all its
    # orientations. Matrix cross products derive the three cyclic signs.
    total_energy = sp.Symbol('sqrt_x', positive=True)
    first_x, second_x, second_y = sp.symbols('p_plus_x p_minus_x p_minus_y', real=True)
    positive = sp.Matrix([first_x, 0, 0])
    negative = sp.Matrix([second_x, second_y, 0])
    neutral = -positive-negative
    # epsilon^{0123}=+1 with three lowered pion momenta gives this minus sign.
    geometric_current = -total_energy*positive.cross(negative)
    pion_energies = sp.symbols('E_plus E_minus E_zero', real=True)
    four_vectors = [sp.Matrix([energy, *momentum]) for energy, momentum in
                    zip(pion_energies, (positive, negative, neutral))]
    metric = sp.diag(1, -1, -1, -1)
    lowered = [metric*vector for vector in four_vectors]
    epsilon_current = sp.Matrix([sum(sp.LeviCivita(mu, first, second, third)
                                     *lowered[0][first]*lowered[1][second]*lowered[2][third]
                                     for first in range(4) for second in range(4) for third in range(4))
                                 for mu in range(4)])
    epsilon_difference = sp.simplify(epsilon_current.subs(pion_energies[2], total_energy-sum(pion_energies[:2]))
                                    - sp.Matrix([0, *geometric_current]))
    cyclic_products = tuple(sp.factor(geometric_current.dot(first.cross(second-third)))
                            for first, second, third in (
                                (positive, neutral, negative),
                                (negative, positive, neutral),
                                (neutral, negative, positive)))
    cyclic_ratios = tuple(sp.factor(product/geometric_current.dot(geometric_current))
                          for product in cyclic_products)

    denominator_constant, denominator_linear = sp.symbols('a b', positive=True)
    polar_integrand = (1-cosine**2)/(denominator_constant+denominator_linear*cosine)
    quotient, remainder = sp.div(1-cosine**2, denominator_constant+denominator_linear*cosine, cosine)
    primitive = sp.integrate(quotient, cosine) + remainder*sp.log(denominator_constant+denominator_linear*cosine)/denominator_linear
    polar_integral = sp.factor(primitive.subs(cosine, 1)-primitive.subs(cosine, -1))
    zero_momentum_limit = sp.simplify(sp.limit(polar_integral, denominator_linear, 0))

    source = sp.MatrixSymbol('c_q^A', 3, 3)
    neutral_generator = sp.MatrixSymbol('T_0', 3, 3)
    electric_charge, alp_decay_constant, pion_decay_constant = sp.symbols('e f_a f_pi', nonzero=True)
    born_coefficient = 4*electric_charge*sp.Trace(source*neutral_generator)/(3*alp_decay_constant*pion_decay_constant)
    # The Born photon vertex contributes 2*p_r; its -q part drops out of
    # J_3pi dot (J_Born cross q). The remaining polar average is integral/2.
    polar_integrals = sp.symbols('I_plus I_minus I_zero')
    overlap = sp.factor(sum(2*born_coefficient*photon_momentum**2/2
                            *normal_magnitude**2*ratio*integral/2
                            for ratio, integral in zip(cyclic_ratios, polar_integrals)))
    tensor_coefficient = sp.factor(overlap/(2*photon_momentum**2))
    reconstruction_coefficient = sp.factor(-3*tensor_coefficient/normal_magnitude**2)
    return dict(azimuth_integrand=azimuth_integrand, azimuth_average=azimuth_average,
                azimuth_ratio=azimuth_ratio, epsilon_difference=epsilon_difference,
                cyclic_ratios=cyclic_ratios,
                polar_integrand=polar_integrand, quotient=quotient, remainder=remainder,
                primitive=primitive, polar_integral=polar_integral,
                zero_momentum_limit=zero_momentum_limit, born_coefficient=born_coefficient,
                overlap=overlap, tensor_coefficient=tensor_coefficient,
                reconstruction_coefficient=reconstruction_coefficient)


def main():
    """Print the derivation with the angular and analytic-continuation conventions."""
    print('Exact isospin; e retained; fixed Dalitz invariants; normalized SO(3) average.')
    print('Q=p+q=(sqrt(x),0); J_3pi spatial=-sqrt(x)*(p_plus cross p_minus).')
    print('C=4e Trace(c_q^A T_0)/(3 f_a f_pi). Cyclic order (r,s,t)=(+,0,-),(-,+,0),(0,-,+).')
    results = derive()
    for name, expression in results.items():
        print(f'{name}:\n{sp.latex(expression)}\n')
    print('a=q^2-2 q^0 E_r; b=2 abs(q) abs(p_r). Integral first derived for a>b>0.')
    print('Continue from that domain with the original D_r+i0; do not choose endpoint logs independently.')
    print('overlap=<J_3pi dot (P0 J_Born cross q)>; B_ij=<J_3pi_i P0 J_Born_j>.')
    print('B_ij=tensor_coefficient epsilon_ijk q_k; reconstructed current=coefficient*(J_3pi cross q).')
    print('The projected time component is zero. This selects the normal J=1 channel only.')


if __name__ == '__main__':
    main()
