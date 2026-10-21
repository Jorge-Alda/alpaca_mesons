"""Derive the invariant a -> V P1 P2 width from the recorded tree amplitude.

Run: .venv/bin/python scripts/derive_alp_vector_two_pseudoscalar_width.py
Natural units; masses and momenta in GeV; metric (+---), epsilon^0123=+1.
m_V is the physical external mass; M_V in the parent amplitude is a separate
Lagrangian parameter. Internal line shapes remain amplitude inputs.
Only |Delta S|=1 entries of c_q^VS enter the even amplitude. No narrow-width
approximation is made for internal lines. See memory/alp/three-meson/
alp-vpp-amplitude.md for the amplitude coefficients and truncation.
"""
import sympy as sp


def derive_width():
    """Derive invariant contractions, phase-space boundaries and normalization."""
    alp_mass, vector_mass, first_mass, second_mass = sp.symbols(
        'm_a m_V m_1 m_2', positive=True)
    pair_mass_squared, second_invariant, third_invariant = sp.symbols('s t u', real=True)
    cosine = sp.Symbol('z', real=True)
    identical_factor = sp.Symbol('S_PP', positive=True)
    odd_coefficient = sp.Function('F_odd')(pair_mass_squared, second_invariant)
    first_coefficient = sp.Function('F_1')(pair_mass_squared, second_invariant)
    second_coefficient = sp.Function('F_2')(pair_mass_squared, second_invariant)
    mandelstam_sum = alp_mass**2 + vector_mass**2 + first_mass**2 + second_mass**2
    third_replacement = sp.solve(sp.Eq(pair_mass_squared + second_invariant +
                                      third_invariant, mandelstam_sum), third_invariant)[0]
    vector_first_product = (third_replacement-vector_mass**2-first_mass**2)/2
    vector_second_product = (second_invariant-vector_mass**2-second_mass**2)/2
    first_second_product = (pair_mass_squared-first_mass**2-second_mass**2)/2
    gram_matrix = sp.Matrix([
        [vector_mass**2, vector_first_product, vector_second_product],
        [vector_first_product, first_mass**2, first_second_product],
        [vector_second_product, first_second_product, second_mass**2]])
    # The Levi-Civita contraction identity gives -J.J = det Gram(k,p1,p2).
    # J.k=J.p1=J.p2=0 by antisymmetry, so the longitudinal polarization
    # term and both odd-even contractions vanish without angular integration.
    odd_contraction = sp.factor(gram_matrix.det())
    even_contractions = sp.Matrix(2, 2, lambda row, column:
        sp.factor(-gram_matrix[row+1,column+1] +
                  gram_matrix[0,row+1]*gram_matrix[0,column+1]/vector_mass**2))
    squared_amplitude = (odd_contraction*sp.Abs(odd_coefficient)**2 +
        even_contractions[0,0]*sp.Abs(first_coefficient)**2 +
        even_contractions[1,1]*sp.Abs(second_coefficient)**2 +
        2*even_contractions[0,1]*sp.re(first_coefficient*sp.conjugate(second_coefficient)))
    kallen = lambda first, second, third: sp.expand(
        first**2+second**2+third**2-2*first*second-2*first*third-2*second*third)
    first_kallen = kallen(alp_mass**2, vector_mass**2, pair_mass_squared)
    second_kallen = kallen(pair_mass_squared, first_mass**2, second_mass**2)
    # In the P1 P2 rest frame, t=(k+p2)^2 and z is their relative cosine.
    pair_energy_vector = (alp_mass**2-vector_mass**2-pair_mass_squared)/(2*sp.sqrt(pair_mass_squared))
    pair_energy_second = (pair_mass_squared+second_mass**2-first_mass**2)/(2*sp.sqrt(pair_mass_squared))
    momentum_product = sp.sqrt(first_kallen)*sp.sqrt(second_kallen)/(4*pair_mass_squared)
    second_invariant_of_cosine = (vector_mass**2+second_mass**2 +
        2*pair_energy_vector*pair_energy_second-2*momentum_product*cosine)
    lower_boundary = sp.simplify(second_invariant_of_cosine.subs(cosine, 1))
    upper_boundary = sp.simplify(second_invariant_of_cosine.subs(cosine, -1))
    jacobian = sp.simplify(-sp.diff(second_invariant_of_cosine, cosine))
    # dPhi3 = ds/(2pi) dPhi2(a;V,Q) dPhi2(Q;P1,P2).
    # Overall orientation and the second decay azimuth have been integrated.
    phase_space_cosine = (sp.sqrt(first_kallen)/(8*sp.pi*alp_mass**2) *
        sp.sqrt(second_kallen)/(16*sp.pi*pair_mass_squared)/(2*sp.pi))
    phase_space_mandelstam = sp.simplify(phase_space_cosine/jacobian)
    width_density = sp.simplify(phase_space_mandelstam/(2*alp_mass*identical_factor))
    # The strangeness charge is fixed in the (u,d,s) basis. Its double
    # commutator selects precisely us,su,ds,sd entries; no isospin limit needed.
    vector_scalar_source = sp.MatrixSymbol('c_q^{VS}', 3, 3)
    strangeness = sp.diag(0, 0, -1)
    first_commutator = strangeness*vector_scalar_source-vector_scalar_source*strangeness
    restricted_source = sp.simplify(strangeness*first_commutator-first_commutator*strangeness)
    return locals()


def show(label, expression):
    """Print each analytical step using SymPy's LaTeX renderer."""
    print(f'{label}:\n{sp.latex(expression)}\n')


def main():
    """Display the reproducible invariant derivation and integration formula."""
    result = derive_width()
    print('Natural units, metric (+---), on-shell external vector; no initial spin average. '
          'S_PP=2! for identical pseudoscalars and 1 otherwise. Coefficients include '
          'their full complex phases and all coherent internal exchanges.')
    for label, name in [
        ('Mandelstam elimination u =', 'third_replacement'),
        ('Gram matrix ordered as k,p1,p2', 'gram_matrix'),
        ('Odd polarization contraction', 'odd_contraction'),
        ('Even polarization contractions ordered as p1,p2', 'even_contractions'),
        ('Polarization-summed squared amplitude', 'squared_amplitude'),
        ('t(s,z), with z the k,p2 cosine in the PP rest frame', 'second_invariant_of_cosine'),
        ('t minimum', 'lower_boundary'), ('t maximum', 'upper_boundary'),
        ('Absolute dt/dz', 'jacobian'), ('dPhi3/(ds dz)', 'phase_space_cosine'),
        ('dPhi3/(ds dt)', 'phase_space_mandelstam'),
        ('Multiplier of squared amplitude in dGamma/(ds dt)', 'width_density'),
        ('c_q^VS restricted to |Delta S|=1', 'restricted_source')]:
        show(label, result[name])
    print('Integrate s from (m_1+m_2)^2 to (m_a-m_V)^2 and t from t minimum '
          'to t maximum. The width vanishes below threshold. Odd-even interference '
          'is zero pointwise by antisymmetry. Retain F_1 F_2* interference. '
          'Insert the restricted source into the two even flavour traces of '
          'alp-vpp-amplitude.md; F_odd is unchanged. A pole on the integration '
          'domain requires an explicitly specified finite line shape for a finite width.')


if __name__ == '__main__':
    main()
