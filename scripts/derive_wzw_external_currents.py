"""Expand the two-current part of Bai:2025fvl Eq. (A11).

Invocation: .venv/bin/python scripts/derive_wzw_external_currents.py

Input: Bai:2025fvl Eq. (A11). Its left/right naming is mapped to Alda
by D U=dU-i r U+i U l; at U=1, u_mu=r-l=2 hat cA da/fa.
The pure meson functional has no ALP; terms with three or more external
currents cannot contribute; the source-only counterterm has no pions.
The remaining four two-current groups are expanded explicitly below.
Natural units; e=f_pi=f_a=1 during coefficient extraction. No mass relation.
"""
import sympy as sp

plus, minus, photon, alp = sp.symbols('pi_plus pi_minus A_gamma a')
charge, axial, generator_plus, generator_minus = sp.symbols(
    r'Q \widehat{c}_q^A lambda_+^dagger lambda_-^dagger', commutative=False)
fields = (plus, minus, photon, alp)


def truncate(expression):
    """Keep only powers that can contribute to one of each external field."""
    expanded=sp.expand(expression)
    return sp.Add(*(term for term in sp.Add.make_args(expanded)
        if all(term.as_powers_dict().get(field,0)<=1 for field in fields)))


def multiply(*forms):
    """Exterior product of matrix-valued forms preserving flavor order.

    SymPy has no matrix-valued exterior product. Its noncommutative algebra
    handles coefficients; this helper only antisymmetrizes coordinate indices.
    """
    result = {(): sp.S.One}
    for form in forms:
        updated = {}
        for indices, coefficient in result.items():
            for other_indices, other_coefficient in form.items():
                combined = indices+other_indices
                if len(set(combined)) != len(combined):
                    continue
                sign = (-1)**sum(combined[i]>combined[j]
                    for i in range(len(combined)) for j in range(i+1,len(combined)))
                ordered = tuple(sorted(combined))
                updated[ordered] = updated.get(ordered,0)+sign*coefficient*other_coefficient
        result = {indices:truncate(value) for indices,value in updated.items()}
    return result


def add(*forms):
    """Add forms with the same ordered coordinate indices."""
    result = {}
    for form in forms:
        for indices,value in form.items():
            result[indices]=result.get(indices,0)+value
    return result


def scale(coefficient, form):
    """Multiply every form coefficient from the left."""
    return {indices:coefficient*value for indices,value in form.items()}


def derivative(expression):
    """Differentiate plane waves with independent vector directions.

    epsilon_gamma=e0, k=e1, p_plus=e2, p_minus=e3,
    p_a=k+p_plus+p_minus; all fields outgoing except the incoming ALP.
    This extracts the epsilon tensor coefficient without fixing invariants.
    """
    momenta=((0,0,1,0),(0,0,0,1),(0,1,0,0),(0,-1,-1,-1))
    return {(index,):sp.expand(sp.I*sum(momentum[index]*field*sp.diff(expression,field)
        for field,momentum in zip(fields,momenta))) for index in range(4)}


def exterior_derivative(form):
    """Differentiate the coefficients of a coordinate form."""
    return add(*(multiply(derivative(value),{indices:1}) for indices,value in form.items()))


def cyclic_trace(expression):
    """Collect cyclic words; SymPy Trace does not collect cyclic rotations."""
    if sp.expand(expression)==0:
        return sp.S.Zero
    collected={}
    for term in sp.Add.make_args(sp.expand(expression)):
        coefficient,word=term.args_cnc()
        word=tuple(item for factor in word for item in
            ([factor.base]*int(factor.exp) if isinstance(factor,sp.Pow) and factor.exp.is_Integer else [factor]))
        representative=min((word[index:]+word[:index] for index in range(len(word))),
            key=lambda entry:tuple(str(item) for item in entry))
        collected[representative]=collected.get(representative,0)+sp.Mul(*coefficient)
    matrices={symbol:sp.MatrixSymbol(str(symbol),3,3)
        for symbol in expression.atoms(sp.Symbol) if not symbol.is_commutative}
    return sp.Add(*(sp.simplify(coefficient)*sp.Trace(sp.Mul(*(matrices[item] for item in word)))
        for word,coefficient in collected.items() if sp.simplify(coefficient)!=0))


def derive_terms(photon_polarization=(1,0,0,0)):
    """Extract one photon, one ALP and both pion fields from Eq. (A11)."""
    phi=plus*generator_plus+minus*generator_minus
    unitary={():1+sp.I*phi-phi**2/2}
    adjoint={():1-sp.I*phi-phi**2/2}
    alpha=multiply(exterior_derivative(unitary),adjoint)
    beta=multiply(adjoint,exterior_derivative(unitary))
    vector={(index,):-photon*charge*value for index,value in enumerate(photon_polarization)}
    axial_current=scale(axial,derivative(alp))
    left=add(vector,axial_current)
    right=add(vector,scale(-1,axial_current))
    derivative_left=exterior_derivative(left)
    derivative_right=exterior_derivative(right)
    terms=[
        scale(-sp.I/2,add(multiply(left,alpha,left,alpha),scale(-1,multiply(right,beta,right,beta)))),
        scale(sp.I,add(multiply(left,unitary,right,adjoint,alpha,alpha),
            scale(-1,multiply(right,adjoint,left,unitary,beta,beta)))),
        scale(sp.I,add(multiply(derivative_right,exterior_derivative(adjoint),left,unitary),
            scale(-1,multiply(derivative_left,exterior_derivative(unitary),right,adjoint)))),
        scale(sp.I,add(multiply(derivative_left,left,alpha),multiply(left,derivative_left,alpha),
            multiply(derivative_right,right,beta),multiply(right,derivative_right,beta))),
    ]
    coefficients=[]
    for term in terms:
        expression=sp.expand(term.get((0,1,2,3),0))
        for field in fields:
            expression=expression.coeff(field,1)
        coefficients.append(sp.expand(expression))
    return coefficients


def compare_pion_normalizations():
    """Derive decay-constant conversion from kinetic and exponent matching.

    Bai Eq. (A10) has f_Bai^2/8 and exponent sqrt(2)*i*pi^a*tau^a/f_Bai.
    The project has f_pi^2/4 and exponent i*Phi/f_pi. Identical canonically
    normalized pion fields identify pi^a*tau^a with the project matrix.
    """
    bai_decay_constant,project_decay_constant=sp.symbols('f_pi_Bai f_pi',positive=True)
    conversion=sp.solve(sp.Eq(bai_decay_constant**2/8,project_decay_constant**2/4),bai_decay_constant)[0]
    print('Decay-constant mapping from equality of kinetic prefactors:')
    print(sp.latex(sp.Eq(bai_decay_constant,conversion)))
    exponent_difference=sp.simplify(sp.sqrt(2)/bai_decay_constant-1/project_decay_constant)
    print('Difference of charged-pion exponent coefficients after conversion:')
    print(sp.latex(sp.simplify(exponent_difference.subs(bai_decay_constant,conversion))))
    return conversion


def derive_gamma_three_pion():
    """Derive the photon-linear cubic-pion coefficient from the same A11.

    At cubic order alpha=beta=i*dPhi/f_pi. Matrix-valued exterior
    multiplication generates all six flavor orderings. The numerical
    alternating trace 2 is the indexed gamma3pi record's normalization
    input, not a new definition of pion generators.
    """
    pion_decay_constant,electric_charge,number_colors=sp.symbols('f_pi e N_c',positive=True)
    gamma_field=sp.Symbol('A_gamma')
    first_pion_derivative,second_pion_derivative,neutral_pion_derivative=sp.symbols('d_pi_plus d_pi_minus d_pi0')
    neutral_generator=sp.Symbol('lambda_0',commutative=False)
    derivative_phi={(1,):first_pion_derivative*generator_minus,
        (2,):second_pion_derivative*generator_plus,(3,):neutral_pion_derivative*neutral_generator}
    alpha=scale(sp.I/pion_decay_constant,derivative_phi)
    photon_source={(0,):-electric_charge*gamma_field*charge}
    cubic_term=add(multiply(photon_source,alpha,alpha,alpha),multiply(photon_source,alpha,alpha,alpha))
    coefficient=sp.expand(cubic_term[(0,1,2,3)])
    for field in (gamma_field,first_pion_derivative,second_pion_derivative,neutral_pion_derivative):
        coefficient=coefficient.coeff(field,1)
    # Extract the common coefficient; flavor order is retained in the polynomial.
    common_coefficient=2*(-electric_charge)*(sp.I/pion_decay_constant)**3
    alternating_words=sp.expand(coefficient/common_coefficient)
    print('Ordered cubic-pion flavor polynomial from A11 (inside the trace):')
    print(sp.latex(alternating_words))
    coefficient_in_project=sp.simplify(number_colors/(48*sp.pi**2)*common_coefficient*2)
    print('gamma3pi derivative coefficient using the indexed alternating trace 2:')
    print(sp.latex(coefficient_in_project))
    orientation=sp.Symbol('orientation',real=True)
    recorded_coefficient=orientation*sp.I*number_colors*electric_charge/(12*sp.pi**2*pion_decay_constant**3)
    print('Ratio of A11 coefficient in chosen exterior orientation to indexed gamma3pi coefficient:')
    print(sp.latex(sp.simplify(coefficient_in_project/recorded_coefficient)))
    print('Both vertices carry the same WZW orientation when matched to project memory.')


def derive_neutral_pion_two_photons():
    """Expand A11 with two photons and one incoming source-basis pion.

    Incoming q1=e1, q2=e3 and p=-(q1+q2); photon polarizations e0,e2.
    Fourier derivatives are +i times incoming momentum. All vectors here
    have covariant indices, hence the exterior coefficient directly uses
    the canonical upper epsilon. The scalar R has covariant epsilon with
    contravariant vectors and is the same contraction after raising.
    """
    pion,first_photon,second_photon=sp.symbols('pi0 A1 A2')
    neutral_generator=sp.Symbol('lambda_0',commutative=False)
    pion_decay_constant,electric_charge,number_colors=sp.symbols('f_pi e N_c',positive=True)
    momenta={pion:(0,-1,0,-1),first_photon:(0,1,0,0),second_photon:(0,0,0,1)}
    def local_derivative(form):
        result={}
        for indices,expression in form.items():
            differentiated={(index,):sp.I*sum(momentum[index]*field*sp.diff(expression,field)
                for field,momentum in momenta.items()) for index in range(4)}
            result=add(result,multiply(differentiated,{indices:1}))
        return result
    unitary={():1+sp.I*pion*neutral_generator/pion_decay_constant}
    adjoint={():1-sp.I*pion*neutral_generator/pion_decay_constant}
    alpha=multiply(local_derivative(unitary),adjoint)
    beta=multiply(adjoint,local_derivative(unitary))
    vector={(0,):-electric_charge*first_photon*charge,(2,):-electric_charge*second_photon*charge}
    derivative_vector=local_derivative(vector)
    terms=[
        scale(-sp.I/2,add(multiply(vector,alpha,vector,alpha),scale(-1,multiply(vector,beta,vector,beta)))),
        scale(sp.I,add(multiply(vector,unitary,vector,adjoint,alpha,alpha),
            scale(-1,multiply(vector,adjoint,vector,unitary,beta,beta)))),
        scale(sp.I,add(multiply(derivative_vector,local_derivative(adjoint),vector,unitary),
            scale(-1,multiply(derivative_vector,local_derivative(unitary),vector,adjoint)))),
        scale(sp.I,add(multiply(derivative_vector,vector,alpha),multiply(vector,derivative_vector,alpha),
            multiply(derivative_vector,vector,beta),multiply(vector,derivative_vector,beta)))
    ]
    coefficients=[]
    for term in terms:
        coefficient=sp.expand(term.get((0,1,2,3),0))
        for field in (pion,first_photon,second_photon):
            coefficient=coefficient.coeff(field,1)
        coefficients.append(coefficient)
    print('Neutral-pion two-photon source vertex per A11 group, before Nc/(48 pi^2):')
    for coefficient in coefficients:
        print(sp.latex(cyclic_trace(coefficient)))
    print('Total neutral-pion source vertex:')
    print(sp.latex(number_colors/(48*sp.pi**2)*cyclic_trace(sum(coefficients))))
    # The plane-wave coordinate choice gives E=+1 and R=-1 explicitly.
    basis=sp.eye(4)
    tensor_e=sp.Matrix.vstack(basis[1,:],basis[0,:],basis[3,:],basis[2,:]).det()
    tensor_r=sp.Matrix.vstack(basis[0,:],basis[2,:],basis[1,:],basis[3,:]).det()
    print('Canonical E and R for the chosen covectors:')
    print(sp.latex(sp.Tuple(tensor_e,tensor_r)))
    flavor_trace,decay_constant_rcht=sp.symbols('F_Q F',real=True,nonzero=True)
    # Indexed field dictionary: Phi_A=sqrt(2) Phi_R, f_pi=F,
    # hence tr(Q^2 lambda0)=F_Q/sqrt(2).
    source_scalar_e=number_colors*flavor_trace/(4*sp.sqrt(2)*sp.pi**2*decay_constant_rcht)
    source_scalar_r=sp.simplify(source_scalar_e*tensor_e/tensor_r)
    # Indexed meson phase: fixed-phase physical meson is minus the source meson.
    fixed_phase_scalar_r=-source_scalar_r
    recorded_scalar=-sp.sqrt(2)*number_colors*flavor_trace/(8*sp.pi**2*decay_constant_rcht)
    print('A11 source scalar in E, then fixed-meson-phase scalar in R:')
    print(sp.latex(sp.Tuple(source_scalar_e,fixed_phase_scalar_r)))
    print('Recorded canonical WZW scalar in R:')
    print(sp.latex(recorded_scalar))
    matched_orientation=sp.simplify(recorded_scalar/fixed_phase_scalar_r)
    print('Common WZW orientation required by the signed neutral-pion comparison:')
    print(sp.latex(matched_orientation))
    return matched_orientation


def main():
    matched_orientation=derive_neutral_pion_two_photons()
    compare_pion_normalizations()
    derive_gamma_three_pion()
    terms=derive_terms()
    print("Natural units; O(e/f_a), quadratic pions; no quark-mass equality.")
    print("A_L=r and A_R=l: D U=dU-i r U+i U l, so u_mu at U=1 is r-l.")
    for number,term in enumerate(terms,1):
        print(f'Two-current group {number}:\n{sp.latex(cyclic_trace(term))}\n')
    number_colors,electric_charge,pion_decay_constant,alp_decay_constant=sp.symbols(
        'N_c e f_pi f_a',positive=True)
    prefactor=number_colors*electric_charge/(48*sp.pi**2*alp_decay_constant*pion_decay_constant**2)
    print('Amplitude coefficient multiplying epsilon_lower epsilon_gamma p_plus p_minus k:')
    print(sp.latex(matched_orientation*prefactor*cyclic_trace(sum(terms))))


if __name__=='__main__':
    main()
