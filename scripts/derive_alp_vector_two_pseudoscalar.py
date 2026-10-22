"""Derive VVP polynomials and audit a -> V P P fit null directions.

Invocation: .venv/bin/python scripts/derive_alp_vector_two_pseudoscalar.py
Alda pseudoscalars: Phi=sum_P lambda_P P, tr(lambda_i lambda_j)=2 delta_ij,
u=exp(i Phi/(2 f_pi)), f_pi=0.093 GeV. Internal RChT formulas use
Phi_R=Phi/sqrt(2), lambda_R=lambda/sqrt(2), F=f_pi. Each meson
insertion therefore carries 1/sqrt(2) relative to a formula in lambda_R.
Vector generators retain tr(T_i T_j)=delta_ij.
Natural units, GeV; metric (+---), epsilon^{0123}=+1. The imported
contractions are the reference tensor calculation; the canonical amplitude
normalization is fixed by memory/rcht/odd-parity/alp-vv-amplitude.md.
This script derives the strong VVP coefficient, its three vector exchanges,
the generic vector-source scalar exchanges after a flavour rotation, and
the fit-null obstruction for electromagnetic conversion. Extra vector-mass
spurion operators and independent VPPP operators are outside its truncation.
"""
from pathlib import Path
import json
import sympy as sp
from derive_radiative_tff_amplitude import contractions, q1, q2, e1, e2, reduce, m, s1, s2

ODD_COEFFICIENTS = sp.symbols('c_1:8 d_1:5', real=True)


def derive():
    """Contract both vector assignments before extracting coefficient rows."""
    first = contractions(q1, q2, e1, e2)
    second = contractions(q2, q1, e2, e1)
    derivative_trace, mass_trace = sp.symbols('D P', real=True)
    vector_scale = sp.Symbol('M_V', positive=True)
    c = ODD_COEFFICIENTS[:7]
    d = ODD_COEFFICIENTS[7:]
    # Native contractions include their Fourier factors. Their sum yields
    # -2 m^2, 2(m^2-s1-s2), and 8 for O1, O3, and O2, respectively.
    # The factor -2 converts their reference normalization to the canonical
    # source-normalized numerator recorded in alp-vv-amplitude.md.
    strong = reduce(-2*(d[0]*derivative_trace*(first[0]+second[0])
                    +d[1]*mass_trace*(first[2]+second[2])
                    +d[2]*derivative_trace*(first[4]+second[4])
                    +d[3]*derivative_trace*(first[6]+second[6])))
    # VJP derivative and mass source normalizations are the same as the
    # canonical no-mixing aVV record. First is one ordered vector assignment.
    ordered = reduce(-2*(derivative_trace*sum(c[index]*first[index]
                    for index in (0,1,4,5,6))+c[2]*mass_trace*first[2])/vector_scale)
    exact_nulls = sp.zeros(11,4)
    exact_nulls[0,0] = exact_nulls[4,0] = 1
    for index in (1,4,5):
        exact_nulls[index,1] = 1
    exact_nulls[6,2] = exact_nulls[10,3] = 1
    strong_row = sp.Matrix([[sp.diff(strong, coefficient)
                            for coefficient in ODD_COEFFICIENTS]])
    ordered_row = sp.Matrix([[sp.diff(ordered, coefficient)
                             for coefficient in ODD_COEFFICIENTS]])
    return {'first':sp.Matrix(first), 'second':sp.Matrix(second),
            'strong':strong,'ordered':ordered,'nulls':exact_nulls,
            'strong_projection':(strong_row*exact_nulls).applyfunc(reduce),
            'conversion_projection':(ordered_row*exact_nulls).applyfunc(reduce),
            'D':derivative_trace,'P':mass_trace,'M_V':vector_scale}


def fitted_null_verification(exact_nulls):
    """Compare the exact four-dimensional space with the actual fit archive."""
    from summarize_odd_preliminary_fit import main as reconstruct_central_fit
    _, fitted_columns = reconstruct_central_fit(write_output=False)
    fitted_nulls = sp.Matrix(fitted_columns)
    annihilator = sp.Matrix.hstack(*exact_nulls.T.nullspace()).T
    discrepancy = annihilator*fitted_nulls
    return max(abs(float(entry)) for entry in discrepancy)


def show(label, expression):
    print(f'{label}:\n{sp.latex(expression)}\n')


def main():
    result = derive()
    print('The coefficient traces are D = Tr({T_i,T_j} U)/2 and '
          'P = Tr({T_i,T_j} Y)/2. For an ALP, U=2 c_q^A/f_a and '
          'Y=-2{c_q^P,chi}/f_a. For a meson, U=-lambda/f_pi '
          'and Y=+{lambda,chi}/f_pi. chi=2 B_0 m_q. '
          'Outgoing charged generators are Hermitian conjugated. Alda normalization: tr(lambda_i lambda_j)=2 delta_ij; f_pi=0.093 GeV. The retained vector basis has tr(T_i T_j)=delta_ij; RChT pseudoscalar generators equal lambda/sqrt(2).')
    show('Ordered tensor contractions O1,...,O7',result['first'])
    show('VVP numerator including both vector assignments',result['strong'])
    show('Ordered VJP numerator',result['ordered'])
    show('Exact four original-coefficient null columns',result['nulls'])
    show('VVP coefficient row times all fit nulls',result['strong_projection'])
    show('Ordered VJP row times all fit nulls',result['conversion_projection'])
    print('The last result has a nonzero first entry for s_2 != 0: '
          'conversion into a massive vector requires an unmeasured '
          'original-operator coefficient. At e=0 the VJP contribution is absent.')
    print('Numerical fit-space discrepancy:',fitted_null_verification(result['nulls']))
    print('Strong exchange topologies: a -> V V* -> V P1 P2; '
          'a -> P1 V* -> P1 V P2; a -> P2 V* -> P2 V P1. '
          'Every allowed flavor generator of the retained vector nonet is summed. '
          'The VPP and aVP tensor currents are proportional to p1 wedge p2 '
          'and p_a wedge p_i. Each contains the internal momentum as a wedge '
          'factor; therefore the above d4 zero applies also to these lines. '
          'Broadness does not decide which internal vectors occur.')
    print('A general c_q^VS additionally requires the even direct aVPP '
          'and pseudoscalar-exchange graphs, derived below by the allowed '
          'flavour rotation. '
          'Independent VPPP odd operators are excluded by the requested '
          'eleven-operator truncation, rather than fixed by the fit.')






def derive_tensor_exchange():
    """Invert the free tensor Hessian on the actual VPP momentum current.

    The source term is J V_{mu nu} C^{mu nu}, C=p1 wedge p2,
    q=p1+p2. Both antisymmetric components are included in this notation.
    The returned response is the stationary classical field; Fourier and
    vertex factors must be included in J before constructing i M.
    """
    from itertools import combinations
    q0, qz, vector_mass, source = sp.symbols('q_0 q_z M_R J', real=True)
    momentum = sp.Matrix([q0,0,0,qz])
    pseudoscalar_momentum = sp.Matrix(sp.symbols('p_1_0:4',real=True))
    metric = sp.diag(1,-1,-1,-1)
    pairs = list(combinations(range(4),2))
    components = sp.Matrix(sp.symbols('v_0:6'))
    tensor = sp.zeros(4)
    for component,(mu,nu) in zip(components,pairs):
        tensor[mu,nu] = component
        tensor[nu,mu] = -component
    divergence = tensor.T*metric*momentum
    quadratic = -(divergence.T*metric*divergence)[0]/2
    quadratic += vector_mass**2*sum(metric[mu,mu]*metric[nu,nu]*tensor[mu,nu]**2
                                    for mu,nu in pairs)/2
    hessian = sp.hessian(quadratic,components)
    second = momentum-pseudoscalar_momentum
    current = pseudoscalar_momentum*second.T-second*pseudoscalar_momentum.T
    linear = 2*source*sum(metric[mu,mu]*metric[nu,nu]*tensor[mu,nu]*current[mu,nu]
                         for mu,nu in pairs)
    response = (-hessian.inv()*sp.Matrix([sp.diff(linear,v) for v in components])).applyfunc(sp.factor)
    expected = sp.Matrix([2*source*current[mu,nu]/(q0**2-qz**2-vector_mass**2)
                          for mu,nu in pairs])
    assert (response-expected).applyfunc(sp.simplify)==sp.zeros(6,1)
    return hessian,response,expected




def derive_odd_exchange_amplitude():
    """Assemble all three strong vector exchanges with explicit flavour traces.

    A real Hermitian orthonormal nonet basis is summed. Charged physical
    states are complex linear combinations of this basis; consequently the
    outgoing pseudoscalar generators are conjugated, while each internal
    real-basis generator occurs unchanged at its two endpoints.
    """
    result = derive()
    external_vector, internal_vector = sp.MatrixSymbol('T_V',3,3),sp.MatrixSymbol('T_R',3,3)
    first_generator,second_generator = sp.MatrixSymbol('lambda_1',3,3),sp.MatrixSymbol('lambda_2',3,3)
    axial_coupling,mass_coupling,quark_mass = [sp.MatrixSymbol(name,3,3)
                                            for name in ('c_q^A','c_q^P','chi')]
    decay_constant,alp_scale,even_coupling = sp.symbols('f_pi f_a G_V',positive=True)
    external_mass,internal_mass,alp_mass,first_mass,second_mass = sp.symbols('M_ext M_R m_a m_1 m_2',positive=True)
    pair_squared,first_squared,second_squared = sp.symbols('s t u',real=True)
    symmetric_vectors=external_vector*internal_vector+internal_vector*external_vector

    def numerator(derivative_source,mass_source,source_mass,internal_squared):
        return result['strong'].subs({result['D']:sp.Trace(symmetric_vectors*derivative_source)/2,
                                     result['P']:sp.Trace(symmetric_vectors*mass_source)/2,
                                     m:source_mass,s1:external_mass**2,s2:internal_squared})

    alp_numerator=numerator(2*axial_coupling/alp_scale,
                           -2*(mass_coupling*quark_mass+quark_mass*mass_coupling)/alp_scale,
                           alp_mass,pair_squared)
    first_outgoing,second_outgoing=sp.Adjoint(first_generator),sp.Adjoint(second_generator)
    first_numerator=numerator(-first_outgoing/decay_constant,
                             (first_outgoing*quark_mass+quark_mass*first_outgoing)/decay_constant,
                             first_mass,second_squared)
    second_numerator=numerator(-second_outgoing/decay_constant,
                              (second_outgoing*quark_mass+quark_mass*second_outgoing)/decay_constant,
                              second_mass,first_squared)
    pair_trace=sp.Trace(internal_vector*(first_outgoing*second_outgoing-second_outgoing*first_outgoing))
    first_trace=sp.Trace(internal_vector*(first_outgoing*axial_coupling-axial_coupling*first_outgoing))
    second_trace=sp.Trace(internal_vector*(second_outgoing*axial_coupling-axial_coupling*second_outgoing))
    pair_exchange=-even_coupling/sp.sqrt(2)*alp_numerator*pair_trace/(decay_constant**2*external_mass*(pair_squared-internal_mass**2))
    first_exchange=sp.sqrt(2)*even_coupling*second_numerator*first_trace/(decay_constant*alp_scale*external_mass*(first_squared-internal_mass**2))
    second_exchange=-sp.sqrt(2)*even_coupling*first_numerator*second_trace/(decay_constant*alp_scale*external_mass*(second_squared-internal_mass**2))
    coefficient=pair_exchange+first_exchange+second_exchange
    null_projection=sp.Matrix([[sp.diff(coefficient,value) for value in ODD_COEFFICIENTS]])*result['nulls']
    assert null_projection==sp.zeros(1,4)
    return coefficient,null_projection




def derive_diagonal_vector_source_identity():
    """Derive the identities for the conserved diagonal vector-source case.

    Diagonality is an explicit scope condition here. It is not imposed on
    the general flavour sources of the odd amplitude above.
    """
    vector_source=sp.diag(*sp.symbols('c_u^VS c_d^VS c_s^VS',real=True))
    quark_mass=sp.diag(*sp.symbols('chi_u chi_d chi_s',real=True))
    neutral_vector=sp.diag(*sp.symbols('T_u T_d T_s',real=True))
    neutral_commutator=vector_source*neutral_vector-neutral_vector*vector_source
    mass_commutator=vector_source*quark_mass-quark_mass*vector_source
    # For a charged complex meson with vector charge h, its aPP vertex
    # contains h*(p_out^2-p_internal^2)/f_a. Multiplication by its
    # propagator cancels the inverse propagator on the on-shell external
    # leg. This is the scalar-current Ward identity needed to cancel the
    # aVPP contact obtained by covariantizing the two VPP derivatives.
    internal_squared,external_mass,vector_charge,alp_scale=sp.symbols('r_squared m_P h f_a',real=True)
    source_vertex=vector_charge*(external_mass**2-internal_squared)/alp_scale
    exchange_factor=sp.cancel(source_vertex/(internal_squared-external_mass**2))
    assert neutral_commutator==sp.zeros(3)
    assert mass_commutator==sp.zeros(3)
    assert exchange_factor==-vector_charge/alp_scale
    return neutral_commutator,mass_commutator,exchange_factor




def derive_general_vector_source_vertices():
    """Rotate away cVS and differentiate the resulting scalar mass source.

    Alda:2026nrd permits cV -> cV+delta, cS -> cS+delta. Starting
    with cV=cVS,cS=0, choose delta=-cVS. Both displayed vertices follow
    from the retained mass terms; their exchange reproduces the contact
    and derivative-source graphs before that invertible field rotation.
    Extra vector-mass spurion operators are excluded by the specified
    even-sector truncation.
    """
    first_field,second_field,alp_field=sp.symbols('P_1 P_2 a')
    first_generator,second_generator=[sp.MatrixSymbol(name,3,3)
                                      for name in ('lambda_1','lambda_2')]
    vector_source,quark_mass=[sp.MatrixSymbol(name,3,3)
                              for name in ('c_q^VS','chi')]
    alp_scale=sp.Symbol('f_a',real=True)
    commutator=vector_source*quark_mass-quark_mass*vector_source
    shifted_mass=quark_mass-sp.I*alp_field*commutator/alp_scale
    pseudoscalar=first_field*sp.Adjoint(first_generator)+second_field*sp.Adjoint(second_generator)
    # Alda L_chi at Phi^2 is -Tr(Phi^2 chi)/4.
    pseudoscalar_lagrangian=-pseudoscalar*pseudoscalar*shifted_mass/4
    def field_coefficient(expression):
        # The expression is quadratic in the two meson fields and linear
        # in a. Its mixed finite difference extracts exactly their product.
        # This avoids unsupported MatrixExpr scalar derivative operations.
        from itertools import product
        return sp.MatAdd(*( (-1)**(3-sum(values))*expression.subs(dict(zip(
            (first_field,second_field,alp_field),values)))
            for values in product((0,1),repeat=3))).doit()
    pseudoscalar_vertex=field_coefficient(pseudoscalar_lagrangian)
    return sp.Trace(pseudoscalar_vertex)


def derive_general_vector_source_exchange():
    """Contract scalar-source and G_V vertices with each scalar propagator.

    i M=(i A) i/(r^2-m_R^2) (i B). The VPP coefficient B uses the
    outgoing external tensor wavefunction -i(k wedge epsilon)/M_ext.
    Scalar products epsilon.p1 and epsilon.p2 remain explicit.
    """
    first,second,internal=[sp.MatrixSymbol(name,3,3)
                           for name in ('lambda_1','lambda_2','lambda_R')]
    vector,source,mass=[sp.MatrixSymbol(name,3,3)
                        for name in ('T_V','c_q^VS','chi')]
    alp_scale,decay_constant,even_coupling,external_mass,internal_mass=sp.symbols('f_a f_pi G_V M_ext m_R',positive=True)
    first_squared,second_squared=sp.symbols('t u',real=True)
    first_polarization,second_polarization=sp.symbols('epsilon_dot_p1 epsilon_dot_p2',real=True)
    commutator=source*mass-mass*source
    contributions=[]
    for first_generator,second_generator,internal_squared,polarization_product in (
            (sp.Adjoint(first),sp.Adjoint(second),first_squared,second_polarization),
            (sp.Adjoint(second),sp.Adjoint(first),second_squared,first_polarization)):
        scalar_vertex=sp.I*sp.Trace((first_generator*internal+internal*first_generator)*commutator)/(4*alp_scale)
        vector_vertex=even_coupling/sp.sqrt(2)*external_mass*polarization_product*sp.Trace(vector*(internal*second_generator-second_generator*internal))/decay_constant**2
        diagram=(sp.I*scalar_vertex)*(sp.I/(internal_squared-internal_mass**2))*(sp.I*vector_vertex)/sp.I
        contributions.append(diagram)
    return sp.Add(*contributions)


def derive_exchange_lorentz_contractions():
    """Reduce the three actual wedge-response epsilon contractions."""
    polarization,vector_momentum,first,second=[sp.Matrix(sp.symbols(name+'_0:4'))
                                              for name in ('eps','k','p1','p2')]
    reference=sp.Matrix.hstack(polarization,first,second,vector_momentum).det()
    pair=sp.Matrix.hstack(polarization,second,vector_momentum,first+second).det()
    first_crossed=sp.Matrix.hstack(polarization,first,vector_momentum,vector_momentum+second).det()
    second_crossed=sp.Matrix.hstack(polarization,second,vector_momentum,vector_momentum+first).det()
    residuals=sp.Matrix([sp.expand(pair-reference),sp.expand(first_crossed+reference),
                         sp.expand(second_crossed-reference)])
    assert residuals==sp.zeros(3,1)
    return residuals


def derive_direct_even_vertex():
    """Extract the original-frame G_V aPP contact from both vielbeins.

    Noncommuting symbols denote flavour matrices inside the final trace;
    Lorentz-index labels denote components, never summed in this function.
    The full Lagrangian contracts this expression with V_mu_nu and traces.
    """
    first_field,second_field,alp_field=sp.symbols('P_1 P_2 a')
    first,second,axial,vector=sp.symbols('lambda_1 lambda_2 c_q^A c_q^VS',commutative=False)
    first_mu,first_nu,second_mu,second_nu,alp_mu,alp_nu=sp.symbols('p1_mu p1_nu p2_mu p2_nu pa_mu pa_nu')
    decay_constant,alp_scale,even_coupling=sp.symbols('f_pi f_a G_V',positive=True)
    meson=first_field*first+second_field*second
    commutator=meson*vector-vector*meson
    vielbein_mu=-(first_field*first_mu*first+second_field*second_mu*second)/decay_constant
    vielbein_nu=-(first_field*first_nu*first+second_field*second_nu*second)/decay_constant
    vielbein_mu+=2*alp_field*alp_mu*axial/alp_scale-sp.I*alp_field*alp_mu*commutator/(decay_constant*alp_scale)
    vielbein_nu+=2*alp_field*alp_nu*axial/alp_scale-sp.I*alp_field*alp_nu*commutator/(decay_constant*alp_scale)
    lagrangian=sp.I*even_coupling*(vielbein_mu*vielbein_nu-vielbein_nu*vielbein_mu)/sp.sqrt(2)
    coefficient=sp.expand(lagrangian).coeff(first_field).coeff(second_field).coeff(alp_field)
    assert not coefficient.has(axial)
    return coefficient




if __name__=='__main__':
    main()
    hessian, response, _ = derive_tensor_exchange()
    show('Free antisymmetric tensor Hessian for q=(q0,0,0,qz)',hessian)
    show('VPP source response, six independent tensor components',response)
    print('For L_source=J V_{mu nu}(p1^mu p2^nu-p2^mu p1^nu), '
          'the stationary tensor is 2J(p1 wedge p2)/(q²-M_R²). '
          'For the Fourier-transformed G_V vertex, '
          'J=-i G_V Tr(T_R[lambda1^dagger,lambda2^dagger])/(2 sqrt(2) f_pi²). '
          'This response includes the local part of the tensor propagator.')
    coefficient,null_projection=derive_odd_exchange_amplitude()
    print('Sum the following coefficient over every real Hermitian normalized '
          'generator T_R of the vector nonet, using its own pole mass M_R. '
          'Multiply the sum by epsilon(epsilon_V,p1,p2,k). '
          's=(p1+p2)^2, t=(k+p2)^2, u=(k+p1)^2; k is the outgoing '
          'vector momentum. chi=2 B_0 m_q. Every source and trace is displayed.')
    show('Three strong vector exchanges: coefficient per internal nonet generator',coefficient)
    show('Complete strong odd exchange coefficient times all four fitted nulls',null_projection)
    neutral_commutator,mass_commutator,exchange_factor=derive_diagonal_vector_source_identity()
    show('For diagonal cVS and final omega/phi: [cVS,T_V]',neutral_commutator)
    show('For diagonal cVS: [cVS,chi]',mass_commutator)
    show('aPP vertex times pseudoscalar propagator on its external mass shell',exchange_factor)
    print('For diagonal cVS, the covariant-vector kinetic aVV flavour trace '
          'vanishes. The remaining even terms cancel by conservation of '
          'the diagonal vector current: covariantizing each derivative in '
          'VPP produces h/f_a times that derivative, while the corresponding '
          'pseudoscalar-exchange contribution gives -h/f_a. '
          'This statement assumes cS=0 in the no-mixing basis and [cVS,chi]=0. '
          'Off-diagonal cVS is treated by the general scalar-source '
          'exchange expression below, subject to its stated mass-operator scope.')
    pseudoscalar_vertex=derive_general_vector_source_vertices()
    show('General cVS scalar aPP Lagrangian coefficient after the allowed vector rotation',pseudoscalar_vertex)
    show('Two generic cVS pseudoscalar exchanges per real pseudoscalar generator lambda_R',derive_general_vector_source_exchange())
    show('Residuals of all three explicit epsilon-contraction reductions',derive_exchange_lorentz_contractions())
    show('Original-frame direct G_V aPP contact inside Tr(T_V ...), per antisymmetric Lorentz pair',derive_direct_even_vertex())
    print('Sum the last expression over the canonically normalized retained '
          'pseudoscalar mass eigenstates, using their individual masses m_R. '
          'epsilon_dot_p1 and epsilon_dot_p2 denote epsilon_V dot p1 and '
          'epsilon_V dot p2. The expression excludes the rotated vector-mass '
          'spurion vertex unless its parent operator is explicitly retained.')
