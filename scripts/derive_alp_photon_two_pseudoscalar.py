"""Derive odd-resonance a -> gamma P1 P2 at leading electromagnetic order.

Run: .venv/bin/python scripts/derive_alp_photon_two_pseudoscalar.py
Tests: .venv/bin/python scripts/derive_alp_photon_two_pseudoscalar.py --test

Alda pseudoscalars: Phi=sum_P lambda_P P, u=exp(i Phi/(2 f_pi)),
tr(lambda_i lambda_j)=2 delta_ij, f_pi=0.093 GeV. RChT input
operators use Phi_R=Phi/sqrt(2), lambda_R=lambda/sqrt(2), F=f_pi;
all displayed pseudoscalar sources and traces use Alda generators.
Natural units, GeV; metric (+---), epsilon^0123=+1. Sources are already
in the no-mixing basis. Retain seven VJP, four VVP operators and G_V.
Electromagnetic compatibility requires [Q,c_q^A]=[Q,c_q^P]=0.
This derives their odd vector exchanges, not independent local WZW terms
or the even vector/scalar-source radiative amplitude. Tensor responses
include the local part of the antisymmetric-tensor propagator.
"""
import argparse
import unittest
import sympy as sp
from derive_alp_vector_two_pseudoscalar import derive as derive_two_vector
from derive_radiative_tff_amplitude import (
    contractions, q1, q2, e1, e2, m, s1, s2,
)


def show(label, expression):
    """Print a labelled equation through SymPy's LaTeX renderer."""
    print(label + ':\n' + sp.latex(expression) + '\n')


def derive():
    """Construct source vertices from the original operator contractions."""
    inherited = derive_two_vector()
    coefficients = sp.symbols('c_1:8 d_1:5', real=True)
    electric_charge, decay_constant, alp_scale, even_coupling = sp.symbols(
        'e f_pi f_a G_V', nonzero=True)
    normalization_mass = inherited['M_V']
    internal_mass, alp_mass, first_mass, second_mass = sp.symbols(
        'M_R m_a m_1 m_2', positive=True)
    pair_squared, first_squared, second_squared = sp.symbols('s t u', real=True)
    internal_vector, photon_response, charge = [sp.MatrixSymbol(name, 3, 3)
        for name in ('T_R', 'W', 'Q')]
    first_generator, second_generator = [sp.MatrixSymbol(name, 3, 3)
        for name in ('lambda_1', 'lambda_2')]
    first_outgoing, second_outgoing = map(sp.Adjoint,
        (first_generator, second_generator))
    axial, pseudoscalar, mass = [sp.MatrixSymbol(name, 3, 3)
        for name in ('c_q^A', 'c_q^P', 'chi')]

    def anticommutator(first, second):
        return first*second+second*first

    def commutator(first, second):
        return first*second-second*first

    def radiative_numerator(derivative_source, mass_source, source_mass,
                           vector_squared, outgoing_generator=None):
        # The external real photon occupies the first VVP tensor slot.
        converted = sp.sqrt(2)*electric_charge*inherited['strong'].subs({
            inherited['D']: sp.Trace(anticommutator(photon_response,
                internal_vector)*derivative_source)/2,
            inherited['P']: sp.Trace(anticommutator(photon_response,
                internal_vector)*mass_source)/2,
            m: source_mass, s1: 0, s2: vector_squared})
        # In VJP the resonance occupies the first slot, photon the second.
        direct = -2*electric_charge*inherited['ordered'].subs({
            inherited['D']: sp.Trace(anticommutator(internal_vector,
                charge)*derivative_source)/2,
            inherited['P']: sp.Trace(anticommutator(internal_vector,
                charge)*mass_source)/2,
            m: source_mass, s1: vector_squared, s2: 0})
        c4_piece = sp.S.Zero
        if outgoing_generator is not None:
            # f_-=-i e [Phi,Q] Fgamma/f_pi; chi_+=2chi.
            # O4=i epsilon Tr(V [f_-,chi_+])/M_V. Its scalar
            # source coefficient is obtained by multiplication below.
            lagrangian_coefficient = (sp.I*(-sp.I
                *electric_charge/decay_constant)*2
                *coefficients[3]/normalization_mass)
            # The canonical contraction is minus four; this follows
            # from the same epsilon/LSZ adapter as the O3 contraction.
            # It reproduces the -4 c4 F_K term in radiative-sm-amplitude.
            c4_piece = sp.simplify(-inherited['first'][2]
                *lagrangian_coefficient)*sp.Trace(
                internal_vector*commutator(commutator(
                    outgoing_generator, charge), mass))
        return sp.expand(converted+direct+c4_piece), converted, direct, c4_piece

    sources = (
        (2*axial/alp_scale,
         -2*anticommutator(pseudoscalar,mass)/alp_scale,
         alp_mass,pair_squared,None),
        (-second_outgoing/decay_constant,
         anticommutator(second_outgoing,mass)/decay_constant,
         second_mass,first_squared,second_outgoing),
        (-first_outgoing/decay_constant,
         anticommutator(first_outgoing,mass)/decay_constant,
         first_mass,second_squared,first_outgoing))
    radiative_pieces = [radiative_numerator(*source) for source in sources]
    traces = (
        sp.Trace(internal_vector*commutator(first_outgoing,second_outgoing)),
        sp.Trace(internal_vector*commutator(first_outgoing,axial)),
        sp.Trace(internal_vector*commutator(second_outgoing,axial)))
    prefactors = (
        -even_coupling/sp.sqrt(2)/(decay_constant**2
            *(pair_squared-internal_mass**2)),
        sp.sqrt(2)*even_coupling/(decay_constant*alp_scale
            *(first_squared-internal_mass**2)),
        -sp.sqrt(2)*even_coupling/(decay_constant*alp_scale
            *(second_squared-internal_mass**2)))
    exchanges = [prefactor*trace*pieces[0]
        for prefactor,trace,pieces in zip(prefactors,traces,radiative_pieces)]
    amplitude_coefficient = sp.Add(*exchanges)
    null_projection = (sp.Matrix([[sp.diff(amplitude_coefficient,coefficient)
        for coefficient in coefficients]])*inherited['nulls']).applyfunc(sp.expand)
    # Normalize rational scalar coefficients without asking SymPy to factor
    # abstract matrix sums inside Trace (unsupported by its polynomial API).
    for index, expression in enumerate(null_projection):
        traces = sorted(expression.atoms(sp.Trace), key=sp.default_sort_key)
        replacements = {trace:sp.Dummy() for trace in traces}
        null_projection[index] = sp.cancel(expression.xreplace(replacements)).xreplace(
            {value:key for key,value in replacements.items()})

    # Derive the signed photon response from its stationary quadratic form.
    residue, pole_mass, response, photon_field = sp.symbols(
        'R_i M_i r_i F_gamma', nonzero=True)
    quadratic = pole_mass**2*response**2/4 - sp.sqrt(2)*electric_charge*residue*response*photon_field/2
    response_solution = sp.solve(sp.diff(quadratic,response),response)[0]
    return locals()


class RadiativeAmplitudeTests(unittest.TestCase):
    """Gauge, normalization and charged-flavour checks independent of assembly."""

    @classmethod
    def setUpClass(cls):
        cls.result = derive()

    def test_alda_even_vertex_conversion(self):
        # Derive the pair and crossed G_V currents directly from the
        # vielbein product, independently of the assembled prefactors.
        result = self.result
        first_derivative, second_derivative, alp_derivative = sp.symbols(
            'dP_1 dP_2 da')
        first, second, axial = sp.symbols('lambda_1 lambda_2 c_q^A',
            commutative=False)
        vielbein_first = -first_derivative*first/result['decay_constant']
        vielbein_second = -second_derivative*second/result['decay_constant']
        vielbein_alp = 2*alp_derivative*axial/result['alp_scale']
        pair = sp.expand(vielbein_first*vielbein_second
                         -vielbein_second*vielbein_first)
        crossed = sp.expand(vielbein_first*vielbein_alp
                            -vielbein_alp*vielbein_first)
        self.assertEqual(pair.coeff(first_derivative).coeff(second_derivative),
            sp.expand((first*second-second*first)/result['decay_constant']**2))
        self.assertEqual(crossed.coeff(first_derivative).coeff(alp_derivative),
            sp.expand(-2*(first*axial-axial*first)
            /(result['decay_constant']*result['alp_scale'])))

    def test_original_fit_null_directions(self):
        self.assertEqual(self.result['null_projection'],sp.zeros(1,4))

    def test_photon_ward_identity_at_operator_level(self):
        self.assertEqual(contractions(q1,q2,e1,q2,False),[0]*7)
        self.assertEqual(contractions(q2,q1,q2,e1,False),[0]*7)

    def test_signed_photon_response(self):
        result = self.result
        self.assertEqual(sp.simplify(result['response_solution']
            -sp.sqrt(2)*result['electric_charge']*result['residue']
            *result['photon_field']/result['pole_mass']**2),0)

    def test_c4_charged_pion_and_kaon_mass_difference(self):
        charge = sp.diag(sp.Rational(2,3),-sp.Rational(1,3),-sp.Rational(1,3))
        light_mass,strange_mass = sp.symbols('chi_l chi_s')
        mass = sp.diag(light_mass,light_mass,strange_mass)
        for column,expected in ((1,0),(2,strange_mass-light_mass)):
            outgoing = sp.zeros(3)
            outgoing[0,column] = 1
            internal = outgoing.T
            first = outgoing*charge-charge*outgoing
            trace = sp.trace(internal*(first*mass-mass*first))
            self.assertEqual(sp.expand(trace+expected),0)

    def test_c4_ward_identity(self):
        # The otherwise omitted O4 is proportional to epsilon V F_gamma.
        polarization,photon_momentum,vector_momentum = [sp.Matrix(
            sp.symbols(name+'_0:4')) for name in ('eps','k','r')]
        self.assertEqual(sp.Matrix.hstack(vector_momentum,polarization,
            photon_momentum,photon_momentum).det(),0)

    def test_complete_charged_radiative_vertex_against_sm_derivation(self):
        from derive_radiative_sm_amplitude import derive as derive_sm
        import derive_radiative_sm_amplitude as sm
        result = self.result
        charge = sp.diag(sp.Rational(2,3),-sp.Rational(1,3),-sp.Rational(1,3))
        outgoing = sp.zeros(3)
        outgoing[0,2] = 1
        vector = outgoing.T
        chi_light,chi_strange,vector_source,pole_mass = sp.symbols(
            'chi_l chi_s F_V M_0', nonzero=True)
        mass = sp.diag(chi_light,chi_light,chi_strange)
        substitutions = {
            result['charge']:charge, result['internal_vector']:vector,
            result['second_generator']:outgoing.T,
            result['mass']:mass,
            result['photon_response']:vector_source*charge/pole_mass**2}
        pieces = [sp.expand(piece.subs(substitutions).doit())
            for piece in result['radiative_pieces'][1][1:]]
        charge_trace = sp.trace((vector*charge+charge*vector)*outgoing)
        mass_trace = sp.trace((vector*charge+charge*vector)
            *(outgoing*mass+mass*outgoing))
        fminus = outgoing*charge-charge*outgoing
        c4_trace = sp.trace(vector*(fminus*mass-mass*fminus))
        direct_numerator,converted_numerator = derive_sm()
        common_substitutions = {sm.fq:charge_trace,sm.fh:mass_trace,
            sm.fk:c4_trace,m:result['second_mass'],
            s1:result['first_squared'],s2:0}
        expected_direct = (2*result['electric_charge']
            /(result['decay_constant']*result['normalization_mass'])
            *direct_numerator.subs(common_substitutions))
        expected_converted = (-2*sp.sqrt(2)*result['electric_charge']*vector_source
            /(result['decay_constant']*pole_mass**2)
            *converted_numerator.subs(common_substitutions))
        self.assertEqual(sp.expand(pieces[1]+pieces[2]-expected_direct),0)
        self.assertEqual(sp.expand(pieces[0]-expected_converted),0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--test',action='store_true')
    arguments = parser.parse_args()
    if arguments.test:
        unittest.main(argv=['derive_alp_photon_two_pseudoscalar.py'])
        return
    result = derive()
    print('Leading e; no-mixing ALP sources. M_odd=epsilon(eps_gamma,p1,p2,k) '
          'times the sum over the real Hermitian normalized vector nonet. '
          's=(p1+p2)^2, t=(k+p2)^2, u=(k+p1)^2; k²=0. '
          'W=sum_i R_i T_i/M_i², with signed source residues R_i. '
          'Tree denominators carry +i0. f_pi=0.093 GeV; tr(lambda_i lambda_j)=2 delta_ij. RChT lambda_R=lambda/sqrt(2), F=f_pi. Each endpoint uses this conversion; vector generators retain tr(T_i T_j)=delta_ij. All masses have GeV units. '
          'Independent local WZW terms and even scalar/vector-source '
          'radiative terms are outside this odd-resonance derivation.')
    show('Quadratic mixing Lagrangian per tensor component',result['quadratic'])
    show('Stationary real-photon tensor response',result['response_solution'])
    for name,pieces in zip(('a;s','P2;t','P1;u'),result['radiative_pieces']):
        for label,piece in zip(('complete','VVP conversion','direct VJP','VJP c4'),pieces):
            show(name+' radiative numerator '+label,piece)
    for name,exchange in zip(('s','t','u'),result['exchanges']):
        show(name+' vector exchange coefficient',exchange)
    show('Complete odd-resonance coefficient per T_R',result['amplitude_coefficient'])
    show('Four original-fit null projections',result['null_projection'])


if __name__=='__main__':
    main()
