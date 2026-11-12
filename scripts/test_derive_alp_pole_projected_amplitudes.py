import random
import unittest
import sympy as sp
from derive_alp_pole_projected_amplitudes import ccwz_vertex, projected_amplitude, schur_propagator
from solve_alp_meson_spectrum import solve_spectrum_exact


class PoleProjectedTests(unittest.TestCase):
    def test_symbolic_schur_inverse_and_normalization(self):
        s, m, a, b, k = sp.symbols('s m a b k', real=True)
        kinetic = sp.Matrix([[1,k],[k,1+k**2]])
        mass = sp.Matrix([[m,b],[b,a]])
        direct = (s*kinetic-mass).inv()
        self.assertEqual(sp.simplify(schur_propagator(kinetic,mass,s)-direct),sp.zeros(2))
        overlap = -(s*k-b)/(s-m)
        denominator = s*(1+k**2)-a-(s*k-b)**2/(s-m)
        normalization = (sp.Matrix([overlap,1]).T*kinetic*sp.Matrix([overlap,1]))[0]
        self.assertEqual(sp.simplify(sp.diff(denominator,s)-normalization),0)

    def test_operator_inverse_and_general_mass_curvature(self):
        from derive_alp_finite_mixing_operator_inverse import (rational_inputs, lo_blocks, generator_components,
            closed_form_overlap, closed_form_alp_inverse_propagator)
        inputs=rational_inputs(77)
        kinetic,mass=lo_blocks(inputs)
        squared_momentum=sp.Rational(7,4)
        inverse=squared_momentum*kinetic-mass
        meson_inverse=inverse[:-1,:-1].inv()
        overlap=-meson_inverse*inverse[:-1,-1:]
        expected=sp.Matrix(generator_components(closed_form_overlap(inputs,squared_momentum)))
        self.assertEqual(sp.simplify(overlap-expected),sp.zeros(9,1))
        denominator=inverse[-1,-1]-(inverse[-1:,:-1]*meson_inverse*inverse[:-1,-1:])[0]
        self.assertEqual(sp.simplify(denominator-closed_form_alp_inverse_propagator(inputs,squared_momentum)),0)
        correction=sp.Rational(3,11)
        mass[9,9] += correction
        updated=schur_propagator(kinetic,mass,squared_momentum)
        self.assertEqual(sp.simplify(updated-(squared_momentum*kinetic-mass).inv()),sp.zeros(10))
        self.assertEqual(sp.simplify(1/updated[-1,-1]-(denominator-correction)),0)

    def test_full_inverse_pole_sum_and_residues(self):
        s=sp.Symbol('s')
        transform=sp.Matrix([[1,sp.Rational(1,2)],[sp.Rational(1,3),1]])
        kinetic=transform.T*transform
        mass=transform.T*sp.diag(1,16)*transform
        spectrum=solve_spectrum_exact(kinetic,mass)
        direct=schur_propagator(kinetic,mass,s)
        self.assertEqual(sp.simplify(direct-spectrum.propagator_from_poles(s)),sp.zeros(2))
        for pole,residue in zip(spectrum.pole_squared_masses,spectrum.residues):
            self.assertEqual((direct*(s-pole)).applyfunc(lambda x: sp.limit(x,s,pole)),residue)

    def test_projected_ccwz_vertex_and_field_change(self):
        z=sp.zeros(3)
        generator=sp.diag(1,0,-1)
        sources=([sp.eye(3),sp.diag(1,-1,0),z,z],[z]*4,[z]*4,[0]*4,sp.Integer(1),sp.Integer(2),0,0)
        transform=sp.Matrix([[1,sp.Rational(1,2)],[sp.Rational(1,3),1]])
        kinetic=transform.T*transform
        mass=transform.T*sp.diag(1,4)*transform
        spectrum=solve_spectrum_exact(kinetic,mass)
        gram=sp.Matrix([[4,-2,-2],[-2,1,1],[-2,1,1]])
        def original(directions,momenta,gram):
            return ccwz_vertex([generator],sources,directions,momenta,gram)
        amplitude=projected_amplitude(spectrum,[1,0,0],gram,original)
        self.assertNotEqual(amplitude,0)
        change=sp.Matrix([[2,1],[1,1]])
        transformed=solve_spectrum_exact(change.T*kinetic*change,change.T*mass*change)
        transported=[change.inv()*spectrum.pole_vector(i) for i in (1,0,0)]
        def primed(directions,momenta,gram):
            return original([change*direction for direction in directions],momenta,gram)
        self.assertEqual(amplitude,projected_amplitude(transformed,[1,0,0],gram,primed,transported))

    def test_four_leg_ccwz_contact_and_exchange(self):
        z=sp.zeros(3)
        generator=sp.diag(1,0,-1)
        sources=([sp.eye(3),sp.diag(1,-1,0),z,z,z],[z]*5,[z]*5,[0]*5,sp.S.One,sp.S(2),0,0)
        transform=sp.Matrix([[1,sp.Rational(1,2)],[sp.Rational(1,3),1]])
        kinetic=transform.T*transform
        mass=transform.T*sp.diag(1,16)*transform
        spectrum=solve_spectrum_exact(kinetic,mass)
        gram=sp.Matrix([[16,-6,-6,-4],[-6,1,sp.Rational(7,2),sp.Rational(3,2)],[-6,sp.Rational(7,2),1,sp.Rational(3,2)],[-4,sp.Rational(3,2),sp.Rational(3,2),1]])
        def original(directions,momenta,gram):
            return ccwz_vertex([generator],sources,directions,momenta,gram)
        amplitude=projected_amplitude(spectrum,[1,0,0,0],gram,original)
        directions=[spectrum.pole_vector(i) for i in (1,0,0,0)]
        momenta=[sp.eye(4)[:,i] for i in range(4)]
        self.assertNotEqual(amplitude,original(directions,momenta,gram))
        change=sp.Matrix([[2,1],[1,1]])
        transformed=solve_spectrum_exact(change.T*kinetic*change,change.T*mass*change)
        transported=[change.inv()*direction for direction in directions]
        def primed(directions,momenta,gram):
            return original([change*direction for direction in directions],momenta,gram)
        self.assertEqual(amplitude,projected_amplitude(transformed,[1,0,0,0],gram,primed,transported))

    def test_complete_ccwz_quadratic_and_four_leg_calculation(self):
        from derive_alp_lo_joint_vertices import quadratic_blocks
        zero=sp.zeros(3)
        generator=sp.diag(1,1,-2)/sp.sqrt(3)
        first_mass_derivative=(sp.diag(1,-sp.Rational(1,2),-sp.Rational(1,2))
                               +sp.I*sp.diag(1,0,sp.Rational(1,2)))
        sources=([sp.eye(3),first_mass_derivative,sp.I*sp.diag(1,0,0),zero,zero],
                 [zero]*5,[zero]*5,[0]*5,sp.S.One,sp.S(2),0,sp.S(16))
        kinetic,mass=quadratic_blocks([generator],*sources)
        self.assertEqual(kinetic,sp.eye(2))
        self.assertEqual(mass,sp.diag(1,16))
        spectrum=solve_spectrum_exact(kinetic,mass)
        gram=sp.Matrix([[16,-6,-6,-4],[-6,1,sp.Rational(7,2),sp.Rational(3,2)],[-6,sp.Rational(7,2),1,sp.Rational(3,2)],[-4,sp.Rational(3,2),sp.Rational(3,2),1]])
        def ccwz(directions,momenta,gram):
            return ccwz_vertex([generator],sources,directions,momenta,gram)
        indices=[1,0,0,0]
        directions=[spectrum.pole_vector(i) for i in indices]
        momenta=[sp.eye(4)[:,i] for i in range(4)]
        contact=ccwz(directions,momenta,gram)
        amplitude=projected_amplitude(spectrum,indices,gram,ccwz)
        self.assertEqual(contact,sp.sqrt(3)/12)
        self.assertEqual(sp.simplify(amplitude),419*sp.sqrt(3)/4928)

    def test_ccwz_vertex_adapter(self):
        z=sp.zeros(3); unit=sp.eye(3)
        gram=sp.Matrix([[1,2,3],[2,4,5],[3,5,6]])
        directions=[sp.Matrix([1,0]),sp.Matrix([1,0]),sp.Matrix([0,1])]
        sources=([unit,z,z,z],[z]*4,[z]*4,[0]*4,1,2,0,0)
        result=ccwz_vertex([unit],sources,directions,[sp.eye(3)[:,i] for i in range(3)],gram)
        self.assertEqual(result,0)
        generator=sp.diag(1,0,-1)
        sources=([unit,sp.diag(1,-1,0),z,z],[z]*4,[z]*4,[0]*4,sp.Integer(1),sp.Integer(2),0,0)
        result=ccwz_vertex([generator],sources,directions,[sp.eye(3)[:,i] for i in range(3)],gram)
        self.assertEqual(result,-sp.Rational(1,4))

    def test_on_shell_contact_exchange_under_random_field_change(self):
        seed=random.Random(419)
        transform=sp.Matrix([[1,sp.Rational(1,2)],[sp.Rational(1,3),1]])
        kinetic=transform.T*transform
        mass=transform.T*sp.diag(1,16)*transform
        spectrum=solve_spectrum_exact(kinetic,mass)
        self.assertTrue(all(spectrum.pole_vector(i)[1] != 0 for i in range(2)))
        change=sp.Matrix([[seed.randint(1,4),seed.randint(1,4)],[seed.randint(1,4),seed.randint(1,4)]])
        while change.det()==0:
            change=sp.Matrix([[seed.randint(1,4),seed.randint(1,4)],[seed.randint(1,4),seed.randint(1,4)]])
        transformed=solve_spectrum_exact(change.T*kinetic*change,change.T*mass*change)
        # Eigenvector signs are arbitrary: transport the original phase choices.
        for i in range(2):
            transported=change.inv()*spectrum.pole_vector(i)
            if transformed.pole_vector(i) == -transported:
                transformed.pole_vectors[:,i] = transported
            self.assertEqual(transformed.pole_vector(i),transported)
        gram=sp.Matrix([[16,-6,-6,-4],[-6,1,sp.Rational(7,2),sp.Rational(3,2)],[-6,sp.Rational(7,2),1,sp.Rational(3,2)],[-4,sp.Rational(3,2),sp.Rational(3,2),1]])
        self.assertEqual(gram*sp.ones(4,1),sp.zeros(4,1))
        from itertools import product
        tensors={n:{ind:sp.Integer(seed.randint(-3,3)) for ind in product(range(2),repeat=n)} for n in (3,4)}
        # Symmetrize to represent derivatives of a local Lagrangian.
        from itertools import permutations
        tensors={n:{ind:sum(tensor[p] for p in permutations(ind))/sp.factorial(n) for ind in tensor} for n,tensor in tensors.items()}
        def evaluate(tensor,directions):
            n=len(directions)
            return sum(tensor[ind]*sp.prod(directions[j][ind[j]] for j in range(n)) for ind in product(range(2),repeat=n))
        def original(directions,momenta,gram):
            return evaluate(tensors[len(directions)],directions)
        transformed_tensors={n:{ind:sum(tensor[old]*sp.prod(change[old[j],ind[j]] for j in range(n)) for old in product(range(2),repeat=n)) for ind in tensor} for n,tensor in tensors.items()}
        def primed(directions,momenta,gram):
            return evaluate(transformed_tensors[len(directions)],directions)
        original_result=projected_amplitude(spectrum,[1,0,0,0],gram,original)
        transformed_result=projected_amplitude(transformed,[1,0,0,0],gram,primed)
        self.assertNotEqual(original_result,0)
        self.assertEqual(sp.simplify(original_result-transformed_result),0)
        self.assertNotEqual(sp.simplify(original_result-original([spectrum.pole_vector(i) for i in [1,0,0,0]],None,None)),0)
        self.assertNotEqual(((sp.Rational(5,2)*kinetic-mass).inv())[0,1],0)

if __name__=='__main__': unittest.main()
