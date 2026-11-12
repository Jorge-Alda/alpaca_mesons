"""Symbolic comparisons for the general-source LO vertices."""
import unittest
import sympy as sp
from derive_alp_lo_joint_vertices import vertex, quadratic_blocks
from derive_alp_lo_neutral_blocks import closed_form_blocks, illustrative_lo_inputs, NEUTRAL_GENERATOR_INDICES
from symbolic_utils import u3_generators
from solve_alp_meson_spectrum import solve_spectrum_numerical
import numpy as np

class JointVertexTests(unittest.TestCase):
    def setUp(self):
        self.z=sp.zeros(3)
        self.fp=sp.Rational(2,3); self.fa=sp.Rational(5,4)
        self.mass=sp.diag(2,3,5)
        self.cv=sp.Matrix([[1,0,0],[0,2,1],[0,1,4]])
        self.ca=sp.Matrix([[2,0,0],[0,1,1],[0,1,3]])
        self.cs=sp.Matrix([[0,0,0],[0,0,2],[0,2,1]])
        self.cp=sp.Matrix([[1,0,0],[0,2,1],[0,1,3]])
        self.chi1=sp.I*(self.cs*self.mass-self.mass*self.cs+self.cp*self.mass+self.mass*self.cp)
        self.gs=u3_generators()

    def calculate(self, gens, alps, p, chi2=None, right1=None, theta2=0):
        n=len(gens)
        pad=[self.z]*(n+1)
        chi=[self.mass,self.chi1,chi2 if chi2 is not None else self.z]+[self.z]*n
        right=[self.cv+self.ca,right1 if right1 is not None else self.z]+[self.z]*n
        left=[self.cv-self.ca]+[self.z]*n
        return vertex(gens,alps,p,chi,right,left,[0,sp.Rational(1,3),theta2]+[0]*n,
                      self.fp,self.fa,sp.Rational(4,5),sp.Rational(7,9))

    def test_source_free_kinetic_normalization(self):
        q=sp.symbols('q0:2',real=True)
        empty=[self.z]*3
        meson=vertex([self.gs[1],self.gs[1]],[0,0],q,empty,empty,empty,[0]*3,
                     self.fp,self.fa,0,0)
        alp=vertex([self.z,self.z],[1,1],q,empty,empty,empty,[0]*3,
                   self.fp,self.fa,0,0)
        self.assertEqual(meson,-q[0]*q[1])
        self.assertEqual(alp,-q[0]*q[1])

    def test_general_blocks_feed_spectrum_solver(self):
        empty=[self.z]*3
        chi=[self.mass,self.z,self.z]
        parameters=(empty,empty,[0,0,0],self.fp,self.fa,0,7)
        kinetic,mass=quadratic_blocks([self.gs[3]],chi,*parameters)
        chi[2]=sp.diag(sp.Rational(1,10),0,0)
        kinetic_corrected,mass_corrected=quadratic_blocks([self.gs[3]],chi,*parameters)
        difference=-self.fp**2*sp.trace(chi[2])/(2*self.fa**2)
        self.assertEqual(kinetic_corrected,kinetic)
        self.assertEqual(mass_corrected-mass,sp.diag(0,difference))
        original=solve_spectrum_numerical(np.array(kinetic,dtype=float),np.array(mass,dtype=float))
        corrected=solve_spectrum_numerical(np.array(kinetic_corrected,dtype=float),np.array(mass_corrected,dtype=float))
        self.assertAlmostEqual(corrected.pole_squared_masses[-1]-original.pole_squared_masses[-1],float(difference))

    def test_quadratic_matches_archived_blocks_and_general_mass_source(self):
        args=illustrative_lo_inputs('neutral')
        self.fp=args['pion_decay_constant']; self.fa=args['alp_decay_constant']
        self.mass=args['quark_spurion']; self.ca=args['axial_source']; self.cp=args['pseudoscalar_source']; self.cv=self.z; self.cs=self.z
        self.chi1=sp.I*(self.cp*self.mass+self.mass*self.cp)
        k,m=closed_form_blocks(**args)
        gens=[self.gs[i] for i in NEUTRAL_GENERATOR_INDICES]+[self.z]
        alp=[0]*5+[1]
        chi=[self.mass,self.chi1,self.z]; currents=[self.ca,self.z,self.z]
        for i in range(6):
            for j in range(i,6):
                value=vertex([gens[i],gens[j]],[alp[i],alp[j]],[1,1],chi,currents,[-x for x in currents],
                    [0,args['gluon_coupling'],0],self.fp,self.fa,args['anomaly_mass_squared'],args['bare_alp_mass_squared'])
                self.assertEqual(sp.simplify(value + k[i,j] + m[i,j]),0, (i,j,value,k[i,j],m[i,j]))
        chi[2]=sp.diag(2,4,6)
        value=vertex([self.z,self.z],[1,1],[1,1],chi,currents,[-x for x in currents],
                    [0,args['gluon_coupling'],0],self.fp,self.fa,args['anomaly_mass_squared'],args['bare_alp_mass_squared'])
        self.assertEqual(sp.simplify(value + k[5,5] + m[5,5] - self.fp**2*sp.trace(chi[2])/(2*self.fa**2)),0)

    def test_linear_two_meson_legacy_and_permutation(self):
        q0,q1,q2=sp.symbols('q0 q1 q2',real=True)
        first,second=self.gs[3],self.gs[6]
        got=self.calculate([self.z,first,second],[1,0,0],[q0,q1,q2])
        # Legacy phi2 expression, differentiated in the two labeled fields.
        t1,t2=sp.symbols('t1 t2',real=True)
        phi=t1*first+t2*second; dphi=-sp.I*(q1*t1*first+q2*t2*second)
        expected=-sp.I*(-sp.I*q0)/(2*self.fa)*sp.trace(dphi*(self.cv*phi-phi*self.cv))
        expected+=-sp.I/(4*self.fa)*sp.trace(phi**2*(self.cs*self.mass-self.mass*self.cs))
        expected=sp.diff(expected,t1,t2).subs({t1:0,t2:0})
        self.assertEqual(sp.simplify(got-expected),0)
        self.assertEqual(sp.simplify(got-self.calculate([second,self.z,first],[0,1,0],[q2,q0,q1])),0)

    def test_linear_three_meson_legacy(self):
        q=sp.symbols('q0:4', real=True); t=sp.symbols('t0:3',real=True)
        gens=[self.gs[3],self.gs[6],self.gs[6]]
        phi=sum((t[i]*gens[i] for i in range(3)),self.z)
        dphi=sum((-sp.I*q[i+1]*t[i]*gens[i] for i in range(3)),self.z)
        comm=lambda a,b:a*b-b*a
        expected=(-sp.I*q[0])/(6*self.fa*self.fp)*sp.trace(self.ca*comm(comm(dphi,phi),phi))
        expected-=sp.trace(phi**3*(self.cp*self.mass+self.mass*self.cp))/(12*self.fa*self.fp)
        expected=sp.diff(expected,*t).subs({x:0 for x in t})
        got=self.calculate([self.z]+gens,[1,0,0,0],q)
        self.assertEqual(sp.simplify(got-expected),0)

    def test_nonlinear_source_and_charged_legs(self):
        q=sp.symbols('q0:4',real=True)
        plus=(self.gs[1]+sp.I*self.gs[2])/sp.sqrt(2)
        minus=(self.gs[1]-sp.I*self.gs[2])/sp.sqrt(2)
        nonlinear=sp.diag(1,2,4)
        base=self.calculate([self.z,self.z,plus,minus],[1,1,0,0],q)
        changed=self.calculate([self.z,self.z,plus,minus],[1,1,0,0],q,chi2=nonlinear,right1=nonlinear,theta2=2)
        self.assertNotEqual(sp.simplify(self.calculate([self.z,self.z,plus,minus],[1,1,0,0],q,chi2=nonlinear)-base),0)
        neutral=[self.z,self.z,self.gs[3]]
        self.assertNotEqual(sp.simplify(self.calculate(neutral,[1,1,0],q[:3],right1=nonlinear)-self.calculate(neutral,[1,1,0],q[:3])),0)
        self.assertNotEqual(sp.simplify(self.calculate([self.z]*3,[1]*3,q[:3],theta2=2)-self.calculate([self.z]*3,[1]*3,q[:3])),0)
        perm=self.calculate([minus,self.z,plus,self.z],[0,1,0,1],[q[3],q[0],q[2],q[1]],chi2=nonlinear,right1=nonlinear,theta2=2)
        self.assertEqual(sp.simplify(changed-perm),0)


class GeneralGeneratorIdentities(unittest.TestCase):
    def test_general_symbolic_quadratic_blocks(self):
        """Compare each quadratic block with the archived closed forms."""
        z=sp.zeros(3)
        def hermitian(prefix):
            diagonal=sp.symbols(f'{prefix}0:3',real=True)
            real=sp.symbols(f'{prefix}r0:3',real=True)
            imag=sp.symbols(f'{prefix}i0:3',real=True)
            return sp.Matrix([[diagonal[0],real[0]+sp.I*imag[0],real[1]+sp.I*imag[1]],
                              [real[0]-sp.I*imag[0],diagonal[1],real[2]+sp.I*imag[2]],
                              [real[1]-sp.I*imag[1],real[2]-sp.I*imag[2],diagonal[2]]])
        first,second=hermitian('g'),hermitian('h')
        axial,pseudoscalar=hermitian('A'),hermitian('P')
        mass=sp.diag(*sp.symbols('m0:3',real=True))
        f_pi,f_a,m0,ma,cg=sp.symbols('f_pi f_a m0 ma cg',positive=True)
        q=sp.symbols('q0:2',real=True)
        chi1=sp.I*(pseudoscalar*mass+mass*pseudoscalar)
        chi=[mass,chi1,z]
        right=[axial,z,z]; left=[-axial,z,z]
        anti=lambda a,b:a*b+b*a
        def calculate(generators,alps):
            return vertex(generators,alps,q,chi,right,left,[0,cg,0],f_pi,f_a,m0,ma)
        mm= -q[0]*q[1]*sp.trace(first*second)/2-sp.trace(anti(first,second)*mass)/4
        mm-=m0*sp.trace(first)*sp.trace(second)/6
        mixed=q[0]*q[1]*f_pi/f_a*sp.trace(axial*first)
        mixed+=f_pi/(2*f_a)*sp.trace(first*anti(pseudoscalar,mass))
        mixed-=f_pi*m0*cg/(3*f_a)*sp.trace(first)
        aa=-q[0]*q[1]*(1+2*f_pi**2/f_a**2*sp.trace(axial*axial))
        aa-=ma+2*f_pi**2*m0*cg**2/(3*f_a**2)
        for actual,expected in ((calculate([first,second],[0,0]),mm),
                                (calculate([first,z],[0,1]),mixed),
                                (calculate([z,z],[1,1]),aa)):
            self.assertEqual(sp.expand(actual-expected),0)

    def test_complex_strangeness_changing_vector_with_charged_mesons(self):
        """Keep both independent components of the Hermitian d-s source."""
        z=sp.zeros(3)
        real,imag=sp.symbols('c_ds_real c_ds_imag',real=True)
        vector=sp.Matrix([[0,0,0],[0,0,real+sp.I*imag],[0,real-sp.I*imag,0]])
        pion=sp.Matrix([[0,1,0],[0,0,0],[0,0,0]])
        kaon=sp.Matrix([[0,0,0],[0,0,0],[1,0,0]])
        f_pi,f_a=sp.symbols('f_pi f_a',positive=True)
        q=sp.symbols('q0:3',real=True)
        got=vertex([z,pion,kaon],[1,0,0],q,[z]*4,
                   [vector,z,z,z],[vector,z,z,z],[0]*4,f_pi,f_a,0,0)
        t1,t2=sp.symbols('t1 t2',real=True)
        phi=t1*pion+t2*kaon
        dphi=-sp.I*(q[1]*t1*pion+q[2]*t2*kaon)
        legacy=-sp.I*(-sp.I*q[0])/(2*f_a)*sp.trace(dphi*(vector*phi-phi*vector))
        expected=sp.diff(legacy,t1,t2).subs({t1:0,t2:0})
        self.assertEqual(sp.expand(got-expected),0)
        self.assertNotEqual(sp.diff(got,real),0)
        self.assertNotEqual(sp.diff(got,imag),0)

    def test_three_meson_linear_sources_with_generic_generators(self):
        """Identity for generic 3x3 generators and Hermitian sources."""
        z=sp.zeros(3)
        generators=[sp.Matrix(3,3,lambda i,j:sp.Symbol(f'g{k}_{i}{j}',real=True)) for k in range(3)]
        mass=sp.diag(*sp.symbols('m0:3',real=True))
        def hermitian(prefix):
            real=sp.Matrix(3,3,lambda i,j:sp.Symbol(f'{prefix}r{i}{j}',real=True))
            imag=sp.Matrix(3,3,lambda i,j:sp.Symbol(f'{prefix}i{i}{j}',real=True))
            return real+real.T+sp.I*(imag-imag.T)
        axial=hermitian('A')
        pseudoscalar=hermitian('P')
        f_pi,f_a=sp.symbols('f_pi f_a',positive=True)
        q=sp.symbols('q0:4',real=True)
        chi1=sp.I*(pseudoscalar*mass+mass*pseudoscalar)
        got=vertex([z]+generators,[1,0,0,0],q,[mass,chi1,z,z,z],
                   [axial,z,z,z,z],[-axial,z,z,z,z],[0]*5,f_pi,f_a,0,0)
        t=sp.symbols('t0:3',real=True)
        phi=sum((t[i]*generators[i] for i in range(3)),z)
        dphi=sum((-sp.I*q[i+1]*t[i]*generators[i] for i in range(3)),z)
        comm=lambda a,b:a*b-b*a
        legacy=(-sp.I*q[0])/(6*f_a*f_pi)*sp.trace(axial*comm(comm(dphi,phi),phi))
        legacy-=sp.trace(phi**3*(pseudoscalar*mass+mass*pseudoscalar))/(12*f_a*f_pi)
        expected=sp.diff(legacy,*t).subs({x:0 for x in t})
        self.assertEqual(sp.expand(got-expected),0)

    def test_general_second_mass_derivative(self):
        """The source correction holds for a generic complex 3x3 derivative."""
        z=sp.zeros(3)
        real=sp.Matrix(3,3,lambda i,j:sp.Symbol(f'r{i}{j}',real=True))
        imag=sp.Matrix(3,3,lambda i,j:sp.Symbol(f'i{i}{j}',real=True))
        chi2=real+sp.I*imag
        f_pi,f_a=sp.symbols('f_pi f_a',positive=True)
        output=vertex([z,z],[1,1],[1,1],[z,z,chi2],[z]*3,[z]*3,[0]*3,f_pi,f_a,0,0)
        self.assertEqual(sp.expand(output+1-f_pi**2*sp.trace(real)/(2*f_a**2)),0)

    def test_two_meson_linear_sources_with_generic_generators(self):
        z=sp.zeros(3)
        generators=[sp.Matrix(3,3,lambda i,j:sp.Symbol(f'h{k}_{i}{j}',real=True)) for k in range(2)]
        mass=sp.diag(*sp.symbols('m0:3',real=True))
        v=sp.Matrix(3,3,lambda i,j:sp.Symbol(f'v{i}{j}',real=True))
        scalar=sp.Matrix(3,3,lambda i,j:sp.Symbol(f's{i}{j}',real=True))
        v=v+v.T; scalar=scalar+scalar.T
        f_pi,f_a=sp.symbols('f_pi f_a',positive=True)
        q=sp.symbols('q0:3',real=True)
        chi1=sp.I*(scalar*mass-mass*scalar)
        got=vertex([z]+generators,[1,0,0],q,[mass,chi1,z,z],
                   [v,z,z,z],[v,z,z,z],[0]*4,f_pi,f_a,0,0)
        t=sp.symbols('t0:2',real=True)
        phi=sum((t[i]*generators[i] for i in range(2)),z)
        dphi=sum((-sp.I*q[i+1]*t[i]*generators[i] for i in range(2)),z)
        legacy=-sp.I*(-sp.I*q[0])/(2*f_a)*sp.trace(dphi*(v*phi-phi*v))
        legacy-=sp.I/(4*f_a)*sp.trace(phi**2*(scalar*mass-mass*scalar))
        expected=sp.diff(legacy,*t).subs({x:0 for x in t})
        self.assertEqual(sp.expand(got-expected),0)

if __name__=='__main__': unittest.main()
