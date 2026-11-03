"""Derive the general finite-mixing LO a -> 3pi amplitude in CCWZ notation.

Run: .venv/bin/python scripts/derive_alp_finite_three_pion.py
Tests: .venv/bin/python scripts/derive_alp_finite_three_pion.py --test
Inputs: memory/alp/shared/finite-alp-meson-mixing.md;
        memory/nlo/lsz-invariant-split.md;
        memory/extracts/Alda:2026nrd/appendix-a-expansion.md.
Natural units, mostly-minus metric, S=1+iT, tr(lambda_i lambda_j)=2 delta_ij.
All flavor matrices remain general and noncommuting. The calculation uses
vacuum derivatives of chi_+ and chi_- after removing the derivative sources.
The exact finite chiral change of variables retains every induced derivative
through fourth order. No numerical source matrix or spectrum is selected.
"""
import itertools
import sys
import unittest
import sympy as sp

f_pi,f_a,B_0=sp.symbols("f_pi f_a B_0",positive=True)


def commutator(first,second):
    return first*second-second*first


def anticommutator(first,second):
    return first*second+second*first


def show(label,expression):
    print(f"{label}:\n{sp.latex(expression)}\n")


def ccwz_vacuum_derivatives(axial,vector,original_plus_derivatives,original_minus_derivatives):
    """Derive vacuum CCWZ derivatives after removing both derivative sources.

    The arguments are derivatives with respect to a, at a=Phi=0. The
    exponential products implement the full finite chiral transformation;
    they retain matrix ordering and arbitrary original higher derivatives.
    """
    transformed_plus=[]
    transformed_minus=[]
    size=original_plus_derivatives[0].rows
    for order in range(len(original_plus_derivatives)):
        first_product=sp.zeros(size)
        second_product=sp.zeros(size)
        for left_order in range(order+1):
            for source_order in range(order-left_order+1):
                right_order=order-left_order-source_order
                multiplicity=sp.factorial(order)/(sp.factorial(left_order)*sp.factorial(source_order)*sp.factorial(right_order))
                first_product += multiplicity*((-sp.I*(vector+axial)/f_a)**left_order
                    *(original_plus_derivatives[source_order]+original_minus_derivatives[source_order])
                    *(sp.I*(vector-axial)/f_a)**right_order)
                second_product += multiplicity*((-sp.I*(vector-axial)/f_a)**left_order
                    *(original_plus_derivatives[source_order]-original_minus_derivatives[source_order])
                    *(sp.I*(vector+axial)/f_a)**right_order)
        transformed_plus.append(sp.expand((first_product+second_product)/2))
        transformed_minus.append(sp.expand((first_product-second_product)/2))
    return transformed_plus,transformed_minus


def linear_source_ccwz_derivatives(axial,pseudoscalar,vector,scalar,quark_mass_matrix,order=4):
    """Evaluate exact derivatives for the linear source coordinates in memory."""
    vacuum_mass_source=2*B_0*quark_mass_matrix
    plus=[2*vacuum_mass_source,2*sp.I*commutator(scalar,vacuum_mass_source)/f_a]
    minus=[sp.zeros(quark_mass_matrix.rows),2*sp.I*anticommutator(pseudoscalar,vacuum_mass_source)/f_a]
    plus += [sp.zeros(quark_mass_matrix.rows)]*(order-1)
    minus += [sp.zeros(quark_mass_matrix.rows)]*(order-1)
    return ccwz_vacuum_derivatives(axial,vector,plus,minus)


def cubic_vertex(first,second,third,chi_plus_first_derivative,
                 chi_minus_second_derivative,chi_plus_third_derivative):
    """Pole=(lambda_i w_i, w_hat_a); derivatives are at a=Phi=0."""
    first_matrix,first_alp=first
    second_matrix,second_alp=second
    third_matrix,third_alp=third
    return (-sp.Trace(chi_plus_first_derivative*(
        first_alp*anticommutator(second_matrix,third_matrix)
        +second_alp*anticommutator(first_matrix,third_matrix)
        +third_alp*anticommutator(first_matrix,second_matrix)))/8
        +f_pi*sp.Trace(chi_minus_second_derivative*(
            first_alp*second_alp*third_matrix+first_alp*third_alp*second_matrix
            +second_alp*third_alp*first_matrix))/(4*sp.I)
        +f_pi**2*first_alp*second_alp*third_alp*sp.Trace(chi_plus_third_derivative)/4)


def cubic_pair_source(first,second,chi_plus_first_derivative,
                      chi_minus_second_derivative,chi_plus_third_derivative):
    """Return the meson and ALP coefficients of a third internal leg."""
    first_matrix,first_alp=first
    second_matrix,second_alp=second
    meson=-(first_alp*anticommutator(chi_plus_first_derivative,second_matrix)
            +second_alp*anticommutator(chi_plus_first_derivative,first_matrix))/4
    meson += f_pi*first_alp*second_alp*chi_minus_second_derivative/(2*sp.I)
    alp=(-sp.Trace(chi_plus_first_derivative*anticommutator(first_matrix,second_matrix))/8
         +f_pi*(first_alp*sp.Trace(chi_minus_second_derivative*second_matrix)
                +second_alp*sp.Trace(chi_minus_second_derivative*first_matrix))/(4*sp.I)
         +f_pi**2*first_alp*second_alp*sp.Trace(chi_plus_third_derivative)/4)
    return meson,alp


def inverse_meson_mass_operator(source,quark_mass_matrix,anomaly_mass_sq,invariant_mass_sq):
    """Solve the matrix form of (M_meson^2-s) without diagonalizing flavor."""
    size=quark_mass_matrix.rows
    identity=sp.eye(size)
    mass_matrix=B_0*quark_mass_matrix-invariant_mass_sq*identity/2
    vector_identity=identity.vec()
    operator=(sp.kronecker_product(identity,mass_matrix)
              +sp.kronecker_product(mass_matrix.T,identity)
              +anomaly_mass_sq*vector_identity*vector_identity.T/size)
    solution=operator.inv()*source.vec()
    return sp.Matrix(size,size,lambda row,column:solution[column*size+row])


def inverse_mixed_mass_operator(pair_source,invariant_coupling,alp_mass_entry,
                                inverse_meson_operator,invariant_mass_sq):
    """Exact inverse of M^2-sK in derivative-source-free coordinates, K=1."""
    meson_source,alp_source=pair_source
    inverse_source=inverse_meson_operator(meson_source)
    inverse_coupling=inverse_meson_operator(invariant_coupling)
    denominator=(alp_mass_entry-invariant_mass_sq
                 -2*f_pi**2*sp.Trace(invariant_coupling*inverse_coupling)/f_a**2)
    alp_field=(alp_source+f_pi*sp.Trace(invariant_coupling*inverse_source)/f_a)/denominator
    meson_field=inverse_source+2*f_pi*inverse_coupling*alp_field/f_a
    return meson_field,alp_field


def quartic_mass_vertex(poles,chi_plus_vacuum,chi_minus_first_derivative,
                        chi_plus_second_derivative,chi_minus_third_derivative,
                        chi_plus_fourth_derivative):
    """Differentiate the CCWZ mass term four times with ordered flavor products."""
    matrices=[pole[0] for pole in poles]
    alp_overlaps=[pole[1] for pole in poles]
    result=sum((sp.Trace(chi_plus_vacuum*sp.MatMul(*[matrices[i] for i in order]))
                for order in itertools.permutations(range(4))),sp.S.Zero)/(96*f_pi**2)
    for alp_leg in range(4):
        remaining=[i for i in range(4) if i != alp_leg]
        result -= alp_overlaps[alp_leg]*sum((sp.Trace(chi_minus_first_derivative
            *sp.MatMul(*[matrices[i] for i in order]))
            for order in itertools.permutations(remaining)),sp.S.Zero)/(24*sp.I*f_pi)
    for alp_legs in itertools.combinations(range(4),2):
        meson_legs=[i for i in range(4) if i not in alp_legs]
        result -= (sp.prod(alp_overlaps[i] for i in alp_legs)
                   *sp.Trace(chi_plus_second_derivative*anticommutator(*[matrices[i] for i in meson_legs]))/8)
    for meson_leg in range(4):
        result += (f_pi*sp.prod(alp_overlaps[i] for i in range(4) if i != meson_leg)
                   *sp.Trace(chi_minus_third_derivative*matrices[meson_leg])/(4*sp.I))
    result += f_pi**2*sp.prod(alp_overlaps)*sp.Trace(chi_plus_fourth_derivative)/4
    return result


def quartic_kinetic_vertex(poles,scalar_products):
    """Expand f_pi^2 tr(u_mu u^mu)/4 with all four momenta incoming."""
    matrices=[pole[0] for pole in poles]
    result=sp.S.Zero
    for first,second,third,fourth in itertools.permutations(range(4)):
        result -= scalar_products[tuple(sorted((second,fourth)))]*sp.Trace(
            commutator(matrices[first],matrices[second])
            *commutator(matrices[third],matrices[fourth]))/(48*f_pi**2)
    return result


def quartic_kinetic_channel(poles,mandelstam_s,mandelstam_t,mandelstam_u):
    first,second,third,fourth=[pole[0] for pole in poles]
    return (-(sp.Trace(commutator(first,second)*commutator(third,fourth))*(mandelstam_t-mandelstam_u)
              +sp.Trace(commutator(first,third)*commutator(second,fourth))*(mandelstam_s-mandelstam_u)
              +sp.Trace(commutator(first,fourth)*commutator(second,third))*(mandelstam_s-mandelstam_t))
            /(24*f_pi**2)).doit()


def main():
    print("All chi_+ and chi_- vacuum derivatives below are evaluated at a=Phi=0 after the exact finite chiral field redefinition removing derivative sources.")
    print("u_mu,chi_+,chi_- transform with h on both sides; X_TC is invariant. Their traces in the LO Lagrangian are therefore unchanged.")
    print("The finite transformation uses exp[-ia(cV+cA)/fa] on the right-flavor side and exp[-ia(cV-cA)/fa] on the left-flavor side.")
    print("Its complete exponentials, including every induced higher-a term, are retained; no new source coefficients are independently chosen.")
    print("The transformed vacuum chi_+ is half the sum of:")
    print("exp[-ia(cV+cA)/fa](chi_++chi_-) exp[ia(cV-cA)/fa],")
    print("exp[-ia(cV-cA)/fa](chi_+-chi_-) exp[ia(cV+cA)/fa].")
    print("The transformed vacuum chi_- is half their difference. Products on the right use original vacuum CCWZ building blocks at the same a.")
    print("X_TC=-tr(Phi)/f_pi-2cGA a/fa in these coordinates; cGA=cG+tr(cA).")
    print("The first derivatives are d_a chi_+=-4i B0[cVS,mq]/fa and d_a chi_-=-4i B0{cAP,mq}/fa.")
    print("Higher derivatives follow by direct differentiation of the displayed exponentials; the program evaluates them for the original linear-source coordinates.")
    axial,pseudoscalar,vector,scalar,quark_masses=[sp.MatrixSymbol(name,3,3)
        for name in ("c_q^A","c_q^P","c_q^V","c_q^S","m_q")]
    plus_derivatives,minus_derivatives=linear_source_ccwz_derivatives(axial,pseudoscalar,vector,scalar,quark_masses)
    for order in range(1,5):
        show(f"Derivative order {order} of chi_+",plus_derivatives[order])
        show(f"Derivative order {order} of chi_-",minus_derivatives[order])
    plus=[sp.MatrixSymbol(r"\left.\partial_a^{"+str(n)+r"}\chi_+\right|_0",3,3) for n in range(5)]
    minus=[sp.MatrixSymbol(r"\left.\partial_a^{"+str(n)+r"}\chi_-\right|_0",3,3) for n in range(5)]
    phi=sp.MatrixSymbol("Phi",3,3)
    alp=sp.Symbol("a",real=True)
    cubic=(-alp*sp.Trace(plus[1]*phi**2)/8+f_pi*alp**2*sp.Trace(minus[2]*phi)/(8*sp.I)
           +f_pi**2*alp**3*sp.Trace(plus[3])/24)
    quartic=(sp.Trace(plus[0]*phi**4)/(96*f_pi**2)-alp*sp.Trace(minus[1]*phi**3)/(24*sp.I*f_pi)
             -alp**2*sp.Trace(plus[2]*phi**2)/16+f_pi*alp**3*sp.Trace(minus[3]*phi)/(24*sp.I)
             +f_pi**2*alp**4*sp.Trace(plus[4])/96)
    show("Cubic part of f_pi^2 tr(chi_+)/4",cubic)
    show("Quartic part of f_pi^2 tr(chi_+)/4",quartic)
    print("V3(v1,v2,v3) = [d/dx1 d/dx2 d/dx3 L3(Phi=sum_r x_r lambda_i w_i^(r), a=sum_r x_r w_hat-a^(r))]_(x=0).")
    print("V4_chi(v1,v2,v3,v4) = [d/dx1 d/dx2 d/dx3 d/dx4 L4_chi(Phi=sum_r x_r lambda_i w_i^(r), a=sum_r x_r w_hat-a^(r))]_(x=0).")
    print("V4_kin(v1,v2,v3,v4) is the corresponding fourfold derivative of L4_kin, with derivative pairs replaced by -k_i.k_j for all incoming momenta.")
    print("V4 = V4_kin + V4_chi. M = [V4 + V3(M^2-q_s^2K)^(-1)V3 + V3(M^2-q_t^2K)^(-1)V3 + V3(M^2-q_u^2K)^(-1)V3] / product of four external LSZ denominators.")
    print("The quartic part of f_pi^2 tr(u_mu u^mu)/4 is tr([Phi,dPhi][Phi,dPhi])/(48f_pi^2).")
    print("C=-B0{mq,cAP}-(m0^2/Nf)cGA 1. K=1 in these exact coordinates.")
    print("Mij^2=B0 tr({lambda_i,lambda_j}mq)/2+m0^2 tr(lambda_i)tr(lambda_j)/(2Nf).")
    print("Mi,hat-a^2=-(f_pi/fa)tr(C lambda_i).")
    print("Mhat-a,hat-a^2=ma0^2-f_pi^2 tr(d_a^2 chi_+)/4+2f_pi^2m0^2cGA^2/(Nf fa^2).")
    print("For each nonzero bare-ALP-overlap neutral pole use Alda's convention w_hat-a^(n)=1:")
    print("w_i^(n)=(f_pi/fa)tr(lambda_i M_mn^(-1)(C)).")
    print("Its LSZ denominator is sqrt[1+(2f_pi^2/fa^2)tr((M_mn^(-1)(C))^2)].")
    print("Here M_mn[X]=B0{mq,X}+(m0^2/Nf)tr(X)1-mn^2 X is Alda's pole operator.")
    print("The pole satisfies Mhat-a,hat-a^2-mn^2-(2f_pi^2/fa^2)tr(C M_mn^(-1)(C))=0.")
    print("Explicit matrix inverse: M_s^(-1)(C) = integral_0^infinity d tau exp[-tau(B0 mq-s1/2)] C exp[-tau(B0 mq-s1/2)]")
    print(" - {(m0^2/Nf) tr(integral_0^infinity d tau exp[-tau(B0 mq-s1/2)] C exp[-tau(B0 mq-s1/2)]) / [1+(m0^2/Nf)tr((2B0 mq-s1)^(-1))]} (2B0 mq-s1)^(-1).")
    print("The integral is first evaluated where it converges and analytically continued; the matrix inverse code gives the same rational continuation directly.")
    print("Pure-meson poles have w_hat-a=0, (M_meson^2-mn^2)w=0 and tr(C lambda_i)w_i=0.")
    poles=[(sp.MatrixSymbol(r"\lambda_i w_i^{("+name+r")}",3,3),sp.Symbol(r"w_{\hat a}^{("+name+r")}"))
           for name in ("a","0")]
    poles += [(sp.MatrixSymbol(r"\lambda_{\pi^+}",3,3),0),(sp.MatrixSymbol(r"\lambda_{\pi^-}",3,3),0)]
    show("V3(v1,v2,v3) =",cubic_vertex(*poles[:3],plus[1],minus[2],plus[3]))
    pair=cubic_pair_source(*poles[:2],plus[1],minus[2],plus[3])
    show("Meson coefficient of internal cubic leg",pair[0])
    show("ALP coefficient of internal cubic leg",pair[1])
    mandelstam_s,mandelstam_t,mandelstam_u=sp.symbols("s t u")
    for channel,channel_poles in (("pi0 pi+ pi-",poles),("3pi0",[poles[0]]+[poles[1]]*3)):
        show(channel+": V4_kin =",quartic_kinetic_channel(channel_poles,mandelstam_s,mandelstam_t,mandelstam_u))
        show(channel+": V4_chi =",quartic_mass_vertex(channel_poles,plus[0],minus[1],plus[2],minus[3],plus[4]))
    inverse_at_alp=sp.MatrixSymbol(r"\mathbb M_{m_a^2}^{-1}(C)",3,3)
    inverse_at_pion=sp.MatrixSymbol(r"\mathbb M_{m_{\pi^0}^2}^{-1}(C)",3,3)
    explicit_poles=[(2*f_pi*inverse_at_alp/f_a,1),(2*f_pi*inverse_at_pion/f_a,1)]+poles[2:]
    alp_denominator=sp.sqrt(1+2*f_pi**2*sp.Trace(inverse_at_alp**2)/f_a**2)
    pion_denominator=sp.sqrt(1+2*f_pi**2*sp.Trace(inverse_at_pion**2)/f_a**2)
    reference_pion=(sp.MatrixSymbol(r"\lambda_{\pi^0}",3,3),0)
    invariant_axial_pseudoscalar=sp.MatrixSymbol(r"c_q^{AP}",3,3)
    invariant_first_minus=-4*sp.I*B_0*anticommutator(invariant_axial_pseudoscalar,quark_masses)/f_a
    for channel,channel_poles in (
        ("pi0 pi+ pi-",[explicit_poles[0],reference_pion]+poles[2:]),
        ("3pi0",[explicit_poles[0]]+[reference_pion]*3)):
        linear_contact=(quartic_kinetic_channel(channel_poles,mandelstam_s,mandelstam_t,mandelstam_u)
                        +quartic_mass_vertex(channel_poles,plus_derivatives[0],invariant_first_minus,
                            sp.ZeroMatrix(3,3),sp.ZeroMatrix(3,3),sp.ZeroMatrix(3,3)))
        show(channel+": M_linear =",linear_contact)
    for channel,channel_poles,denominator in (
        ("pi0 pi+ pi-",explicit_poles,alp_denominator*pion_denominator),
        ("3pi0",[explicit_poles[0]]+[explicit_poles[1]]*3,alp_denominator*pion_denominator**3)):
        contact=(quartic_kinetic_channel(channel_poles,mandelstam_s,mandelstam_t,mandelstam_u)
                 +quartic_mass_vertex(channel_poles,plus_derivatives[0],minus_derivatives[1],
                     plus_derivatives[2],minus_derivatives[3],plus_derivatives[4]))/denominator
        show(channel+" fully substituted finite contact",contact)
    for first_index,second_index in ((0,1),(2,3),(0,2),(1,3),(0,3),(1,2),(1,1)):
        meson_coefficient,alp_coefficient=cubic_pair_source(explicit_poles[first_index],explicit_poles[second_index],
            plus_derivatives[1],minus_derivatives[2],plus_derivatives[3])
        show(f"Fully substituted internal meson cubic coefficient for external legs {first_index},{second_index}",meson_coefficient)
        show(f"Fully substituted internal ALP cubic coefficient for external legs {first_index},{second_index}",alp_coefficient)
    print("Charged poles are (lambda_pi+,0) and (lambda_pi-,0); an electrically neutral ALP cannot mix bilinearly with them.")
    print("The full amplitude is the four-leg contact plus the three pairings of cubic vertices contracted with (M^2-q^2K)^(-1).")
    print("For 3pi0 the pairings use (a,0)|(0,0) at each of s,t,u; for pi0 pi+ pi- they are (a,0)|(+,-), (a,+)|(0,-), (a,-)|(0,+).")
    print("The exchange sign is positive because the inverse uses M^2-q^2K; all vertices and pole components are determined by the invariant CCWZ derivatives and C.")
    print("Quantified comparison, at common Mandelstam variables: the linear incoming pole is (2f_pi M_ma^(-1)(C)/fa,1), and the reference outgoing pion is (lambda_pi0,0).")
    print("The incoming finite pole is this same column divided by sqrt[1+2f_pi^2 tr((M_ma^(-1)(C))^2)/fa^2].")
    print("The outgoing finite pion differs from the reference by ( (2f_pi/fa)M_mpi^(-1)(C)/sqrt[1+2f_pi^2 tr((M_mpi^(-1)(C))^2)/fa^2]-lambda_pi0, 1/sqrt[1+2f_pi^2 tr((M_mpi^(-1)(C))^2)/fa^2] ).")
    print("For pi0 pi+ pi-, M_finite-M_linear equals (1/incoming-denominator-1)M_linear plus V4(linear incoming, the displayed outgoing correction,lambda_pi+,lambda_pi-)/incoming-denominator plus all normalized exchanges.")
    print("For 3pi0, expand the full V4 in the three outgoing corrections at their labeled momenta: three terms with one correction, three terms with two corrections, one term with three corrections, all divided by the incoming denominator, plus its normalization correction and all exchanges.")
    incoming_column=r"\left(\frac{2f_\pi}{f_a}\mathbb M_{m_a^2}^{-1}(C),1\right)"
    incoming_norm=r"\sqrt{1+\frac{2f_\pi^2}{f_a^2}\operatorname{tr}[(\mathbb M_{m_a^2}^{-1}(C))^2]}"
    pion_norm=r"\sqrt{1+\frac{2f_\pi^2}{f_a^2}\operatorname{tr}[(\mathbb M_{m_{\pi^0}^2}^{-1}(C))^2]}"
    pion_correction=(r"\left(\frac{2f_\pi\mathbb M_{m_{\pi^0}^2}^{-1}(C)}{f_a"+pion_norm
                     +r"}-\lambda_{\pi^0},\frac{1}{"+pion_norm+r"}\right)")
    neutral_reference=r"(\lambda_{\pi^0},0)"
    print(r"\mathcal M_{0+-}^{\rm finite}-\mathcal M_{0+-}^{\rm linear}="
          +r"\left(\frac1{"+incoming_norm+r"}-1\right)\mathcal M_{0+-}^{\rm linear}"
          +r"+\frac{\mathcal V^{(4)}["+incoming_column+","+pion_correction
          +r",(\lambda_{\pi^+},0),(\lambda_{\pi^-},0)]}{"+incoming_norm+r"}"
          +r"+\mathcal M_{0+-}^{\rm exchange}.")
    neutral_corrections=[]
    for corrected_legs in ((0,),(1,),(2,),(0,1),(0,2),(1,2),(0,1,2)):
        outgoing=[pion_correction if leg in corrected_legs else neutral_reference for leg in range(3)]
        neutral_corrections.append(r"\mathcal V^{(4)}["+incoming_column+","+",".join(outgoing)+"]")
    print(r"\mathcal M_{000}^{\rm finite}-\mathcal M_{000}^{\rm linear}="
          +r"\left(\frac1{"+incoming_norm+r"}-1\right)\mathcal M_{000}^{\rm linear}"
          +r"+\frac{"+"+".join(neutral_corrections)+r"}{"+incoming_norm+r"}+\mathcal M_{000}^{\rm exchange}.")
    print("Away from degeneracy, incoming normalization-1 and the outgoing meson correction are O((f_pi/fa)^2); the outgoing bare-ALP projection is O(f_pi/fa). The linear amplitude is O(f_pi/fa), and each finite-minus-linear contribution starts at O((f_pi/fa)^3), giving relative O((f_pi/fa)^2). This counting fails near a mixed pole degeneracy.")
    print("The small-mixing comparison uses the external-pion phase chosen continuously so its normalized meson projection tends to lambda_pi0.")


class CcwzThreePionTests(unittest.TestCase):
    def test_first_derivatives(self):
        axial=sp.Matrix([[1,2],[2,3]])
        pseudoscalar=sp.Matrix([[2,1],[1,-1]])
        vector=sp.Matrix([[0,sp.I],[-sp.I,4]])
        scalar=sp.Matrix([[1,-2],[-2,0]])
        masses=sp.Matrix([[3,1],[1,7]])
        plus,minus=linear_source_ccwz_derivatives(axial,pseudoscalar,vector,scalar,masses)
        self.assertEqual(sp.simplify(plus[1]+4*sp.I*B_0*commutator(vector-scalar,masses)/f_a),sp.zeros(2))
        self.assertEqual(sp.simplify(minus[1]+4*sp.I*B_0*anticommutator(axial-pseudoscalar,masses)/f_a),sp.zeros(2))

    def test_higher_original_derivative_is_retained(self):
        vacuum=sp.Matrix([[2,1],[1,5]])
        second=sp.Matrix([[3,1],[1,4]])
        plus,minus=ccwz_vacuum_derivatives(sp.zeros(2),sp.zeros(2),[vacuum,sp.zeros(2),second],[sp.zeros(2)]*3)
        self.assertEqual(plus[2],second)
        self.assertEqual(minus[2],sp.zeros(2))

    def test_finite_derivatives_match_exponential_series(self):
        axial=sp.Matrix([[1,2],[2,-1]])
        vector=sp.Matrix([[0,1],[1,3]])
        original_plus=[sp.Matrix([[2,1],[1,5]]),sp.Matrix([[3,sp.I],[-sp.I,1]])]+[sp.zeros(2)]*3
        original_minus=[sp.zeros(2),sp.I*sp.Matrix([[1,2],[2,3]])]+[sp.zeros(2)]*3
        alp=sp.Symbol("a")
        first_left=sum(((-sp.I*alp*(vector+axial)/f_a)**n/sp.factorial(n) for n in range(5)),sp.zeros(2))
        first_right=sum(((sp.I*alp*(vector-axial)/f_a)**n/sp.factorial(n) for n in range(5)),sp.zeros(2))
        second_left=sum(((-sp.I*alp*(vector-axial)/f_a)**n/sp.factorial(n) for n in range(5)),sp.zeros(2))
        second_right=sum(((sp.I*alp*(vector+axial)/f_a)**n/sp.factorial(n) for n in range(5)),sp.zeros(2))
        first=sp.expand(first_left*(original_plus[0]+alp*(original_plus[1]+original_minus[1]))*first_right)
        second=sp.expand(second_left*(original_plus[0]+alp*(original_plus[1]-original_minus[1]))*second_right)
        plus,minus=ccwz_vacuum_derivatives(axial,vector,original_plus,original_minus)
        for order in range(5):
            for expanded,actual in (((first+second)/2,plus[order]),((first-second)/2,minus[order])):
                coefficient=expanded.applyfunc(lambda entry:sp.expand(entry).coeff(alp,order)*sp.factorial(order))
                self.assertEqual(sp.simplify(coefficient-actual),sp.zeros(2))

    def test_cubic_internal_coefficient(self):
        matrices=[sp.Matrix([[n,1],[1,-n]]) for n in (1,2,3)]
        plus1=sp.Matrix([[2,sp.I],[-sp.I,3]])
        minus2=sp.I*sp.Matrix([[1,2],[2,-1]])
        plus3=sp.diag(3,5)
        overlaps=sp.symbols("w0:3")
        poles=list(zip(matrices,overlaps))
        vertex=cubic_vertex(*poles,plus1,minus2,plus3).doit()
        meson,alp=cubic_pair_source(*poles[:2],plus1,minus2,plus3)
        paired=(sp.Trace(meson*matrices[2])/2+alp*overlaps[2]).doit()
        self.assertEqual(sp.simplify(vertex-paired),0)

    def test_quartic_mass_polarization(self):
        labels=sp.symbols("x0:4")
        matrices=[sp.Matrix([[n,1],[1,-n]]) for n in (1,2,3,4)]
        overlaps=sp.symbols("w0:4")
        phi=sum((labels[i]*matrices[i] for i in range(4)),sp.zeros(2))
        alp=sum(labels[i]*overlaps[i] for i in range(4))
        derivatives=[sp.Matrix([[n+1,2],[2,-n]]) for n in range(5)]
        potential=(sp.Trace(derivatives[0]*phi**4)/(96*f_pi**2)-alp*sp.Trace(derivatives[1]*phi**3)/(24*sp.I*f_pi)
                   -alp**2*sp.Trace(derivatives[2]*phi**2)/16+f_pi*alp**3*sp.Trace(derivatives[3]*phi)/(24*sp.I)
                   +f_pi**2*alp**4*sp.Trace(derivatives[4])/96).doit()
        coefficient=sp.Poly(sp.expand(potential),*labels).coeff_monomial(sp.prod(labels))
        vertex=quartic_mass_vertex(list(zip(matrices,overlaps)),*derivatives).doit()
        self.assertEqual(sp.simplify(coefficient-vertex),0)

    def test_kinetic_trace_reduction(self):
        matrices=[sp.Matrix([[n,1],[1,-n]]) for n in (1,2,3,4)]
        poles=[(matrix,0) for matrix in matrices]
        mandelstam_s,mandelstam_t,mandelstam_u=sp.symbols("s t u")
        masses=sp.symbols("m0:4")
        scalar_products={(0,1):(mandelstam_s-masses[0]-masses[1])/2,(2,3):(mandelstam_s-masses[2]-masses[3])/2,
                         (0,2):(mandelstam_t-masses[0]-masses[2])/2,(1,3):(mandelstam_t-masses[1]-masses[3])/2,
                         (0,3):(mandelstam_u-masses[0]-masses[3])/2,(1,2):(mandelstam_u-masses[1]-masses[2])/2}
        self.assertEqual(sp.simplify(quartic_kinetic_vertex(poles,scalar_products).doit()
                         -quartic_kinetic_channel(poles,mandelstam_s,mandelstam_t,mandelstam_u).doit()),0)

    def test_neutral_kinetic_zero(self):
        incoming=sp.MatrixSymbol("A",3,3)
        pion=sp.MatrixSymbol("P",3,3)
        self.assertEqual(quartic_kinetic_channel([(incoming,1)]+[(pion,1)]*3,*sp.symbols("s t u")),0)

    def test_noncommuting_meson_inverse(self):
        masses=sp.Matrix([[7,1],[1,4]])
        source=sp.Matrix([[1,sp.I],[-sp.I,3]])
        anomaly=sp.Rational(3,2)
        invariant=sp.Rational(1,3)
        field=inverse_meson_mass_operator(source,masses,anomaly,invariant)
        recovered=B_0*anticommutator(masses,field)+anomaly*sp.trace(field)*sp.eye(2)/2-invariant*field
        self.assertEqual(sp.simplify(recovered-source),sp.zeros(2))

    def test_mixed_inverse(self):
        masses=sp.Matrix([[7,1],[1,4]])
        source=sp.Matrix([[1,sp.I],[-sp.I,3]])
        coupling=sp.Matrix([[1,2],[2,-1]])
        anomaly=sp.Rational(3,2)
        invariant=sp.Rational(1,3)
        alp_entry=sp.Rational(9,2)
        inverse=lambda matrix:inverse_meson_mass_operator(matrix,masses,anomaly,invariant)
        meson,alp=inverse_mixed_mass_operator((source,sp.Rational(2,3)),coupling,alp_entry,inverse,invariant)
        meson=meson.doit()
        alp=alp.doit()
        recovered=B_0*anticommutator(masses,meson)+anomaly*sp.trace(meson)*sp.eye(2)/2-invariant*meson-2*f_pi*coupling*alp/f_a
        self.assertEqual(sp.simplify(recovered-source),sp.zeros(2))
        self.assertEqual(sp.simplify((alp_entry-invariant)*alp-f_pi*sp.trace(coupling*meson)/f_a-sp.Rational(2,3)),0)


if __name__ == "__main__":
    if "--test" in sys.argv:
        unittest.main(argv=[sys.argv[0]])
    else:
        main()
