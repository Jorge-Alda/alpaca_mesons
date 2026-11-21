"""Derive aligned mass derivatives; run .venv/bin/python scripts/derive_aligned_mass_anticommutators.py.

Dimensionless x=a/f_a; mass matrices have dimension one; no expansion in f_pi/f_a.
"""
import sympy as sp


def scalar_pseudoscalar_derivatives():
    x = sp.Symbol('x', real=True)
    anticommutator_map = sp.Symbol(r'\mathcal A', real=True)
    scalar = sp.symbols(r'\widehat S_0:5', real=True)
    pseudoscalar = sp.symbols(r'\widehat P_0:5', real=True)
    scalar_source = sum(x**n * scalar[n] / sp.factorial(n) for n in range(5))
    pseudoscalar_source = sum(x**n * pseudoscalar[n] / sp.factorial(n) for n in range(5))
    aligned = sp.series(sp.exp(sp.I*x*anticommutator_map)*(scalar_source+sp.I*pseudoscalar_source), x, 0, 5).removeO().expand(complex=True)
    return anticommutator_map, scalar, pseudoscalar, [
        (sp.expand(sp.diff(sp.re(aligned), x, n).subs(x, 0)),
         sp.expand(sp.diff(sp.im(aligned), x, n).subs(x, 0))) for n in range(5)]


def derive():
    x = sp.Symbol('x', real=True)
    axial = sp.Symbol(r'\widehat c_q^A', commutative=False)
    mass = sp.symbols(r'\widehat M_0:5', commutative=False)
    exponential = sum((sp.I*x*axial)**k / sp.factorial(k) for k in range(5))
    source = sum(x**k * mass[k] / sp.factorial(k) for k in range(5))
    aligned = sp.expand(exponential * source * exponential)
    derivatives = [sp.expand(sp.diff(aligned, x, n).subs(x, 0)) for n in range(5)]
    return axial, mass, derivatives


def nested_anticommutator(axial, matrix, order):
    for _ in range(order):
        matrix = sp.expand(axial*matrix + matrix*axial)
    return matrix


def shifted_derivatives():
    x = sp.Symbol('x', real=True)
    axial = sp.Symbol(r'\widehat c_q^A', commutative=False)
    mass = sp.symbols(r'\widehat M_0:5', commutative=False)
    source = sum(x**n * mass[n] / sp.factorial(n) for n in range(5))
    results = [source]
    for _ in range(4):
        previous = results[-1]
        results.append(sp.expand(sp.diff(previous, x) + sp.I * (axial * previous + previous * axial)))
    return [sp.expand(result.subs(x, 0)) for result in results]


def main():
    _, _, _, scalar_pseudoscalar = scalar_pseudoscalar_derivatives()
    for n, (scalar, pseudoscalar) in enumerate(scalar_pseudoscalar):
        print(f'Scalar derivative order {n}:\n{sp.latex(scalar)}\n')
        print(f'Pseudoscalar derivative order {n}:\n{sp.latex(pseudoscalar)}\n')
    axial, mass, derivatives = derive()
    for n, result in enumerate(shifted_derivatives()):
        print(f'Derivative plus anticommutator, order {n}:\n{sp.latex(result)}\n')
        assert sp.expand(result - derivatives[n]) == 0
    print('Aligned source: exp(i x axial) mass(x) exp(i x axial); x=a/f_a')
    for n, derivative in enumerate(derivatives):
        print(f'Order {n}:\n{sp.latex(derivative)}\n')
        expected = sum(sp.binomial(n,k)*sp.I**k*nested_anticommutator(axial,mass[n-k],k) for k in range(n+1))
        assert sp.expand(derivative-expected) == 0
    x = sp.Symbol('x')
    anticommutator_map = sp.Symbol(r'\mathcal A_{\widehat c_q^A}')
    mass_coefficients = sp.symbols(r'\widehat M_0:5')
    source = sum(x**n * mass_coefficients[n] / sp.factorial(n) for n in range(5))
    generating_series = sp.series(sp.exp(sp.I*x*anticommutator_map)*source, x, 0, 5).removeO().expand()
    print('Formal operator exponential times the mass-source series:')
    print(sp.latex(generating_series))
    for n in range(5):
        coefficient = sp.expand(sp.diff(generating_series, x, n).subs(x, 0))
        print(f'Formal derivative order {n}:\n{sp.latex(coefficient)}\n')
        assert coefficient == sp.expand(sum(sp.binomial(n,k)*sp.I**k*anticommutator_map**k*mass_coefficients[n-k] for k in range(n+1)))


if __name__ == '__main__':
    main()
