"""Compare Gauss-Legendre and adaptive 3P/VPP widths; all units are GeV.

Tests: .venv/bin/python -m unittest scripts.benchmark_three_body_quadrature
Timing: .venv/bin/python -m scripts.benchmark_three_body_quadrature --repeats 5
Illustrative sources and fitted hadronic inputs follow the scalar benchmark.
"""
import argparse
from collections import Counter
from math import factorial
import json
import math
import statistics
import time
import unittest
from unittest.mock import patch

import numpy as np
from scipy.integrate import quad
from alpaca_chiral import (
    ALPcouplings, three_pseudoscalar_amplitude, three_pseudoscalar_width,
    vector_two_pseudoscalar_amplitude,
    vector_two_pseudoscalar_width,
)
from alpaca_chiral.three_pseudoscalar import _pseudoscalars
from scripts.benchmark_three_body_scalar import cases as scalar_cases, coupling_inputs


def adaptive_reference_three_pseudoscalar_width(ma: float, fa: float, couplings: ALPcouplings,
                              meson_1: str, meson_2: str, meson_3: str) -> float:
    """Return the scalar-prepared Gamma(a -> P1 P2 P3) in GeV.

    Identical particles contribute 1/prod(n!). Numerical integration uses
    relative tolerance 1e-5. Thresholds use observed pseudoscalar masses.
    """
    names = meson_1, meson_2, meson_3
    particles = _pseudoscalars()
    first_mass, second_mass, third_mass = (particles[name][0] for name in names)
    if ma <= first_mass+second_mass+third_mass:
        return 0.0
    amplitude = three_pseudoscalar_amplitude(ma, fa, couplings, *names)

    first_squared, second_squared, third_squared = (
        mass*mass for mass in (first_mass, second_mass, third_mass))
    parent_squared = ma*ma
    pair_threshold = (first_mass+second_mass)**2
    pair_difference = (first_mass-second_mass)**2
    upper_invariant = (ma-third_mass)**2
    upper_sum = (ma+third_mass)**2

    def integrate_second_invariant(s12):
        """Integrate over the second Dalitz invariant at fixed s12."""
        center = second_squared+third_squared+(
            (s12+second_squared-first_squared)*(parent_squared-s12-third_squared))/(2*s12)
        product = ((s12-pair_threshold)*(s12-pair_difference)
                   *(upper_invariant-s12)*(upper_sum-s12))
        half_width = np.sqrt(max(0., product))/(2*s12)
        return quad(lambda s23: abs(amplitude(s12, s23))**2,
                    center-half_width, center+half_width, epsabs=0., epsrel=1e-5, limit=200)[0]

    integral = quad(integrate_second_invariant, (first_mass+second_mass)**2,
                    (ma-third_mass)**2, epsabs=0., epsrel=1e-5, limit=200)[0]
    symmetry = np.prod([factorial(number) for number in Counter(names).values()])
    return float(integral/(256*np.pi**3*ma**3*symmetry))


def adaptive_reference_vector_two_pseudoscalar_width(ma: float, fa: float, couplings: ALPcouplings,
                                  vector: str, meson_1: str, meson_2: str) -> float:
    """Return scalar-prepared Gamma(a -> VPP), in GeV.

    Uses relative tolerance 1e-5, internal pole subdivisions and the
    identical-pseudoscalar factor 1/2!. Fixed-s coefficients are reused
    by all evaluations of the inner quadrature.
    """
    amplitude = vector_two_pseudoscalar_amplitude(
        ma, fa, couplings, vector, meson_1, meson_2)
    if not amplitude.terms:
        return 0.
    vector_mass, first_mass, second_mass = amplitude.masses
    total_squared = amplitude.total_squared
    crossed_poles = tuple((invariant, mass*mass)
                          for _, invariant, mass, _, _, _ in amplitude.terms
                          if invariant in (1, 2))

    def integrate_t(s):
        """Integrate a prepared scalar polarization sum over t at fixed s."""
        lower, upper = amplitude.t_limits(s)
        points = set()
        for invariant, mass_squared in crossed_poles:
            pole = mass_squared if invariant == 1 else total_squared-s-mass_squared
            if lower < pole < upper:
                points.add(pole)
        return quad(amplitude.at_fixed_s(s, lower, upper), lower, upper,
                    points=sorted(points), epsabs=0., epsrel=1e-5, limit=150)[0]

    lower, upper = (first_mass+second_mass)**2, (ma-vector_mass)**2
    points = sorted({mass*mass for _, invariant, mass, _, _, _ in amplitude.terms
                     if invariant == 0 and lower < mass*mass < upper})
    integral = quad(integrate_t, lower, upper, points=points,
                    epsabs=0., epsrel=1e-5, limit=150)[0]
    identical_factor = 2 if amplitude.final_state[1] == amplitude.final_state[2] else 1
    return integral/(256*math.pi**3*ma**3*identical_factor)


def cases():
    """Physical thresholds, poles, kaons, repeated states and complex ds sources."""
    diagonal = coupling_inputs()[0]
    source = np.zeros((3, 3), dtype=complex)
    source[0, 1], source[1, 0] = .2+.1j, .2-.1j
    # Down-sector matrix indices 0,1 map to physical d,s.
    zero = np.zeros((3, 3))
    changing = ALPcouplings(.1, np.diag([.8, 0., 0.]), zero, source, source)
    return scalar_cases()+[
        ('3P', .415, diagonal, ('pi0', 'pi+', 'pi-')),
        ('3P', 1.5, changing, ('pi+', 'pi-', 'KS')),
        ('VPP', 1.063, diagonal, ('omega', 'pi+', 'pi-')),
    ]


def functions(family):
    """Retained adaptive benchmark reference and official Gauss-Legendre width."""
    return ((adaptive_reference_three_pseudoscalar_width, three_pseudoscalar_width)
            if family == '3P' else
            (adaptive_reference_vector_two_pseudoscalar_width, vector_two_pseudoscalar_width))


def tightened_reference(family, mass, scale, couplings, *state):
    """Independent scalar nested integration, pole splits, tolerance 2e-7."""
    if family == '3P':
        amplitude = three_pseudoscalar_amplitude(mass, scale, couplings, *state)
        first, second, spectator = (_pseudoscalars()[name][0] for name in state)
        poles = amplitude.invariant_poles
        total = mass**2+first**2+second**2+spectator**2
        def limits(s):
            center = second**2+spectator**2+(s+second**2-first**2)*(mass**2-s-spectator**2)/(2*s)
            product = (s-(first+second)**2)*(s-(first-second)**2)
            product *= ((mass-spectator)**2-s)*((mass+spectator)**2-s)
            half = math.sqrt(max(0., product))/(2*s)
            return center-half, center+half
        squared = lambda s, t: abs(amplitude(s, t))**2
        symmetry = math.prod(math.factorial(number) for number in Counter(state).values())
    else:
        amplitude = vector_two_pseudoscalar_amplitude(mass, scale, couplings, *state)
        if not amplitude.terms:
            return 0.
        spectator, first, second = amplitude.masses
        total = amplitude.total_squared
        poles = tuple(tuple(pole for _, pole, _, _ in group) for group in amplitude.prepared_terms)
        limits = amplitude.t_limits
        squared = amplitude.squared
        symmetry = 2 if state[1] == state[2] else 1
    lower, upper = (first+second)**2, (mass-spectator)**2
    def integrate_inner(s):
        inner_lower, inner_upper = limits(s)
        points = sorted({pole.real for pole in poles[1] if inner_lower < pole.real < inner_upper}
                        |{total-s-pole.real for pole in poles[2] if inner_lower < total-s-pole.real < inner_upper})
        return quad(lambda t: squared(s,t), inner_lower, inner_upper, points=points,
                    epsabs=0., epsrel=2e-7, limit=250)[0]
    points = sorted({pole.real for pole in poles[0] if lower < pole.real < upper})
    integral = quad(integrate_inner, lower, upper, points=points,
                    epsabs=0., epsrel=2e-7, limit=250)[0]
    return integral/(256*math.pi**3*mass**3*symmetry)


class ThreeBodyQuadratureTests(unittest.TestCase):
    """Complete widths, independent tightened integration, order and zeros."""
    def test_widths_and_order_convergence(self):
        for family, mass, couplings, state in cases():
            old, new = functions(family)
            arguments = mass, 100., couplings, *state
            with self.subTest(family=family, mass=mass, state=state):
                original = old(*arguments)
                candidate = new(*arguments)
                tightened = tightened_reference(family, *arguments)
                refined = new(*arguments, order=32, epsrel=2e-6)
                np.testing.assert_allclose(candidate, original, rtol=3e-5, atol=0.)
                np.testing.assert_allclose(candidate, tightened, rtol=1e-5, atol=0.)
                np.testing.assert_allclose(candidate, refined, rtol=1e-5, atol=0.)

    def test_zeros(self):
        diagonal = coupling_inputs()[0]
        zero = np.zeros((3,3))
        absent = ALPcouplings(0.,zero,zero,zero,zero)
        for family,mass,couplings,state in [
                ('3P',.3,diagonal,('pi0','pi+','pi-')),
                ('3P',1.3,absent,('pi0','pi+','pi-')),
                ('VPP',1.,diagonal,('omega','pi+','pi-')),
                ('VPP',1.8,diagonal,('omega','pi0','pi0')),
                ('VPP',1.8,absent,('omega','pi+','pi-'))]:
            for function in functions(family):
                self.assertEqual(function(mass,100.,couplings,*state),0.)

    def test_exchange(self):
        diagonal = coupling_inputs()[0]
        for family,state,exchanged in [
                ('3P',('pi0','pi+','pi-'),('pi+','pi0','pi-')),
                ('VPP',('omega','pi+','pi-'),('omega','pi-','pi+'))]:
            new = functions(family)[1]
            np.testing.assert_allclose(new(1.5,100.,diagonal,*state),
                                       new(1.5,100.,diagonal,*exchanged),rtol=1e-5,atol=0.)

    def test_fallback(self):
        for family,state in [('3P',('pi0','pi+','pi-')),('VPP',('omega','pi+','pi-'))]:
            arguments = 1.5,100.,coupling_inputs()[0],*state
            with patch('alpaca_chiral._three_body_quadrature.quad',wraps=quad) as adaptive:
                value = functions(family)[1](*arguments,order=1,epsrel=1e-8)
            self.assertGreater(adaptive.call_count,1)
            np.testing.assert_allclose(value,tightened_reference(family,*arguments),rtol=1e-7,atol=0.)


def timed(function, arguments, repeats):
    """Median full-call timing, including amplitude and quadrature preparation."""
    durations=[]
    for _ in range(repeats):
        start=time.perf_counter()
        width=function(*arguments)
        durations.append(time.perf_counter()-start)
    return width,statistics.median(durations)


def benchmark(repeats=5):
    """Machine-specific full-call medians and width discrepancies as JSON."""
    records=[]
    for family,mass,couplings,state in cases():
        old,new=functions(family)
        arguments=mass,100.,couplings,*state
        old(*arguments);new(*arguments)
        original,old_seconds=timed(old,arguments,repeats)
        candidate,new_seconds=timed(new,arguments,repeats)
        records.append(dict(family=family,ma_GeV=mass,final_state=state,width_GeV=original,
            relative_difference=abs(candidate/original-1) if original else 0.,
            adaptive_seconds=old_seconds,gauss_legendre_seconds=new_seconds,
            speedup=old_seconds/new_seconds))
    return dict(repeats=repeats,initial_order=16,epsrel=1e-5,results=records)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repeats',type=int,default=5)
    arguments=parser.parse_args()
    print(json.dumps(benchmark(arguments.repeats),indent=2))
