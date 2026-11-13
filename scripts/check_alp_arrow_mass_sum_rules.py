#!/usr/bin/env python3
"""Exact LO arrow-spectrum sum rules. Run: .venv/bin/python scripts/check_alp_arrow_mass_sum_rules.py

All dimensional inputs are in GeV and GeV squared; the SU(3)-symmetric
spurion makes the meson block diagonal and permits closed-form exact poles.
See memory/alp/shared/alp-lo-mass-sum-rules.md and alp-lo-arrow-matrix-spectrum.md.
"""
import sympy as sp
from derive_alp_lo_mass_sum_rules import (rational_sources, blocks_with_second_order_scalar,
    trace_sum_rule_closed_form, product_rule_closed_form)
from derive_alp_meson_mass_shifts_exact import canonical_mass_matrix_from_sources
from solve_alp_meson_spectrum import solve_spectrum_exact


def check_sum_rules():
    sources = rational_sources()
    sources.update(quark_spurion=sp.eye(3) * sp.Rational(2, 5),
                   axial_source=sp.eye(3) * sp.Rational(1, 7),
                   pseudoscalar_source=sp.eye(3) * sp.Rational(1, 5),
                   gluon_coupling=sp.Rational(1, 4),
                   anomaly_mass_squared=sp.Rational(3, 5),
                   bare_alp_mass_squared=sp.Rational(4, 5),
                   pion_decay_constant=sp.Rational(1, 10),
                   alp_decay_constant=sp.Rational(1, 2))
    second_order = dict(condensate=sp.Rational(2), second_order_scalar_trace=sp.Rational(1, 50))
    kinetic, mass = blocks_with_second_order_scalar(sources, **second_order)
    arrow = canonical_mass_matrix_from_sources(sources, **second_order)
    assert arrow[:5, :5].is_diagonal()
    assert any(arrow[i, 5] != 0 for i in range(5))
    arrow_poles = arrow.eigenvals()
    spectrum = solve_spectrum_exact(kinetic, mass)
    arrow_sum = sum(pole * multiplicity for pole, multiplicity in arrow_poles.items())
    arrow_product = sp.prod(pole ** multiplicity for pole, multiplicity in arrow_poles.items())
    solver_sum = sum(spectrum.pole_squared_masses)
    solver_product = sp.prod(spectrum.pole_squared_masses)
    expected_sum = trace_sum_rule_closed_form(**sources, **second_order)
    expected_product = product_rule_closed_form(**sources, **second_order)
    return {name: sp.simplify(value - expected) == 0 for name, value, expected in (
        ('arrow trace', arrow_sum, expected_sum), ('arrow product', arrow_product, expected_product),
        ('solver trace', solver_sum, expected_sum), ('solver product', solver_product, expected_product))}


if __name__ == '__main__':
    results = check_sum_rules()
    for label, passed in results.items():
        print(f"{label}: {'PASS' if passed else 'FAIL'}")
    if not all(results.values()):
        raise SystemExit(1)
