#!/usr/bin/env python3
"""Loop-to-tree crossing values f_pi/f_a for ALP lines with physical inputs.

Recomputes the regime map (section L4) of memory/alp/loop-corrections/finite-mixing/power-counting-and-validity/loop-counting/alp-loop-power-counting.md
with physical inputs, keeping or dropping the MSbar logarithms.

Conventions (Alda:2026nrd; dressed frame and unmixed approximation of the power-counting
record, imported from derive_alp_loop_power_counting.py):
  * Units: GeV.  r = f_pi/f_a.  B_0 below is always the low-energy constant
    B_0 = m_pi0^2/(m_u + m_d) (Gasser:1983yg); the dimensionless MSbar bubble function is
    called `bubble_function` in the code and "the bubble function" in words.
  * Meson-line masses from quark masses (not from the physical kaon mass):
        m_qq^2 = 2 B_0 m_q,   m_qr^2 = B_0 (m_q + m_r).
  * Tadpole function  A_0(m^2) = m^2/(16 pi^2) [1 - ln(m^2/mu^2)],  mu = 0.77 GeV.
    Size estimate: |A_0(m^2)| -> m^2/(16 pi^2).
  * Bubble function  B(s; m_1, m_2) = -int_0^1 dx ln[(x m_1^2 + (1-x) m_2^2 - x(1-x) s - i0)/mu^2]
    (MSbar, dimensionless).  Evaluation point (choice): the threshold of the intermediate
    state, s = (m_1 + m_2)^2, where Im B = 0 and
        B_threshold = 2 - 2 [m_1 ln(m_1/mu) + m_2 ln(m_2/mu)]/(m_1 + m_2).
    Size estimate: |B| -> 1.
  * Source convention: c_q^P = 0 and either c_q^A = c * identity ("universal") or
    c_u^A = c, c_d^A = c_s^A = 0 ("up_only").  kappa_q = c_q^A (c_q^A - 2 c_q^P),
    kappa4_q = (c_q^A)^3 (c_q^A - 4 c_q^P), c_q^AP = c_q^A - c_q^P.
  * Loop classes (ratio |loop/tree| = 1 solved for r):
      (i)   ALP line in the meson entry qq at fixed (tuned) m_a:
            ratio = 2 kappa_q r^2 |A_0(m_a^2)|/f_pi^2,  independent of B_0;
      (ii)  aP bubble in P_q P_q -> P_q P_q:
            ratio = 4 (c_q^AP)^2 r^2 m_qq^2 |B(s; m_a, m_qq)|/(16 pi^2 f_pi^2);
      (iii) aa bubble in P_q P_q -> P_q P_q:
            ratio = 4 kappa_q^2 r^4 m_qq^2 |B(s; m_a, m_a)|/(16 pi^2 f_pi^2);
      (iv)  ALP line in M^2_aa relative to m_a^2:
            ratio = V_aaaa |A_0(m_a^2)|/(2 m_a^2)
                  = 4 r^4 sum_q kappa4_q m_qq^2 |1 - ln(m_a^2/mu^2)|/(16 pi^2 f_pi^2);
      (v)   meson-only loops: ratio m_P^2/(16 pi^2 f_pi^2) (times |1 - ln(m^2/mu^2)| for the
            tadpole or |B_threshold| for the PP bubble), independent of r: no crossing.
    Untuned comparison of (i): m_a^2 = M^2_aa|quark = 2 r^2 sum_q kappa_q m_qq^2.

Minimal invocation:
    .venv/bin/python scripts/derive_alp_loop_crossing_table.py
Tests:
    .venv/bin/python -m unittest discover -s scripts -p test_derive_alp_loop_crossing_table.py
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

from derive_alp_loop_power_counting import loop_function_A0

FLAVOURS = ('u', 'd', 's')
PHYSICAL_INPUTS = {
    'neutral_pion_mass': 0.13498,
    'pion_decay_constant': 0.09207,
    'charged_kaon_mass': 0.49368,   # only for the comparison 2 m_K^2 - m_pi^2
    'renormalization_scale': 0.77,
}
QUARK_MASS_SETS = {
    'primary': (0.0022, 0.0047, 0.095),      # charge_model_parameters default
    'alternative': (0.00216, 0.00470, 0.095),  # alp-lo-partial-wave-unitarity.md (m_s added)
}
ALP_MASSES = (0.1, 0.3, 0.5, 1.0, 2.0)
COUPLING_VALUES = (1, 3)
LOOP_FACTOR = 16*np.pi**2


def source_coefficients(convention, coupling):
    """(c_q^A, c_q^P) for q = u, d, s in the two source conventions (c_q^P = 0)."""
    if convention == 'universal':
        return (coupling,)*3, (0.0,)*3
    if convention == 'up_only':
        return (coupling, 0.0, 0.0), (0.0,)*3
    raise ValueError(convention)


def low_energy_constant_B0(quark_masses, pion_mass=PHYSICAL_INPUTS['neutral_pion_mass']):
    """B_0 = m_pi0^2/(m_u + m_d) (GeV)."""
    return pion_mass**2/(quark_masses[0]+quark_masses[1])


def meson_line_squared_masses(quark_masses, low_energy_B0):
    """m_qr^2 = B_0 (m_q + m_r) for all flavour pairs (diagonal: 2 B_0 m_q)."""
    return {(q, r): low_energy_B0*(quark_masses[q]+quark_masses[r])
            for q in range(3) for r in range(3)}


def tadpole_size(squared_mass, scale, keep_logarithm):
    """|A_0(m^2)| with the MSbar logarithm, or the size estimate m^2/(16 pi^2)."""
    if keep_logarithm:
        return abs(loop_function_A0(float(squared_mass), scale**2))
    return squared_mass/LOOP_FACTOR


def bubble_function_at_threshold(first_mass, second_mass, scale):
    """Real MSbar bubble function at s = (m_1 + m_2)^2 (closed form, Im = 0 there)."""
    return 2-2*(first_mass*np.log(first_mass/scale)
                + second_mass*np.log(second_mass/scale))/(first_mass+second_mass)


def bubble_size(first_mass, second_mass, scale, keep_logarithm):
    """|bubble function| at threshold, or the size estimate 1."""
    if keep_logarithm:
        return abs(bubble_function_at_threshold(first_mass, second_mass, scale))
    return 1.0


# ------------------------------------------------------------------ ratios (loop/tree)

def ratio_alp_line_in_meson_entry(ratio_f, axial, pseudoscalar, alp_mass, f_pi, scale,
                                  keep_logarithm):
    """(i) 2 kappa_q r^2 |A_0(m_a^2)|/f_pi^2 for a diagonal entry."""
    kappa = axial*(axial-2*pseudoscalar)
    return abs(2*kappa*ratio_f**2*tadpole_size(alp_mass**2, scale, keep_logarithm))/f_pi**2


def ratio_alp_meson_bubble(ratio_f, axial, pseudoscalar, meson_squared_mass, alp_mass, f_pi,
                           scale, keep_logarithm):
    """(ii) 4 (c^AP)^2 r^2 m_qq^2 |B|/(16 pi^2 f_pi^2), B at s = (m_a + m_qq)^2."""
    bubble = bubble_size(alp_mass, np.sqrt(meson_squared_mass), scale, keep_logarithm)
    return 4*(axial-pseudoscalar)**2*ratio_f**2*meson_squared_mass*bubble/(LOOP_FACTOR*f_pi**2)


def ratio_alp_pair_bubble(ratio_f, axial, pseudoscalar, meson_squared_mass, alp_mass, f_pi,
                          scale, keep_logarithm):
    """(iii) 4 kappa^2 r^4 m_qq^2 |B|/(16 pi^2 f_pi^2), B at s = 4 m_a^2."""
    kappa = axial*(axial-2*pseudoscalar)
    bubble = bubble_size(alp_mass, alp_mass, scale, keep_logarithm)
    return 4*kappa**2*ratio_f**4*meson_squared_mass*bubble/(LOOP_FACTOR*f_pi**2)


def ratio_alp_line_in_alp_entry(ratio_f, axial_values, pseudoscalar_values, meson_squared_masses,
                                alp_mass, f_pi, scale, keep_logarithm):
    """(iv) V_aaaa |A_0(m_a^2)|/(2 m_a^2), V_aaaa = 16 f_pi^2 B_0 sum m_q kappa4_q/f_a^4."""
    weighted = sum(a**3*(a-4*p)*meson_squared_masses[(q, q)]
                   for q, (a, p) in enumerate(zip(axial_values, pseudoscalar_values)))
    alp_vertex = 8*weighted*ratio_f**4/f_pi**2  # 16 f^2 B_0 sum m_q kappa4/f_a^4
    return abs(alp_vertex*tadpole_size(alp_mass**2, scale, keep_logarithm)/(2*alp_mass**2))


def untuned_alp_squared_mass(ratio_f, axial_values, pseudoscalar_values, meson_squared_masses):
    """M^2_aa|quark = 2 r^2 sum_q kappa_q m_qq^2 (no bare or topological ALP mass)."""
    return 2*ratio_f**2*sum(a*(a-2*p)*meson_squared_masses[(q, q)]
                            for q, (a, p) in enumerate(zip(axial_values, pseudoscalar_values)))


def ratio_untuned_alp_line_in_meson_entry(ratio_f, flavour, axial_values, pseudoscalar_values,
                                          meson_squared_masses, f_pi, scale, keep_logarithm):
    """(i) with m_a^2 = M^2_aa|quark(r)."""
    alp_squared_mass = untuned_alp_squared_mass(ratio_f, axial_values, pseudoscalar_values,
                                                meson_squared_masses)
    return ratio_alp_line_in_meson_entry(ratio_f, axial_values[flavour],
                                         pseudoscalar_values[flavour], np.sqrt(alp_squared_mass),
                                         f_pi, scale, keep_logarithm)


def meson_loop_ratios(meson_squared_mass, f_pi, scale):
    """(v) meson-only loop/tree: size estimate, tadpole with logarithm, PP bubble at threshold."""
    size = meson_squared_mass/(LOOP_FACTOR*f_pi**2)
    return {'size_estimate': size,
            'tadpole_with_logarithm': size*abs(1-np.log(meson_squared_mass/scale**2)),
            'bubble_at_threshold': size*abs(bubble_function_at_threshold(
                np.sqrt(meson_squared_mass), np.sqrt(meson_squared_mass), scale))}


# ------------------------------------------------------------------ crossings (closed forms)

def crossing_from_ratio(ratio_at_unit_r, power):
    """Solve ratio_at_unit_r * r^power = 1; infinity when the coefficient vanishes."""
    if ratio_at_unit_r <= 0:
        return np.inf
    return ratio_at_unit_r**(-1/power)


def crossing_alp_line_in_meson_entry(axial, pseudoscalar, alp_mass, f_pi, scale, keep_logarithm):
    """(i) r_* = f_pi/sqrt(2 kappa |A_0|); size estimate 4 pi f_pi/(sqrt(2 kappa) m_a)."""
    return crossing_from_ratio(ratio_alp_line_in_meson_entry(
        1.0, axial, pseudoscalar, alp_mass, f_pi, scale, keep_logarithm), 2)


def crossing_alp_meson_bubble(axial, pseudoscalar, meson_squared_mass, alp_mass, f_pi, scale,
                              keep_logarithm):
    """(ii) r_* = 4 pi f_pi/(2 |c^AP| m_qq sqrt|B|)."""
    return crossing_from_ratio(ratio_alp_meson_bubble(
        1.0, axial, pseudoscalar, meson_squared_mass, alp_mass, f_pi, scale, keep_logarithm), 2)


def crossing_alp_pair_bubble(axial, pseudoscalar, meson_squared_mass, alp_mass, f_pi, scale,
                             keep_logarithm):
    """(iii) r_* = [16 pi^2 f_pi^2/(4 kappa^2 m_qq^2 |B|)]^(1/4)."""
    return crossing_from_ratio(ratio_alp_pair_bubble(
        1.0, axial, pseudoscalar, meson_squared_mass, alp_mass, f_pi, scale, keep_logarithm), 4)


def crossing_alp_line_in_alp_entry(axial_values, pseudoscalar_values, meson_squared_masses,
                                   alp_mass, f_pi, scale, keep_logarithm):
    """(iv) r_* = [16 pi^2 f_pi^2/(4 sum kappa4 m_qq^2 |1 - ln(m_a^2/mu^2)|)]^(1/4)."""
    return crossing_from_ratio(ratio_alp_line_in_alp_entry(
        1.0, axial_values, pseudoscalar_values, meson_squared_masses, alp_mass, f_pi, scale,
        keep_logarithm), 4)


def crossing_untuned_alp_line_in_meson_entry(flavour, axial_values, pseudoscalar_values,
                                             meson_squared_masses, f_pi, scale, keep_logarithm,
                                             largest_ratio=30.0):
    """Untuned (i).  Size estimate: r_* = [16 pi^2 f_pi^2/(4 kappa_q sum kappa m^2)]^(1/4).
    With the logarithm: smallest root of ratio = 1 on (0, largest_ratio) (grid + brentq)."""
    def excess(ratio_f):
        return ratio_untuned_alp_line_in_meson_entry(
            ratio_f, flavour, axial_values, pseudoscalar_values, meson_squared_masses, f_pi,
            scale, keep_logarithm)-1
    if not keep_logarithm:
        return crossing_from_ratio(excess(1.0)+1, 4)
    grid = np.linspace(1e-3, largest_ratio, 30000)
    values = [excess(point) for point in grid]
    for left, right, left_value, right_value in zip(grid[:-1], grid[1:], values[:-1], values[1:]):
        if left_value < 0 <= right_value:
            return brentq(excess, left, right, xtol=1e-12)
    return np.inf


# ------------------------------------------------------------------ table

def crossing_table(quark_masses, convention, coupling, keep_logarithm,
                   scale=PHYSICAL_INPUTS['renormalization_scale'], alp_masses=ALP_MASSES):
    """Rows: class label -> list of crossings over alp_masses."""
    f_pi = PHYSICAL_INPUTS['pion_decay_constant']
    low_energy_B0 = low_energy_constant_B0(quark_masses)
    squared = meson_line_squared_masses(quark_masses, low_energy_B0)
    axial, pseudoscalar = source_coefficients(convention, coupling)
    rows = {}
    for q, label in ((0, 'u'), (2, 's')):
        rows[f'(i) ALP line, {label}{label} entry'] = [crossing_alp_line_in_meson_entry(
            axial[q], pseudoscalar[q], m_a, f_pi, scale, keep_logarithm) for m_a in alp_masses]
    for q, label in ((0, 'u'), (2, 's')):
        rows[f'(ii) aP bubble, {label}'] = [crossing_alp_meson_bubble(
            axial[q], pseudoscalar[q], squared[(q, q)], m_a, f_pi, scale, keep_logarithm)
            for m_a in alp_masses]
    for q, label in ((0, 'u'), (2, 's')):
        rows[f'(iii) aa bubble, {label}'] = [crossing_alp_pair_bubble(
            axial[q], pseudoscalar[q], squared[(q, q)], m_a, f_pi, scale, keep_logarithm)
            for m_a in alp_masses]
    rows['(iv) ALP line in M2_aa'] = [crossing_alp_line_in_alp_entry(
        axial, pseudoscalar, squared, m_a, f_pi, scale, keep_logarithm) for m_a in alp_masses]
    return rows


def format_table(rows, alp_masses=ALP_MASSES):
    header = f'{"class":28s}' + ''.join(f'{m:>9.1f}' for m in alp_masses)
    lines = [header]
    for label, values in rows.items():
        lines.append(f'{label:28s}' + ''.join(
            f'{"none":>9s}' if not np.isfinite(v) else f'{v:9.3f}' for v in values))
    return '\n'.join(lines)


def main():
    f_pi = PHYSICAL_INPUTS['pion_decay_constant']
    scale = PHYSICAL_INPUTS['renormalization_scale']
    for name, masses in QUARK_MASS_SETS.items():
        low_energy_B0 = low_energy_constant_B0(masses)
        squared = meson_line_squared_masses(masses, low_energy_B0)
        print(f'quark masses {name} {masses} GeV: B_0 = {low_energy_B0:.4f} GeV,'
              f' m_uu^2 = {squared[(0, 0)]:.5f}, m_dd^2 = {squared[(1, 1)]:.5f},'
              f' m_ss^2 = {squared[(2, 2)]:.5f} GeV^2')
    kaon_based = 2*PHYSICAL_INPUTS['charged_kaon_mass']**2-PHYSICAL_INPUTS['neutral_pion_mass']**2
    print(f'2 m_K^2 - m_pi0^2 = {kaon_based:.5f} GeV^2 (comparison for m_ss^2)')
    primary = QUARK_MASS_SETS['primary']
    alternative = QUARK_MASS_SETS['alternative']
    print(f'A_0 vanishes at m_a = sqrt(e) mu = {np.sqrt(np.e)*scale:.4f} GeV;'
          f' B(4 m_a^2) vanishes at m_a = e mu = {np.e*scale:.4f} GeV')
    for convention in ('universal', 'up_only'):
        for coupling in COUPLING_VALUES:
            for keep_logarithm in (False, True):
                tag = 'MSbar logarithms kept' if keep_logarithm else 'size estimate'
                print(f'\n== crossing f_pi/f_a, convention {convention}, c = {coupling}, {tag},'
                      f' mu = {scale} GeV; columns m_a/GeV')
                rows = crossing_table(primary, convention, coupling, keep_logarithm)
                print(format_table(rows))
                changed = crossing_table(alternative, convention, coupling, keep_logarithm)
                shifts = {label: [(b/a-1) if np.isfinite(a) and np.isfinite(b) else 0.0
                                  for a, b in zip(rows[label], changed[label])] for label in rows}
                print('relative change with (2.16, 4.70, 95) MeV: max |shift| per class',
                      {k: f'{max(abs(x) for x in v):.4f}' for k, v in shifts.items()})
    print('\n== untuned case (i): m_a^2 = M^2_aa|quark, u entry')
    for convention in ('universal', 'up_only'):
        for coupling in COUPLING_VALUES:
            axial, pseudoscalar = source_coefficients(convention, coupling)
            for name, masses in QUARK_MASS_SETS.items():
                squared = meson_line_squared_masses(masses, low_energy_constant_B0(masses))
                size = crossing_untuned_alp_line_in_meson_entry(
                    0, axial, pseudoscalar, squared, f_pi, scale, False)
                full = crossing_untuned_alp_line_in_meson_entry(
                    0, axial, pseudoscalar, squared, f_pi, scale, True)
                alp_mass = np.sqrt(untuned_alp_squared_mass(size, axial, pseudoscalar, squared))
                print(f'{convention:9s} c={coupling} {name:11s}: size estimate {size:.3f}'
                      f' (m_a = {alp_mass:.3f} GeV), with logarithm {full:.3f}')
    print('\n== (v) meson-only loops (no crossing, independent of f_pi/f_a), primary set')
    squared = meson_line_squared_masses(primary, low_energy_constant_B0(primary))
    for label, value in (('uu', squared[(0, 0)]), ('dd', squared[(1, 1)]),
                         ('ss', squared[(2, 2)]),
                         ('pi0', PHYSICAL_INPUTS['neutral_pion_mass']**2)):
        print(label, {k: round(v, 4) for k, v in meson_loop_ratios(value, f_pi, scale).items()})


if __name__ == '__main__':
    main()
