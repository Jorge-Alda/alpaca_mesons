"""Switchable charged/neutral K*(892) treatment of the joint vector likelihood.

Run: .venv/bin/python scripts/joint_vector_kstar_charge_options.py
(prints coordinate names, starting values and bounds for every option setting).

The baseline likelihood in fit_babar_kkpi_global_vector_sector.py and the
earlier variations in joint_vector_likelihood_options.py are left
byte-identical, because the saved fit ledgers verify their SHA-256. The class
below inherits from VariedLikelihood; with the three options below at their
defaults it calls the inherited methods unchanged, so residuals and
Jacobians are bitwise equal to the baseline.

Baseline K*(892) treatment (all options off):
  * one K* mass m_K* (the mass-spectrum value; its likelihood row is the PDG
    charged K* mass) in the vertex normalization and in every BaBar K*
    propagator;
  * K*+ and K*0 total widths predicted by K pi saturation,
    Gamma = G_V^2 * sum over the two K pi modes of fit_vector_gv.width_coefficient
    (physical parent and daughter masses in phase space, m_K* in the vertex);
  * every BaBar K* propagator uses the mean of the two predicted widths.

Options:
  separate_kstar_charge_states (default False): the BaBar K* propagators use
      charge-specific masses and widths. Channel assignment, read off the
      phase-space construction of phase_points: in K+ K- pi0 both t and u are
      K(+-) pi0 invariants (charged K*); in K_S K(+-) pi(-+) the integrated
      invariant t = (p_K+- + p_pi)^2 is the neutral K*, and
      u = q + 2 MK^2 + MPI^2 - s - t the charged K*. The charged K* has mass
      m_K* and width Gamma(K*+); the neutral K* has mass m_K* + Delta m and
      width Gamma(K*0). In the K pi saturation prediction of Gamma(K*0) the
      vertex mass becomes m_K* + Delta m as well. Delta m is fixed at the PDG
      2024 difference of the physical masses already stored in the two K*
      width rows (0.89555 - 0.89167 = 0.00388 GeV), unless the next option
      makes it a coordinate.
  free_kstar_mass_difference (default False; requires the previous option):
      Delta m = m(K*0) - m(K*+) in GeV becomes a fitted coordinate.
      Start: the PDG 2024 difference above. Bounds: [-0.02, 0.02] GeV
      (assumption: five times the PDG difference; no likelihood row measures
      the K*0 mass, so Delta m is constrained by the BaBar bins only).
  kstar_width_treatment (default 'kpi_saturation'): 'free_total_widths'
      replaces the K pi saturation prediction of the K*+ and K*0 total widths
      by two total-width pulls, Gamma = value + sigma * pull, exactly as the
      baseline treats rho0, rho+, omega, phi and eta'. The widths then no
      longer depend on G_V. Start: pull 0 (the PDG 2024 measured widths of
      the likelihood rows). Bounds: lower -0.99 value/sigma, no upper bound,
      as for the baseline pulls.

Added coordinates are appended after the 32 baseline coordinates, in the
order of ADDED_COORDINATE_NAMES restricted to the selected options, so the
baseline indices are unchanged. Units: GeV for masses, widths and energies;
GeV^2 for invariants; residuals are dimensionless. Conventions as in
Alda:2026nrd.
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fit_babar_kkpi_global_vector_sector as joint_fit
import fit_vector_gv as gv_fit
import joint_vector_likelihood_options as options_module
from fit_babar_kkpi_g_coefficients import CHANNELS

BASELINE_COORDINATE_COUNT = 32
KSTAR_WIDTH_TREATMENTS = ('kpi_saturation', 'free_total_widths')
KSTAR_NAMES = ('Kstarplus', 'Kstar0')
ADDED_COORDINATE_NAMES = ('total_width_pull_Kstarplus', 'total_width_pull_Kstar0',
                          'kstar_neutral_minus_charged_mass_GeV')
MASS_DIFFERENCE_BOUNDS_GEV = (-.02, .02)
RUNNING_PROPAGATOR = joint_fit.running_propagator
# K pi modes of each K* charge state: (first daughter, second daughter, flavor factor key).
KSTAR_MODES = {
    'Kstarplus': (('K_charged', 'pi0', 'kstar_charged_K+_pi0'),
                  ('K_neutral', 'pi_charged', 'kstar_charged_K0_pi+')),
    'Kstar0': (('K_neutral', 'pi0', 'kstar_neutral_K0_pi0'),
               ('K_charged', 'pi_charged', 'kstar_neutral_K+_pi-'))}


class KstarChargeLikelihood(options_module.VariedLikelihood):
    """VariedLikelihood with switchable charged/neutral K* masses and widths.

    See the module docstring for the definition, starting values and bounds
    of each option. All VariedLikelihood options remain available.
    """

    def __init__(self, energy_order=2, pair_order=16, angle_order=16, bin_block_size=None,
                 separate_kstar_charge_states=False, free_kstar_mass_difference=False,
                 kstar_width_treatment='kpi_saturation', **variation_options):
        super().__init__(energy_order, pair_order, angle_order, bin_block_size, **variation_options)
        if kstar_width_treatment not in KSTAR_WIDTH_TREATMENTS:
            raise ValueError(f'unknown K* width treatment {kstar_width_treatment!r}')
        if free_kstar_mass_difference and not separate_kstar_charge_states:
            raise ValueError('a free K* mass difference requires separate K* charge states')
        self.kstar_options = dict(separate_kstar_charge_states=separate_kstar_charge_states,
                                  free_kstar_mass_difference=free_kstar_mass_difference,
                                  kstar_width_treatment=kstar_width_treatment)
        self.free_kstar_widths = kstar_width_treatment == 'free_total_widths'
        rows = dict(zip(KSTAR_NAMES, self.gv_data['rows'][:2]))
        self.pdg_kstar_mass_difference_gev = (rows['Kstar0']['physical_vector_mass_GeV']
                                              - rows['Kstarplus']['physical_vector_mass_GeV'])
        names = []
        if self.free_kstar_widths:
            names += list(ADDED_COORDINATE_NAMES[:2])
        if free_kstar_mass_difference:
            names.append(ADDED_COORDINATE_NAMES[2])
        self.added_coordinate_names = tuple(names)
        self.added_index = {name: BASELINE_COORDINATE_COUNT+index for index, name in enumerate(names)}
        self.kstar_default = not (separate_kstar_charge_states or self.free_kstar_widths)
        self.last_added_coordinates = None

    @property
    def parameter_names(self):
        """Baseline coordinate names followed by the added ones."""
        return joint_fit.PARAMETER_NAMES+self.added_coordinate_names

    def parameter_bounds(self):
        """Baseline optimizer bounds extended by the added coordinates."""
        lower, upper = joint_fit.parameter_bounds(self)
        added_lower, added_upper = [], []
        for name in self.added_coordinate_names:
            if name == ADDED_COORDINATE_NAMES[2]:
                added_lower.append(MASS_DIFFERENCE_BOUNDS_GEV[0])
                added_upper.append(MASS_DIFFERENCE_BOUNDS_GEV[1])
            else:
                width = self.odd_data['total_widths'][name.removeprefix('total_width_pull_')]
                added_lower.append(-.99*width['value_gev']/width['sigma_gev'])
                added_upper.append(np.inf)
        return np.r_[lower, added_lower], np.r_[upper, added_upper]

    def added_start(self):
        """Starting values of the added coordinates (PDG 2024 likelihood rows)."""
        return np.array([self.pdg_kstar_mass_difference_gev if name == ADDED_COORDINATE_NAMES[2]
                         else 0. for name in self.added_coordinate_names])

    def extend(self, baseline_values):
        """Append the added starting values to a 32-coordinate point."""
        return np.r_[np.asarray(baseline_values, dtype=float)[:BASELINE_COORDINATE_COUNT],
                     self.added_start()]

    def initial(self):
        """Baseline initial point extended by the added starting values."""
        return self.extend(super().initial())

    def kstar_mass_difference(self, values):
        """Delta m = m(K*0) - m(K*+) in GeV used with separate charge states."""
        if self.kstar_options['free_kstar_mass_difference']:
            return values[self.added_index[ADDED_COORDINATE_NAMES[2]]]
        return self.pdg_kstar_mass_difference_gev

    def unpack(self, values):
        """Inherited unpack, then the selected K* width treatment."""
        point, sector, mk, gv, totals = super().unpack(values)
        if self.kstar_default:
            return point, sector, mk, gv, totals
        if self.free_kstar_widths:
            for name in KSTAR_NAMES:
                width = self.odd_data['total_widths'][name]
                totals[name] = (width['value_gev']
                                + width['sigma_gev']*values[self.added_index['total_width_pull_'+name]])
        else:
            neutral_mass = mk+self.kstar_mass_difference(values)
            row = self.gv_data['rows'][1]
            physical = self.alda['physical_masses_gev']
            totals['Kstar0'] = gv**2*sum(gv_fit.width_coefficient(
                row['physical_vector_mass_GeV'], neutral_mass, physical[first], physical[second],
                self.kaon_flavors[flavor]) for first, second, flavor in KSTAR_MODES['Kstar0'])
        if min(totals.values()) <= 0:
            raise ValueError('nonpositive total width')
        return point, sector, mk, gv, totals

    def _babar_blocks_full(self, values, state, derivatives=False, channels=CHANNELS):
        """Charge-specific K* propagators when separate charge states are selected.

        The BaBar propagator cache of the baseline is keyed on the first 12
        and the four pole-residue coordinates only, so it is cleared whenever
        an added coordinate changes.
        """
        if self.kstar_default:
            return super()._babar_blocks_full(values, state, derivatives, channels)
        added = tuple(values[BASELINE_COORDINATE_COUNT:])
        if added != self.last_added_coordinates:
            self.babar_cache_key = None
            self.last_added_coordinates = added
        if not self.kstar_options['separate_kstar_charge_states']:
            return super()._babar_blocks_full(values, state, derivatives, channels)
        totals = state[4]
        mass_difference = self.kstar_mass_difference(values)

        def charge_specific(name, invariant, beta, mass, shared_width, reference):
            if name != 'kstar':
                return RUNNING_PROPAGATOR(name, invariant, beta, mass, shared_width, reference)
            # The only neutral K* invariant is the integrated t of K_S K pi.
            if invariant is self.phase['kskpi']['t']:
                neutral_mass = mass+mass_difference
                return RUNNING_PROPAGATOR(name, invariant, beta, neutral_mass, totals['Kstar0'],
                                          dict(reference, m_Kstar_GeV=neutral_mass))
            return RUNNING_PROPAGATOR(name, invariant, beta, mass, totals['Kstarplus'], reference)

        joint_fit.running_propagator = charge_specific
        try:
            return super()._babar_blocks_full(values, state, derivatives, channels)
        finally:
            joint_fit.running_propagator = RUNNING_PROPAGATOR

    def jacobian(self, values, finite_difference_steps=None):
        """Inherited Jacobian plus forward differences for added coordinates."""
        if not self.added_coordinate_names:
            return super().jacobian(values, finite_difference_steps)
        baseline_jacobian = super().jacobian(values, finite_difference_steps)
        reference = self.residual(values)
        added_columns = []
        for index in range(BASELINE_COORDINATE_COUNT, len(values)):
            step = (np.sqrt(np.finfo(float).eps)*max(1., abs(values[index]))
                    if finite_difference_steps is None else finite_difference_steps[index])
            shifted = values.copy()
            shifted[index] += step
            added_columns.append((self.residual(shifted)-reference)/step)
        return np.column_stack([baseline_jacobian, *added_columns])


def main():
    settings = (dict(), dict(separate_kstar_charge_states=True),
                dict(separate_kstar_charge_states=True, free_kstar_mass_difference=True),
                dict(kstar_width_treatment='free_total_widths'),
                dict(separate_kstar_charge_states=True, free_kstar_mass_difference=True,
                     kstar_width_treatment='free_total_widths'))
    report = []
    for options in settings:
        model = KstarChargeLikelihood(1, 3, 4, bin_block_size=4, **options)
        lower, upper = model.parameter_bounds()
        report.append(dict(options=model.kstar_options, coordinate_count=len(model.parameter_names),
                           added_coordinates=[dict(name=name, start=float(start), lower=float(low),
                                                   upper=float(high))
                                              for name, start, low, high in zip(
                                                  model.added_coordinate_names, model.added_start(),
                                                  lower[BASELINE_COORDINATE_COUNT:],
                                                  upper[BASELINE_COORDINATE_COUNT:])]))
    print(json.dumps(report, indent=1))


if __name__ == '__main__':
    main()
