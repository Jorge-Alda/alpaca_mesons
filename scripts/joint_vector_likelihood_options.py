"""Switchable variations of the joint vector-sector likelihood.

Run: .venv/bin/python scripts/joint_vector_likelihood_options.py
(prints the declared variation list and the row and coordinate counts of each).

The baseline likelihood in fit_babar_kkpi_global_vector_sector.py is left
byte-identical, because the saved fit ledgers verify its SHA-256. The class
below inherits from it; with all options at their defaults it calls the
inherited methods unchanged, so residuals and Jacobians are bitwise equal.

Units: GeV for masses, widths and energies; GeV^2 for invariants; nb for
cross sections; residuals are dimensionless. The 32 optimizer coordinates and
their meaning are those of GlobalLikelihood. Conventions as in Alda:2026nrd.
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import calibrate_odd_low_energy as calibration
import fit_babar_kkpi_global_vector_sector as joint_fit
import fit_vector_gv as gv_fit
from fit_babar_kkpi_g_coefficients import CHANNELS, F, MK, MPI, triangle

ROOT = joint_fit.ROOT
BLOCKS = (('odd', 0, 12), ('masses_and_leptonic', 12, 19),
          ('shared_widths', 19, 26), ('phi_KK', 26, 28),
          ('BaBar', 28, 73), ('normalization_constraints', 73, 75))
BASELINE_REGULATOR_GEV = 2.5
BASELINE_UPPER_ENERGY_GEV = 1.78
SOURCE_HYPOTHESES = ('single_trace_negative', 'single_trace_positive', 'universal')
TENSOR_COORDINATES = (24, 25)
SOURCE_SLOPE_COORDINATE = 5
# Stage-3 scenario matching the preliminary point (PDG charged rho and Kstar
# masses, phi electronic width OUR EVALUATION); used only to build starts.
PDG_POLES = 'data/vector_sector/pdg2024_excited_poles.json'
# Hard-coded baseline excited poles of fit_babar_kkpi_g_coefficients.propagator:
# PDG state -> (propagator species, excitation slot 0 or 1, mass GeV, width GeV).
# These are the values used by the baseline, which differ slightly from the
# PDG 2024 central values for K*(1680) (width) and phi(2170) (mass, width).
EXCITED_POLE_BASELINES = {
    'rho(1450)': ('rho', 0, 1.465, .400), 'rho(1700)': ('rho', 1, 1.720, .250),
    'K*(1410)': ('kstar', 0, 1.414, .232), 'K*(1680)': ('kstar', 1, 1.718, .322),
    'omega(1420)': ('omega', 0, 1.410, .290), 'omega(1650)': ('omega', 1, 1.670, .315),
    'phi(1680)': ('phi', 0, 1.680, .150), 'phi(2170)': ('phi', 1, 2.162, .100)}
# a2(1320) pole in the tensor term of the K_S K pi amplitude
# (fit_babar_kkpi_global_vector_sector.py, GlobalLikelihood._babar_blocks_full).
TENSOR_POLE_BASELINE = ('a2(1320)', 1.3182, .107)
TENSOR_CHANNEL = 'kskpi'
STAGE3_SCENARIO = {'rho_input': 'rho_charged', 'kstar_input': 'kstar_charged',
                   'phi_input': 'phi_evaluation'}


REFERENCE_PROPAGATOR = joint_fit.reference_propagator


def excited_pole_propagator(species, invariant, mass, width):
    """Single excited-pole propagator 1/(M^2 - s - i M Gamma(s)), s in GeV^2.

    Same running widths as fit_babar_kkpi_g_coefficients.propagator: P-wave
    pi pi phase space times sqrt(s)/M for rho states, P-wave K pi phase space
    times M^4/s^2 for K* states, constant width for omega and phi states.
    """
    invariant = np.asarray(invariant)
    if species == 'rho':
        pion_phase = np.maximum(1-4*MPI**2/np.maximum(invariant, 1e-20), 0)**1.5
        running = width*np.sqrt(np.maximum(invariant, 0))/mass*pion_phase/(1-4*MPI**2/mass**2)**1.5
    elif species == 'kstar':
        ratio = np.maximum(triangle(invariant, MK**2, MPI**2), 0)/triangle(mass**2, MK**2, MPI**2)
        running = np.where(invariant > (MK+MPI)**2,
                           width*mass**4/np.maximum(invariant, 1e-20)**2*ratio**1.5, 0)
    else:
        running = width
    return 1/(mass**2-invariant-1j*mass*running)


def pole_shift_variations():
    """The 36 +-1 sigma mass and width shifts of the fixed poles.

    Uncertainties are read from PDG_POLES (PDG 2024); asymmetric errors are
    applied as quoted (+ uses the upper, - the lower uncertainty).
    """
    states = json.loads((ROOT/PDG_POLES).read_text())['states']
    baselines = {state: values[2:] for state, values in EXCITED_POLE_BASELINES.items()}
    baselines[TENSOR_POLE_BASELINE[0]] = TENSOR_POLE_BASELINE[1:]
    variations = []
    for state, (mass, width) in baselines.items():
        slug = state.replace('*', 'star').replace('(', '_').replace(')', '')
        for quantity, baseline in (('mass', mass), ('width', width)):
            for sign, label, key in ((1, 'plus', 'uncertainty_plus_GeV'), (-1, 'minus', 'uncertainty_minus_GeV')):
                sigma = states[state][quantity][key]
                variations.append(dict(
                    name=f'pole_{slug}_{quantity}_{label}_one_sigma',
                    options=dict(pole_shift=dict(state=state, quantity=quantity, shift_gev=sign*sigma)),
                    baseline_value_gev=baseline, pdg_uncertainty_gev=sigma,
                    definition=(f'{state} {quantity} {baseline} GeV shifted by {sign*sigma:+} GeV '
                                f'(PDG 2024 {states[state][quantity]["pdg_value_type"]} uncertainty).'),
                    justification='Fixed pole input of the BaBar amplitude; tests the dependence on it.'))
    return variations


class VariedLikelihood(joint_fit.GlobalLikelihood):
    """GlobalLikelihood with explicit, individually switchable variations.

    Options (defaults reproduce the baseline):
      babar_upper_energy_gev: drop BaBar bins whose upper edge exceeds it.
      regulator_gev: Lambda in the amplitude regulator exp(-q^6/Lambda^6),
          entering the bin weights squared as exp(-2 q^6/Lambda^6), q = sqrt(s).
      source_hypothesis: domain and source matrix of Alda:2026nrd stage 3;
          'universal' sets z_V = 0 identically.
      inverse_mass_branch: 'near_ideal' or 'second'; the second branch
          recentres (Delta_V, kappa_V) on the stage-3 second root.
      phi_kk_correlation: False sets the -0.72 PDG correlation to zero.
      phi_electronic_width: 'evaluation' (PDG OUR EVALUATION) or 'average'
          (PDG OUR AVERAGE, pdg2024.json alternatives.phi_ee_average).
      babar_fractional_systematic: f added in quadrature, sigma^2 + (f*obs)^2,
          uncorrelated between bins (assumption).
      excluded_blocks: names in BLOCKS removed from the residual.
      pole_shift: None, or dict(state, quantity, shift_gev) moving the mass or
          the width of one fixed excited vector pole (EXCITED_POLE_BASELINES)
          or of the a2(1320) tensor pole by shift_gev from its hard-coded
          baseline value. The running-width shapes of the propagator are kept.
    """

    def __init__(self, energy_order=2, pair_order=16, angle_order=16, bin_block_size=None,
                 babar_upper_energy_gev=BASELINE_UPPER_ENERGY_GEV,
                 regulator_gev=BASELINE_REGULATOR_GEV,
                 source_hypothesis='single_trace_negative', inverse_mass_branch='near_ideal',
                 phi_kk_correlation=True, phi_electronic_width='evaluation',
                 babar_fractional_systematic=0., excluded_blocks=(), pole_shift=None):
        super().__init__(energy_order, pair_order, angle_order, bin_block_size)
        if source_hypothesis not in SOURCE_HYPOTHESES:
            raise ValueError(f'source hypothesis {source_hypothesis!r} is not implemented')
        if inverse_mass_branch not in ('near_ideal', 'second'):
            raise ValueError('inverse_mass_branch must be near_ideal or second')
        if phi_electronic_width not in ('evaluation', 'average'):
            raise ValueError('phi_electronic_width must be evaluation or average')
        unknown = set(excluded_blocks) - {name for name, _, _ in BLOCKS}
        if unknown:
            raise ValueError(f'unknown blocks {sorted(unknown)}')
        self.options = dict(babar_upper_energy_gev=babar_upper_energy_gev,
                            regulator_gev=regulator_gev, source_hypothesis=source_hypothesis,
                            inverse_mass_branch=inverse_mass_branch,
                            phi_kk_correlation=phi_kk_correlation,
                            phi_electronic_width=phi_electronic_width,
                            babar_fractional_systematic=babar_fractional_systematic,
                            excluded_blocks=tuple(excluded_blocks), pole_shift=pole_shift)
        self.excited_pole = self.tensor_pole = None
        if pole_shift is not None:
            state, quantity = pole_shift['state'], pole_shift['quantity']
            if quantity not in ('mass', 'width'):
                raise ValueError('pole_shift quantity must be mass or width')
            if state == TENSOR_POLE_BASELINE[0]:
                mass, width = TENSOR_POLE_BASELINE[1:]
                mass, width = ((mass+pole_shift['shift_gev'], width) if quantity == 'mass'
                               else (mass, width+pole_shift['shift_gev']))
                self.tensor_pole = (mass, width)
            elif state in EXCITED_POLE_BASELINES:
                species, slot, mass, width = EXCITED_POLE_BASELINES[state]
                shifted = ((mass+pole_shift['shift_gev'], width) if quantity == 'mass'
                           else (mass, width+pole_shift['shift_gev']))
                self.excited_pole = (species, slot, (mass, width), shifted)
            else:
                raise ValueError(f'unknown pole {state!r}')
        self.custom_unpack = (source_hypothesis != 'single_trace_negative'
                              or inverse_mass_branch != 'near_ideal')
        self.regulator_factor_exponent = (None if regulator_gev == BASELINE_REGULATOR_GEV else
                                          1/regulator_gev**6 - 1/BASELINE_REGULATOR_GEV**6)
        if inverse_mass_branch == 'second':
            branch = self.stage3_branches()[1]
            self.point = self.point.copy()
            self.point[2:4] = branch['delta_gev2'], branch['kappa']
            midpoint = 4*(self.point[1]**2-self.point[0]**2)/(
                3*(self.assumptions['xl_gev2']-self.assumptions['xs_gev2']))
            self.branch_sign = np.sign(self.point[3]-midpoint)
        if not phi_kk_correlation:
            covariance = np.array(self.gv_data['covariance'])[2:, 2:]
            self.phi_cholesky = np.linalg.cholesky(np.diag(np.diag(covariance)))
        if phi_electronic_width == 'average':
            alternative = self.pdg['alternatives']['phi_ee_average']
            self.even_target = self.even_target.copy()
            self.even_sigma = self.even_sigma.copy()
            self.even_target[6] = alternative['value_kev']*1e-6
            self.even_sigma[6] = alternative['sigma_kev']*1e-6
        if babar_fractional_systematic:
            for channel in CHANNELS:
                data = self.babar[channel]
                data['errors'] = np.sqrt(data['errors']**2
                                         + (babar_fractional_systematic*data['observed'])**2)
        mask = np.ones(75, dtype=bool)
        offset = 28
        for channel in CHANNELS:
            upper_edges = self.babar[channel]['edges'][:, 1]
            mask[offset:offset+len(upper_edges)] = upper_edges <= babar_upper_energy_gev + 1e-9
            offset += len(upper_edges)
        for name, start, stop in BLOCKS:
            if name in excluded_blocks:
                mask[start:stop] = False
        self.row_mask = mask
        self.all_rows = bool(mask.all())
        self.return_full_rows = False

    def stage3_branches(self):
        """Return the stage-3 inverse-mass roots, near-ideal one first."""
        scenario = next(item for item in self.stage3['scenarios']
                        if item['rho_input'] == STAGE3_SCENARIO['rho_input']
                        and item['kstar_input'] == STAGE3_SCENARIO['kstar_input'])
        return sorted(scenario['extended_mass_branches'], key=lambda item: abs(item['theta_deg']))

    def stage3_source(self, hypothesis, branch):
        """Stage-3 conditional F_V and z_V for a hypothesis on a mass root."""
        scenario = next(item for item in self.stage3['scenarios']
                        if item['rho_input'] == STAGE3_SCENARIO['rho_input']
                        and item['kstar_input'] == STAGE3_SCENARIO['kstar_input'])
        widths = next(item for item in scenario['width_tests']
                      if item['phi_input'] == STAGE3_SCENARIO['phi_input'])
        branches = scenario['extended_mass_branches']
        index = branches.index(self.stage3_branches()[0 if branch == 'near_ideal' else 1])
        case = next(item for item in widths['fixed_angle_cases']
                    if item['case'] == f'extended_mass_branch_{index}')
        key = {'universal': 'universal', 'single_trace_positive': 'single_trace_spurion',
               'single_trace_negative': 'single_trace_spurion_sign_alternative'}[hypothesis]
        return case[key]['FV_gev'], case[key].get('z_gev_inverse', 0.)

    def fixed_coordinates(self, tensor_products_zero=False):
        """Coordinates held fixed by a variation, mapped to their fixed values."""
        fixed = {}
        if self.options['source_hypothesis'] == 'universal':
            # z_V is set to zero inside unpack; the coordinate has no effect.
            fixed[SOURCE_SLOPE_COORDINATE] = float(
                -self.point[SOURCE_SLOPE_COORDINATE]/joint_fit.EVEN_SCALES[SOURCE_SLOPE_COORDINATE])
        if tensor_products_zero:
            fixed.update({index: 0. for index in TENSOR_COORDINATES})
        return fixed

    def unpack(self, values):
        """Baseline unpack, or the same steps with the selected source domain."""
        if not self.custom_unpack:
            return super().unpack(values)
        assumptions = self.assumptions
        point = self.point+joint_fit.EVEN_SCALES*values[:6]
        hypothesis = self.options['source_hypothesis']
        if hypothesis == 'universal':
            point[5] = 0.
        sector, mk = calibration.J.action_sector(point, hypothesis, assumptions)
        midpoint = 4*(point[1]**2-point[0]**2)/(3*(assumptions['xl_gev2']-assumptions['xs_gev2']))
        if np.sign(point[3]-midpoint) != self.branch_sign:
            raise ValueError('inverse-mass branch crossed')
        light = point[4]+2*point[5]*assumptions['xl_gev2']
        strange = point[4]+2*point[5]*assumptions['xs_gev2']
        if light <= 0:
            raise ValueError('selected source branch crossed')
        if hypothesis == 'single_trace_negative' and strange >= 0:
            raise ValueError('selected source branch crossed')
        if hypothesis == 'single_trace_positive' and strange <= 0:
            raise ValueError('selected source branch crossed')
        gv = values[6]
        if gv <= 0 or np.any(abs(values[26:30]) > 5):
            raise ValueError('selected GV/pole-residue domain crossed')
        totals = {name: self.odd_data['total_widths'][name]['value_gev']
                  + self.odd_data['total_widths'][name]['sigma_gev']*values[7+index]
                  for index, name in enumerate(joint_fit.FREE_TOTALS)}
        physical = self.alda['physical_masses_gev']
        for row, name in zip(self.gv_data['rows'][:2], ('Kstarplus', 'Kstar0')):
            modes = ([('K_charged', 'pi0', 'kstar_charged_K+_pi0'),
                      ('K_neutral', 'pi_charged', 'kstar_charged_K0_pi+')]
                     if name == 'Kstarplus' else
                     [('K_neutral', 'pi0', 'kstar_neutral_K0_pi0'),
                      ('K_charged', 'pi_charged', 'kstar_neutral_K+_pi-')])
            totals[name] = gv**2*sum(gv_fit.width_coefficient(
                row['physical_vector_mass_GeV'], mk, physical[first], physical[second],
                self.kaon_flavors[flavor]) for first, second, flavor in modes)
        if min(totals.values()) <= 0:
            raise ValueError('nonpositive total width')
        return point, sector, mk, gv, totals

    def shifted_reference_propagator(self, name, invariant, parameters, beta=0.):
        """Reference propagator with one excited pole moved.

        Adds beta_slot*(P_shifted - P_baseline) for the shifted pole, where P
        is the single excited-pole propagator of excited_pole_propagator.
        """
        propagator = REFERENCE_PROPAGATOR(name, invariant, parameters, beta)
        species, slot, baseline, shifted = self.excited_pole
        if name != species:
            return propagator
        coefficient = (beta if np.ndim(beta) else (beta, 0.))[slot]
        return propagator+coefficient*(excited_pole_propagator(species, invariant, *shifted)
                                       - excited_pole_propagator(species, invariant, *baseline))

    def shift_tensor_pole(self, values, result, channels):
        """Replace the a2(1320) pole in the K_S K pi amplitude and its prediction."""
        residual, predictions = result
        if TENSOR_CHANNEL not in channels:
            return result
        phase = self.phase[TENSOR_CHANNEL]
        q, s, t = phase['q'], phase['s'], phase['t']
        u = q+2*MK**2+MPI**2-s-t
        source_rho = self.babar_cached[TENSOR_CHANNEL][4]
        _, baseline_mass, baseline_width = TENSOR_POLE_BASELINE
        mass, width = self.tensor_pole
        amplitude = self.current_amplitudes[TENSOR_CHANNEL]+(t-u)*(
            1/(mass**2-s-1j*mass*width)-1/(baseline_mass**2-s-1j*baseline_mass*baseline_width))*(
            -2*values[24]/F**3+np.sqrt(2)*values[25]*source_rho/F**3)
        self.current_amplitudes[TENSOR_CHANNEL] = amplitude
        predictions[TENSOR_CHANNEL] = np.sum(phase['weight']*abs(amplitude)**2, axis=(1, 2, 3))
        residual = np.concatenate([
            (predictions[channel]*(1+joint_fit.RELATIVE_SYSTEMATICS[channel]*values[30+CHANNELS.index(channel)])
             - self.babar[channel]['observed'])/self.babar[channel]['errors'] for channel in channels])
        return residual, predictions

    def _babar_blocks_full(self, values, state, derivatives=False, channels=CHANNELS):
        """Apply the regulator rescaling and the pole shifts when selected."""
        if self.regulator_factor_exponent is None and self.options['pole_shift'] is None:
            return super()._babar_blocks_full(values, state, derivatives, channels)
        original = self.phase
        if self.regulator_factor_exponent is not None:
            self.phase = {channel: dict(phase, weight=phase['weight']*np.exp(
                -2*phase['q']**3*self.regulator_factor_exponent)) for channel, phase in original.items()}
        if self.excited_pole is not None:
            joint_fit.reference_propagator = self.shifted_reference_propagator
        try:
            result = super()._babar_blocks_full(values, state, derivatives, channels)
            if derivatives or self.tensor_pole is None:
                return result
            return self.shift_tensor_pole(values, result, channels)
        finally:
            self.phase = original
            joint_fit.reference_propagator = REFERENCE_PROPAGATOR

    def full_residual(self, values):
        """All 75 rows of the varied likelihood, including removed rows."""
        return super().residual(values)

    def residual(self, values):
        """Rows retained by the energy window and block exclusions."""
        full = super().residual(values)
        if self.all_rows or self.return_full_rows:
            return full
        return full[self.row_mask]

    def jacobian(self, values, finite_difference_steps=None):
        """Inherited derivatives evaluated on all rows, then restricted."""
        self.return_full_rows = True
        try:
            jacobian = super().jacobian(values, finite_difference_steps)
        finally:
            self.return_full_rows = False
        return jacobian if self.all_rows else jacobian[self.row_mask]


def main():
    model = VariedLikelihood(1, 3, 4, bin_block_size=4)
    print(json.dumps(dict(row_count=int(model.row_mask.sum()),
                          options=model.options, blocks=BLOCKS), indent=1))


if __name__ == '__main__':
    main()
