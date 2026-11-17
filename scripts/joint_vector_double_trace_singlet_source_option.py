"""Switchable double-trace singlet-source structure y_V in the joint vector likelihood.

Run: .venv/bin/python scripts/joint_vector_double_trace_singlet_source_option.py
(prints the added coordinate, its start and bounds, and the photon residues and
BaBar source factors at the best accepted endpoint for y_V = 0 and 0.022 GeV^-1).

No fit is run by this module. The baseline likelihood in
fit_babar_kkpi_global_vector_sector.py and the variations in
joint_vector_likelihood_options.py are left unchanged (the saved ledgers check
their SHA-256); the class below inherits from VariedLikelihood.

Operator (memory/rcht/odd-parity/joint-vector-double-trace-source-mass-and-kkpi.md,
scripts/derive_double_trace_vector_mass_and_kkpi.py):
    L_y = y_V/(2 sqrt2) <V_{mu nu}> <chi_+ f_+^{mu nu}>,
which adds y_V <chi_+ Q> 1 to the photon source matrix,
    S_gamma = F_V Q + z_V {chi_+, Q} + y_V <chi_+ Q> 1.
At order e^0 it gives no vector mass term. With chi_+ = diag(x_l, x_l, x_s):
  * photon residues: R_rho = F_l/sqrt2, R_omega_I = (F_l + 2 y_V (x_l - x_s))/(3 sqrt2),
    R_phi_I = -(F_s - y_V (x_l - x_s))/3, with F_l = F_V + 2 z_V x_l,
    F_s = F_V + 2 z_V x_s. They enter the omega, phi leptonic widths and the
    photon-insertion parts of the odd radiative rows through
    vector_sector.source_residues (the 'double_trace' hypothesis of
    derive_odd_low_energy_jacobians.action_sector, with point[6] = y_V);
  * BaBar KK pi current: F_l -> F_l + 2 y_V (x_l - x_s) in the omega_I insertion
    and F_s -> F_s - y_V (x_l - x_s) in the phi_I insertion; the rho^0
    insertion keeps F_l. The baseline uses one factor `light` for both the
    rho^0 and the omega_I insertions; with the option on, the BaBar amplitude
    is evaluated by a copy of the baseline method in which the isovector
    factor (source_rho, also used by the a2 tensor term) and the isoscalar
    factors (source_light, source_strange) are separate.

Option:
  double_trace_singlet_source (default False): appends the coordinate
      ADDED_COORDINATE_NAMES = (y_V in GeV^-1,) after the 32 baseline
      coordinates. Start 0; bounds [-COEFFICIENT_BOUND_GEV_INVERSE,
      +COEFFICIENT_BOUND_GEV_INVERSE] (reason in
      memory/rcht/odd-parity/joint-vector-double-trace-singlet-source-fit-plan.md);
      outside the bounds unpack raises and the residual becomes the baseline
      penalty vector. Nominal degrees of freedom: 75 residuals minus 33
      coordinates = 42 (baseline: 43).
  With the option off, every method is the inherited one, so residuals and
  Jacobians are bitwise equal to the baseline. With the option on and y_V = 0
  the residual and the first 32 Jacobian columns are also bitwise equal.

Assumptions (stated, inherited from the 4.1 record):
  * O(e^2) electromagnetic mass shifts are omitted for every source term.
  * Exact strong isospin x_u = x_d = x_l, as in the baseline.
  * The excited vector poles of the BaBar current keep their fixed couplings.
  * The baseline source-branch check (F_l > 0, F_s < 0) is kept unchanged; it
    acts on F_V and z_V only, not on the y_V-shifted isoscalar factors.
  * Only the default VariedLikelihood options that act inside unpack and the
    BaBar amplitude (source hypothesis, mass branch, regulator, pole shifts)
    are compatible with the option; other combinations are refused.

Units: GeV for masses, widths, energies and F_V; GeV^2 for invariants and x_q;
GeV^-1 for z_V and y_V; residuals are dimensionless. Conventions of Alda:2026nrd.
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import calibrate_odd_low_energy as calibration  # noqa: E402
import fit_babar_kkpi_global_vector_sector as joint_fit  # noqa: E402
import joint_vector_likelihood_options as options_module  # noqa: E402
from fit_babar_kkpi_g_coefficients import CHANNELS, F, MK, MPI, MV  # noqa: E402

BASELINE_COORDINATE_COUNT = 32
COORDINATE_INDEX = 32
# Captured at import: drivers temporarily replace joint_fit.PARAMETER_NAMES to label
# 33 Jacobian columns in the baseline quadrature audit.
BASELINE_PARAMETER_NAMES = joint_fit.PARAMETER_NAMES
ADDED_COORDINATE_NAMES = ('double_trace_singlet_source_y_V_GeV_inverse',)
COEFFICIENT_BOUND_GEV_INVERSE = .5
# Magnitude of the stage-3 near-ideal double-trace value y_V (data/vector_sector/stage3_results.json,
# rho_charged/kstar_charged/phi_evaluation, extended_mass_branch_1), rounded; used for offset starts.
STAGE3_SCALE_GEV_INVERSE = .022


def isoscalar_source_factors(light, strange, coefficient, light_spurion, strange_spurion):
    """Isoscalar photon-insertion factors of the BaBar current with y_V.

    light = F_l and strange = F_s in GeV, coefficient = y_V in GeV^-1,
    light_spurion = x_l and strange_spurion = x_s in GeV^2. Returns
    (F_l + 2 y_V (x_l - x_s), F_s - y_V (x_l - x_s)), the factors of the
    omega_I and phi_I insertions; the rho^0 insertion keeps F_l.
    """
    difference = light_spurion-strange_spurion
    return light+2*coefficient*difference, strange-coefficient*difference


class DoubleTraceSingletSourceLikelihood(options_module.VariedLikelihood):
    """VariedLikelihood with the switchable double-trace singlet-source coordinate y_V.

    See the module docstring for the residue shifts, the BaBar factors, the
    coordinate, its bounds and the assumptions.
    """

    def __init__(self, energy_order=2, pair_order=16, angle_order=16, bin_block_size=None,
                 double_trace_singlet_source=False, **variation_options):
        super().__init__(energy_order, pair_order, angle_order, bin_block_size, **variation_options)
        self.double_trace_singlet_source = bool(double_trace_singlet_source)
        if self.double_trace_singlet_source and (
                self.custom_unpack or self.regulator_factor_exponent is not None
                or self.options['pole_shift'] is not None):
            raise ValueError('double_trace_singlet_source is implemented for the baseline source hypothesis, '
                             'mass branch, regulator and fixed poles only')
        self.added_coordinate_names = ADDED_COORDINATE_NAMES if self.double_trace_singlet_source else ()

    @property
    def parameter_names(self):
        """Baseline coordinate names followed by the added one."""
        return BASELINE_PARAMETER_NAMES+self.added_coordinate_names

    def parameter_bounds(self):
        """Baseline optimizer bounds extended by the y_V bounds in GeV^-1."""
        lower, upper = joint_fit.parameter_bounds(self)
        count = len(self.added_coordinate_names)
        return (np.r_[lower, np.full(count, -COEFFICIENT_BOUND_GEV_INVERSE)],
                np.r_[upper, np.full(count, COEFFICIENT_BOUND_GEV_INVERSE)])

    def extend(self, baseline_values, coefficient=0.):
        """Append y_V in GeV^-1 to a 32-coordinate point (option on only)."""
        if not self.double_trace_singlet_source:
            raise ValueError('the coordinate y_V exists only with double_trace_singlet_source=True')
        return np.r_[np.asarray(baseline_values, dtype=float)[:BASELINE_COORDINATE_COUNT], float(coefficient)]

    def initial(self):
        """Baseline initial point, extended by y_V = 0 when the option is on."""
        values = super().initial()
        return self.extend(values) if self.double_trace_singlet_source else values

    def check_length(self, values):
        """Reject points whose length does not match the coordinate names."""
        if len(values) != len(self.parameter_names):
            raise ValueError(f'expected {len(self.parameter_names)} coordinates, got {len(values)}')

    def unpack(self, values):
        """Inherited unpack; with the option on, the sector is rebuilt with y_V in the residues."""
        state = super().unpack(values)
        if not self.double_trace_singlet_source:
            return state
        coefficient = values[COORDINATE_INDEX]
        if abs(coefficient) > COEFFICIENT_BOUND_GEV_INVERSE:
            raise ValueError('double-trace singlet-source domain crossed')
        point, _, mk, gv, totals = state
        sector, _ = calibration.J.action_sector(np.r_[point, coefficient], 'double_trace', self.assumptions)
        return point, sector, mk, gv, totals

    def _babar_blocks_full(self, values, state, derivatives=False, channels=CHANNELS):
        """BaBar current; with the option on, separate isovector and isoscalar source factors.

        Option off: the inherited method. Option on: a copy of
        GlobalLikelihood._babar_blocks_full in which source_rho uses F_l and
        source_light, source_strange use isoscalar_source_factors; the
        propagator cache key also contains y_V. At y_V = 0 the arithmetic is
        identical to the baseline (tested bitwise).
        """
        if not self.double_trace_singlet_source:
            return super()._babar_blocks_full(values, state, derivatives, channels)
        point, sector, mk, gv, totals = state
        a, b, p, c3, c4, d1, d2, d3, g, g2, g4, bg = (np.eye(12).reshape(12, 12, 1, 1, 1, 1) if derivatives
                                                      else values[12:24])
        cm = p+8*c3
        dm = d1+8*d2-d3
        light_spurion, strange_spurion = self.assumptions['xl_gev2'], self.assumptions['xs_gev2']
        light = point[4]+2*point[5]*light_spurion
        strange = point[4]+2*point[5]*strange_spurion
        light_isoscalar, strange_isoscalar = isoscalar_source_factors(
            light, strange, values[COORDINATE_INDEX], light_spurion, strange_spurion)
        reference = dict(self.preliminary, m_rho_GeV=point[0], m_Kstar_GeV=mk,
                         Gamma_total_omega_GeV=totals['omega'], Gamma_total_phi_GeV=totals['phi'])
        light_generator = np.diag([1, 1, 0])/np.sqrt(2)
        strange_generator = np.diag([0, 0, 1])
        rotation = np.array([[np.trace(generator@basis).real for generator in sector.generators[1:]]
                             for basis in (light_generator, strange_generator)])

        def excited(name, invariant, beta):
            ground, with_excitations = joint_fit.reference_propagator_pair(name, invariant, beta, reference)
            return with_excitations-ground

        def neutral(invariant, beta):
            omega = 1/(sector.masses[1]**2-invariant-1j*sector.masses[1]*totals['omega'])
            phi = 1/(sector.masses[2]**2-invariant-1j*sector.masses[2]*totals['phi'])
            ll = rotation[0, 0]**2*omega+rotation[0, 1]**2*phi+excited('omega', invariant, beta)
            ss = rotation[1, 0]**2*omega+rotation[1, 1]**2*phi+excited('phi', invariant, beta)
            ls = rotation[0, 0]*rotation[1, 0]*omega+rotation[0, 1]*rotation[1, 1]*phi
            return ll, ss, ls

        def running(name, invariant, beta):
            mass = point[0] if name == 'rho' else mk
            shared = totals['rho0'] if name == 'rho' else (totals['Kstar0']+totals['Kstarplus'])/2
            return joint_fit.running_propagator(name, invariant, beta, mass, shared, reference)
        cache_key = tuple(values[:12])+tuple(values[26:30])+(values[COORDINATE_INDEX],)
        if getattr(self, 'babar_cache_key', None) != cache_key:
            self.babar_cache_key = cache_key
            self.babar_cached = {}
        if not derivatives:
            self.current_amplitudes = {}
        predictions = {}
        for channel in channels:
            index = CHANNELS.index(channel)
            phase = self.phase[channel]
            q, s, t = phase['q'], phase['s'], phase['t']
            u = q+2*MK**2+MPI**2-s-t
            beta = (values[26+index], values[28+index])
            if channel not in self.babar_cached:
                rq = running('rho', q, beta)
                rs = running('rho', s, beta)
                llq, ssq, lsq = neutral(q, beta)
                lls, _, lss = neutral(s, beta)
                lls = lls-np.sqrt(2)*lss
                kt = running('kstar', t, beta)
                ku = running('kstar', u, beta)
                source_rho = light*rq
                source_light = light_isoscalar*llq-np.sqrt(2)*strange_isoscalar*lsq
                source_strange = strange_isoscalar*ssq-light_isoscalar*lsq/np.sqrt(2)
                self.babar_cached[channel] = (rs, lls, kt, ku, source_rho, source_light, source_strange)
            rs, lls, kt, ku, source_rho, source_light, source_strange = self.babar_cached[channel]

            def cr(pair, mass2):
                return a*q-b*pair+cm*mass2

            def crk(pair):
                return cr(pair, MK**2)+24*c4*(MK**2-MPI**2)

            def dr(pair, mass2):
                return dm*mass2+d3*(q+pair)
            gr1 = g*(q+s-4*MK**2-MPI**2)+4*bg*MK**2-4*g2*(q-2*MK**2-MPI**2)-4*g4*(MK**2-MPI**2)
            gr2 = g*(q-s-MPI**2)+2*bg*MPI**2-2*g2*(q-2*MK**2-MPI**2)+4*g4*(MK**2-MPI**2)
            contact = 2*np.sqrt(2)/(F**3*MV)
            direct = 2*np.sqrt(2)*gv/(3*F**3*MV)
            exchange = 2*gv/F**3
            if channel == 'kkpi0':
                amplitude = ((0. if derivatives else 3/(12*np.pi**2*F**3))
                             + contact*((source_rho+source_light/3)*gr1-2*source_strange*gr2/3))
                amplitude = amplitude-direct*((rs+3*lls)*cr(s, MPI**2)+crk(t)*kt+crk(u)*ku)
                amplitude = amplitude+exchange*(2*dr(s, MPI**2)*(source_rho*lls+source_light*rs/3)
                                                + (source_rho+source_light/3-2*source_strange/3)
                                                * (dr(t, MK**2)*kt+dr(u, MK**2)*ku))
            else:
                amplitude = contact*(g*(t-u)*source_rho-source_light*gr1/3+2*source_strange*gr2/3)
                amplitude = amplitude+direct*(rs*cr(s, MPI**2)-2*cr(t, MK**2)*kt+crk(u)*ku)
                amplitude = amplitude-exchange*(2*source_light*rs*dr(s, MPI**2)/3
                                                - (source_rho-source_light/3+2*source_strange/3)*dr(t, MK**2)*kt
                                                + (source_rho+source_light/3-2*source_strange/3)*dr(u, MK**2)*ku)
                if not derivatives:
                    amplitude = amplitude+(t-u)/(1.3182**2-s-1j*1.3182*.107)*(
                        -2*values[24]/F**3+np.sqrt(2)*values[25]*source_rho/F**3)
            if derivatives:
                predictions[channel] = np.sum(phase['weight'][None, ...]*2*np.real(
                    self.current_amplitudes[channel][None, ...].conj()*amplitude), axis=(2, 3, 4)).T
            else:
                self.current_amplitudes[channel] = amplitude
                predictions[channel] = np.sum(phase['weight']*abs(amplitude)**2, axis=(1, 2, 3))
        if derivatives:
            return np.vstack([predictions[channel]
                              * (1+joint_fit.RELATIVE_SYSTEMATICS[channel]*values[30+CHANNELS.index(channel)])
                              / self.babar[channel]['errors'][:, None] for channel in channels])
        residual = np.concatenate([
            (predictions[channel]*(1+joint_fit.RELATIVE_SYSTEMATICS[channel]*values[30+CHANNELS.index(channel)])
             - self.babar[channel]['observed'])/self.babar[channel]['errors'] for channel in channels])
        return residual, predictions

    def residual(self, values):
        """Inherited residual after a check of the coordinate count."""
        self.check_length(values)
        return super().residual(values)

    def jacobian(self, values, finite_difference_steps=None):
        """Inherited Jacobian plus a forward-difference column for y_V."""
        self.check_length(values)
        baseline_jacobian = super().jacobian(values, finite_difference_steps)
        if not self.double_trace_singlet_source:
            return baseline_jacobian
        reference = self.residual(values)
        step = (np.sqrt(np.finfo(float).eps)*max(1., abs(values[COORDINATE_INDEX]))
                if finite_difference_steps is None else finite_difference_steps[COORDINATE_INDEX])
        shifted = np.array(values, dtype=float)
        shifted[COORDINATE_INDEX] += step
        return np.column_stack([baseline_jacobian, (self.residual(shifted)-reference)/step])


def main():
    import refit_joint_vector_sensitivity as sensitivity
    model = DoubleTraceSingletSourceLikelihood(1, 3, 4, bin_block_size=4, double_trace_singlet_source=True)
    lower, upper = model.parameter_bounds()
    best = sensitivity.accepted_context()['accepted'][0]
    point = model.point+joint_fit.EVEN_SCALES*np.asarray(best['endpoint'])[:6]
    light_spurion, strange_spurion = model.assumptions['xl_gev2'], model.assumptions['xs_gev2']
    light = point[4]+2*point[5]*light_spurion
    strange = point[4]+2*point[5]*strange_spurion
    evaluations = {}
    for coefficient in (0., STAGE3_SCALE_GEV_INVERSE):
        _, sector, _, _, _ = model.unpack(model.extend(best['endpoint'], coefficient))
        evaluations[f'{coefficient:g}'] = dict(
            physical_residues_rho_omega_phi_gev=list(sector.residues),
            isoscalar_source_factors_gev=list(isoscalar_source_factors(
                light, strange, coefficient, light_spurion, strange_spurion)),
            isovector_source_factor_gev=light)
    print(json.dumps(dict(
        coordinate_count=len(model.parameter_names), nominal_degrees_of_freedom=75-len(model.parameter_names),
        added_coordinate=dict(name=ADDED_COORDINATE_NAMES[0], start=float(model.initial()[COORDINATE_INDEX]),
                              lower=float(lower[COORDINATE_INDEX]), upper=float(upper[COORDINATE_INDEX])),
        endpoint=best['id'], evaluations_by_y_V_gev_inverse=evaluations), indent=1))


if __name__ == '__main__':
    main()
