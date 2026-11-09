"""Switchable isospin-breaking VPP vertex in the joint vector likelihood.

Run: .venv/bin/python scripts/joint_vector_isospin_breaking_vpp_option.py
(prints the coordinate names, starting values and bounds, and the per-mode
vertex factors T(r)/T(0) at the scan point r = (-0.350, 1.855, 0.729) GeV^-2).

No fit is run by this module. The baseline likelihood in
fit_babar_kkpi_global_vector_sector.py and the variations in
joint_vector_likelihood_options.py are left unchanged; the class below
inherits from VariedLikelihood. With the option off (default) it calls the
inherited methods unchanged, so residuals and Jacobians are bitwise equal to
the baseline. With the option on and r = (0, 0, 0) the residual and the first
32 Jacobian columns are also bitwise equal to the baseline.

Model (memory/rcht/odd-parity/joint-vector-isospin-breaking-vpp-operators.md,
scripts/derive_isospin_breaking_vpp_operators.py): the three O(p^4) operators
    i lambda_1 <V_{mu nu} {u^mu u^nu, chi_+}> + i lambda_2 <V_{mu nu} u^mu chi_+ u^nu>
    + i lambda_3 <V_{mu nu}> <u^mu u^nu chi_+>
keep the LO VPP momentum structure; the flavor trace of each mode becomes
    T(r) = tr(V[P1,P2]) + r1 tr(V{[P1,P2],X}) + r2 tr(V(P1 X P2 - P2 X P1))
           + r3 tr(V) tr([P1,P2] X),
with X = chi_+ at zero fields = 4 diag(B0 m_u, B0 m_d, B0 m_s) (Alda:2026nrd
values in GeV^2), r_i = sqrt2 lambda_i / G_V in GeV^-2, the published pi0
diagonals of Alda:2026nrd, and the flavor matrices of fit_vector_gv.flavor_factors.
Each VPP coupling of a mode is multiplied by its factor T(r)/T(0); each
partial width by the square of that factor.

Option:
  isospin_breaking_vpp (default False): adds the three coordinates
      ADDED_COORDINATE_NAMES = r1, r2, r3 in GeV^-2, appended after the 32
      baseline coordinates. Start (0, 0, 0); bounds [-5, 5] GeV^-2 each (the
      range of scripts/scan_isospin_breaking_vpp_simultaneous.py). Nominal
      degrees of freedom: 73 measurements + 2 normalization constraints = 75
      residuals minus 35 coordinates = 40 (baseline: 43).
  vpp_kstar_excitation_scaling (default 'ground_only'; only with the option
      on): the K* propagators kt, ku of the BaBar amplitude contain the
      K*(892) pole plus beta-weighted K*(1410), K*(1680) terms. 'ground_only'
      multiplies only the K*(892) pole by the mode factor; 'whole_propagator'
      multiplies the whole propagator. This choice is open in
      joint-vector-kstar-kpi-saturation-freedom.md ("Open choice: excited K*
      terms"); the default is an assumption, not a result.

Where the factors enter (implemented):
  * K*+ and K*0 total-width rows: each K pi mode of
    Gamma = G_V^2 sum width_coefficient is multiplied by its squared factor
    (K*+: K+ pi0, K0 pi+; K*0: K0 pi0, K+ pi-).
  * K*0 -> K0 gamma and K*+ -> K+ gamma branching-fraction rows: through
    their denominators, the modified K* total widths. The radiative
    amplitudes themselves contain no VPP vertex and are unchanged.
  * BaBar K* propagator widths: through the baseline charge-mean width
    (Gamma(K*0) + Gamma(K*+))/2 built from the modified totals.
  * BaBar K*-channel terms (`direct` and `exchange` times kt or ku): G_V is
    replaced by G_V T(r)/T(0) of the K* -> K pi vertex of that term,
    implemented by multiplying the K* propagator, exactly as the sqrt(kappa_q)
    rule of joint-vector-kstar-kpi-saturation-freedom.md (Step 3), but with
    one factor per decay mode instead of one per charge state:
      K+ K- pi0, t and u (K+- pi0): 'kstar_charged_K+_pi0';
      K_S K pi, t (K+- pi-+, neutral K*): 'kstar_neutral_K+_pi-';
      K_S K pi, u (K_S pi+-, charged K*): 'kstar_charged_K0_pi+'.
    C-conjugate modes have the same factor (the operators are C-even).
  * phi -> K+K- and phi -> K_L K_S branching rows: the prediction is
    multiplied by the squared factor of phi -> K+K- and phi -> K0 K0bar. The
    phi flavor matrix is the fitted phi generator sector.generators[2], the
    same matrix whose LO trace the baseline rows use (choice stated here: the
    scan script uses the fixed provisional mixing angle instead). The phi
    total width is a free pull in the baseline and is not changed.

NOT implemented (left at the baseline):
  * BaBar s-channel terms with rs, lls (rho, omega, phi -> K Kbar vertices,
    including the charged rho -> K Kbar vertex of K_S K pi). The baseline
    writes these vertices in the light/strange basis of the rotated neutral
    propagator, and r1, r2, r3 act differently on its light and strange parts
    and, through r3, on the trace of V. The records give the factors only
    for the physical rho0, omega -> K Kbar modes at ideal mixing and not for
    the charged rho, so this substitution is not unambiguous and is not made.
  * Running-width shapes: unchanged (on-shell widths only are rescaled).
  * K* masses, rho/omega/eta' widths (free pulls), mass and leptonic rows.

Units: GeV for masses, widths and energies; GeV^2 for invariants and X;
GeV^-2 for r_i; residuals are dimensionless. Conventions as in Alda:2026nrd.
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fit_babar_kkpi_global_vector_sector as joint_fit  # noqa: E402
import fit_vector_gv as gv_fit  # noqa: E402
import joint_vector_likelihood_options as options_module  # noqa: E402
from fit_babar_kkpi_g_coefficients import CHANNELS  # noqa: E402

BASELINE_COORDINATE_COUNT = 32
ADDED_COORDINATE_NAMES = ('isospin_breaking_vpp_r1_GeV_minus2', 'isospin_breaking_vpp_r2_GeV_minus2',
                          'isospin_breaking_vpp_r3_GeV_minus2')
COEFFICIENT_BOUND_GEV_MINUS2 = 5.
EXCITATION_SCALINGS = ('ground_only', 'whole_propagator')
KSTAR_NAMES = ('Kstarplus', 'Kstar0')
# K pi modes of each K* charge state: (first daughter, second daughter, mode key).
KSTAR_MODES = {
    'Kstarplus': (('K_charged', 'pi0', 'kstar_charged_K+_pi0'),
                  ('K_neutral', 'pi_charged', 'kstar_charged_K0_pi+')),
    'Kstar0': (('K_neutral', 'pi0', 'kstar_neutral_K0_pi0'),
               ('K_charged', 'pi_charged', 'kstar_neutral_K+_pi-'))}
# phi -> K Kbar likelihood rows in the order of gv_measurements.json rows 2 and 3.
PHI_MODES = ('phi_K+_K-', 'phi_K0_K0bar')
SCAN_POINT_GEV_MINUS2 = (-.350, 1.855, .729)
RUNNING_PROPAGATOR = joint_fit.running_propagator


def matrix_unit(row, column):
    """3x3 flavor matrix with one unit entry (numpy version of fit_vector_gv.matrix_unit)."""
    result = np.zeros((3, 3))
    result[row, column] = 1.
    return result


def mode_matrices(pion_zero_diagonal, phi_generator=None):
    """(V, P1, P2) flavor matrices of the modes used here, as in flavor_factors.

    pion_zero_diagonal: the three lambda diagonals of the pi0 (the matrix is
    diag(...)/sqrt2). phi_generator: 3x3 phi flavor matrix; phi modes are
    omitted when it is None.
    """
    unit = matrix_unit
    pion_zero = np.diag(pion_zero_diagonal)/np.sqrt(2)
    matrices = {
        'kstar_charged_K+_pi0': (unit(2, 0), unit(0, 2), pion_zero),
        'kstar_charged_K0_pi+': (unit(2, 0), unit(1, 2), unit(0, 1)),
        'kstar_neutral_K0_pi0': (unit(2, 1), unit(1, 2), pion_zero),
        'kstar_neutral_K+_pi-': (unit(2, 1), unit(0, 2), unit(1, 0)),
        'rho0_pi+_pi-': (np.diag([1., -1., 0.])/np.sqrt(2), unit(0, 1), unit(1, 0)),
        'rho+_pi+_pi0': (unit(1, 0), unit(0, 1), pion_zero)}
    if phi_generator is not None:
        matrices['phi_K+_K-'] = (phi_generator, unit(0, 2), unit(2, 0))
        matrices['phi_K0_K0bar'] = (phi_generator, unit(1, 2), unit(2, 1))
    return matrices


def effective_trace(vector, first, second, spurion, coefficients):
    """T(r) of one mode; spurion is the 3x3 X in GeV^2, coefficients (r1, r2, r3) in GeV^-2."""
    commutator = first@second-second@first
    r1, r2, r3 = coefficients
    return (np.trace(vector@commutator)
            + r1*np.trace(vector@(commutator@spurion+spurion@commutator))
            + r2*np.trace(vector@(first@spurion@second-second@spurion@first))
            + r3*np.trace(vector)*np.trace(commutator@spurion))


def vertex_factors(coefficients, spurion_diagonal, pion_zero_diagonal, phi_generator=None):
    """T(r)/T(0) for every mode of mode_matrices (dimensionless, signed).

    coefficients: (r1, r2, r3) in GeV^-2; spurion_diagonal: X = 4 B0 m_q in GeV^2.
    """
    coefficients = np.asarray(coefficients, dtype=float)
    if coefficients.shape != (3,):
        raise ValueError('three coefficients r1, r2, r3 are required')
    spurion = np.diag(np.asarray(spurion_diagonal, dtype=float))
    zero = np.zeros(3)
    return {name: float(np.real(effective_trace(*matrices, spurion, coefficients)
                                / effective_trace(*matrices, spurion, zero)))
            for name, matrices in mode_matrices(pion_zero_diagonal, phi_generator).items()}


class IsospinBreakingVppLikelihood(options_module.VariedLikelihood):
    """VariedLikelihood with the switchable isospin-breaking VPP vertex.

    See the module docstring for the vertex factors, the rows and BaBar terms
    they enter, the terms left at the baseline, the coordinates, starting
    values, bounds and degrees of freedom. All VariedLikelihood options remain
    available.
    """

    def __init__(self, energy_order=2, pair_order=16, angle_order=16, bin_block_size=None,
                 isospin_breaking_vpp=False, vpp_kstar_excitation_scaling='ground_only',
                 **variation_options):
        super().__init__(energy_order, pair_order, angle_order, bin_block_size, **variation_options)
        if vpp_kstar_excitation_scaling not in EXCITATION_SCALINGS:
            raise ValueError(f'unknown K* excitation scaling {vpp_kstar_excitation_scaling!r}')
        if vpp_kstar_excitation_scaling != 'ground_only' and not isospin_breaking_vpp:
            raise ValueError('the K* excitation scaling requires isospin_breaking_vpp=True')
        self.isospin_breaking_vpp = bool(isospin_breaking_vpp)
        self.excitation_scaling = vpp_kstar_excitation_scaling
        self.added_coordinate_names = ADDED_COORDINATE_NAMES if self.isospin_breaking_vpp else ()
        self.spurion_diagonal = [4*value for value in self.alda['B0mq_gev2']]
        self.pion_zero_diagonal = self.alda['published_lambda_diagonals']['pi0']
        self.last_added_coordinates = None

    @property
    def parameter_names(self):
        """Baseline coordinate names followed by the added ones."""
        return joint_fit.PARAMETER_NAMES+self.added_coordinate_names

    def parameter_bounds(self):
        """Baseline optimizer bounds extended by [-5, 5] GeV^-2 for each r_i."""
        lower, upper = joint_fit.parameter_bounds(self)
        count = len(self.added_coordinate_names)
        return (np.r_[lower, np.full(count, -COEFFICIENT_BOUND_GEV_MINUS2)],
                np.r_[upper, np.full(count, COEFFICIENT_BOUND_GEV_MINUS2)])

    def extend(self, baseline_values, coefficients=(0., 0., 0.)):
        """Append (r1, r2, r3) in GeV^-2 to a 32-coordinate point (option on only)."""
        if not self.isospin_breaking_vpp:
            raise ValueError('coordinates r1, r2, r3 exist only with isospin_breaking_vpp=True')
        coefficients = np.asarray(coefficients, dtype=float)
        if coefficients.shape != (3,):
            raise ValueError('three coefficients r1, r2, r3 are required')
        return np.r_[np.asarray(baseline_values, dtype=float)[:BASELINE_COORDINATE_COUNT], coefficients]

    def initial(self):
        """Baseline initial point, extended by r = (0, 0, 0) when the option is on."""
        values = super().initial()
        return self.extend(values) if self.isospin_breaking_vpp else values

    def coefficients(self, values):
        """(r1, r2, r3) in GeV^-2, or None when the option is off or all r_i vanish."""
        if not self.isospin_breaking_vpp:
            return None
        coefficients = np.asarray(values[BASELINE_COORDINATE_COUNT:BASELINE_COORDINATE_COUNT+3], dtype=float)
        return coefficients if np.any(coefficients) else None

    def factors(self, coefficients, phi_generator=None):
        """Per-mode T(r)/T(0) with the Alda:2026nrd X and published pi0."""
        return vertex_factors(coefficients, self.spurion_diagonal, self.pion_zero_diagonal, phi_generator)

    def check_length(self, values):
        """Reject points whose length does not match the coordinate names."""
        if len(values) != len(self.parameter_names):
            raise ValueError(f'expected {len(self.parameter_names)} coordinates, got {len(values)}')

    def unpack(self, values):
        """Inherited unpack, then K* total widths with the squared mode factors."""
        state = super().unpack(values)
        if self.isospin_breaking_vpp and np.any(
                np.abs(values[BASELINE_COORDINATE_COUNT:]) > COEFFICIENT_BOUND_GEV_MINUS2):
            raise ValueError('isospin-breaking VPP coefficient domain crossed')
        coefficients = self.coefficients(values)
        if coefficients is None:
            return state
        point, sector, mk, gv, totals = state
        factors = self.factors(coefficients)
        physical = self.alda['physical_masses_gev']
        for row, name in zip(self.gv_data['rows'][:2], KSTAR_NAMES):
            totals[name] = gv**2*sum(gv_fit.width_coefficient(
                row['physical_vector_mass_GeV'], mk, physical[first], physical[second],
                self.kaon_flavors[mode]*factors[mode]**2) for first, second, mode in KSTAR_MODES[name])
        if min(totals.values()) <= 0:
            raise ValueError('nonpositive total width')
        return point, sector, mk, gv, totals

    def old_blocks(self, values, state):
        """Inherited 28 rows; the two phi -> K Kbar rows get the squared phi factors."""
        residual, odd_jacobian, details = super().old_blocks(values, state)
        coefficients = self.coefficients(values)
        if coefficients is None:
            return residual, odd_jacobian, details
        factors = self.factors(coefficients, np.real(state[1].generators[2]))
        predictions = [prediction*factors[mode]**2
                       for prediction, mode in zip(details['phi_branching_predictions'], PHI_MODES)]
        phi_residual = np.linalg.solve(self.phi_cholesky, np.array(predictions)-self.phi_target)
        return (np.r_[residual[:26], phi_residual], odd_jacobian,
                dict(details, phi_branching_predictions=predictions))

    def _babar_blocks_full(self, values, state, derivatives=False, channels=CHANNELS):
        """BaBar K*-channel terms with the K* -> K pi mode factor of each vertex.

        The baseline propagator cache is keyed on the first 12 and the four
        pole-residue coordinates only, so it is cleared whenever r changes.
        The mode is read off the invariant: the integrated t of each channel
        is identified by identity with self.phase[channel]['t']; the next K*
        call of the same channel is its u (baseline order kt, then ku).
        """
        if not self.isospin_breaking_vpp:
            return super()._babar_blocks_full(values, state, derivatives, channels)
        added = tuple(values[BASELINE_COORDINATE_COUNT:])
        if added != self.last_added_coordinates:
            self.babar_cache_key = None
            self.last_added_coordinates = added
        coefficients = self.coefficients(values)
        if coefficients is None:
            return super()._babar_blocks_full(values, state, derivatives, channels)
        factors = self.factors(coefficients)
        u_modes = {'kkpi0': 'kstar_charged_K+_pi0', 'kskpi': 'kstar_charged_K0_pi+'}
        t_modes = {'kkpi0': 'kstar_charged_K+_pi0', 'kskpi': 'kstar_neutral_K+_pi-'}
        current_channel = []

        def mode_scaled(name, invariant, beta, mass, shared_width, reference):
            full = RUNNING_PROPAGATOR(name, invariant, beta, mass, shared_width, reference)
            if name != 'kstar':
                return full
            channel = next((c for c in CHANNELS if invariant is self.phase[c]['t']), None)
            if channel is not None:
                current_channel[:] = [channel]
                mode = t_modes[channel]
            elif current_channel:
                mode = u_modes[current_channel[0]]
            else:
                raise RuntimeError('K* propagator called before the t invariant of its channel')
            factor = factors[mode]
            if self.excitation_scaling == 'whole_propagator':
                return factor*full
            ground = RUNNING_PROPAGATOR(name, invariant, (0., 0.), mass, shared_width, reference)
            return full+(factor-1)*ground

        joint_fit.running_propagator = mode_scaled
        try:
            return super()._babar_blocks_full(values, state, derivatives, channels)
        finally:
            joint_fit.running_propagator = RUNNING_PROPAGATOR

    def residual(self, values):
        """Inherited residual after a check of the coordinate count."""
        self.check_length(values)
        return super().residual(values)

    def jacobian(self, values, finite_difference_steps=None):
        """Inherited Jacobian plus forward differences for r1, r2, r3."""
        self.check_length(values)
        baseline_jacobian = super().jacobian(values, finite_difference_steps)
        if not self.added_coordinate_names:
            return baseline_jacobian
        reference = self.residual(values)
        added_columns = []
        for index in range(BASELINE_COORDINATE_COUNT, len(values)):
            step = (np.sqrt(np.finfo(float).eps)*max(1., abs(values[index]))
                    if finite_difference_steps is None else finite_difference_steps[index])
            shifted = np.array(values, dtype=float)
            shifted[index] += step
            added_columns.append((self.residual(shifted)-reference)/step)
        return np.column_stack([baseline_jacobian, *added_columns])


def main():
    model = IsospinBreakingVppLikelihood(1, 3, 4, bin_block_size=4, isospin_breaking_vpp=True)
    lower, upper = model.parameter_bounds()
    start = model.initial()
    report = dict(
        coordinate_count=len(model.parameter_names), nominal_degrees_of_freedom=75-len(model.parameter_names),
        added_coordinates=[dict(name=name, start=float(start[index]), lower=float(lower[index]),
                                upper=float(upper[index]))
                           for index, name in enumerate(model.parameter_names)
                           if index >= BASELINE_COORDINATE_COUNT],
        spurion_diagonal_gev2=model.spurion_diagonal,
        kstar_mode_factors_at_scan_point=model.factors(SCAN_POINT_GEV_MINUS2))
    print(json.dumps(report, indent=1))


if __name__ == '__main__':
    main()
