"""Isospin-breaking VPP likelihood with reuse of repeated calculations (bitwise-identical outputs).

Run (prints the bitwise comparison with IsospinBreakingVppLikelihood at the
restart endpoint and at r = 0, cheap quadrature by default):
  .venv/bin/python scripts/joint_vector_isospin_breaking_vpp_reuse.py
Tests: .venv/bin/python -m unittest scripts.test_joint_vector_isospin_breaking_vpp_reuse

No change of numerics. IsospinBreakingVppLikelihood
(joint_vector_isospin_breaking_vpp_option.py, not modified) repeats three
calculations; the subclass below performs each once and returns the same
floating-point numbers:

1. Zero-excitation running propagator. With the option on, every K* call and
   the s-channel rho call of mode_scaled is followed by a second call
   RUNNING_PROPAGATOR(name, invariant, (0, 0), mass, shared_width, reference)
   at the same invariant. In the reference propagator
   (fit_babar_kkpi_g_coefficients.propagator) the ground pole is computed once
   on the invariant and every excitation choice adds beta * (excited poles), so
   with beta = 0 both choices equal the ground pole of the first call (for
   finite excited poles). The second call is therefore
       breit_wigner + (ground - ground),
   with breit_wigner = 1/(m^2 - invariant - i m Gamma(invariant)) and ground
   from the first call. GroundPropagatorMemory keeps these two arrays from the
   first call and returns that expression for the second call when the name,
   the invariant and reference objects (identity) and the mass and shared width
   (equality) agree; otherwise it evaluates the full expression of
   fit_babar_kkpi_global_vector_sector.running_propagator, written out
   operation by operation.
2. Duplicate reference residual in the Jacobian. The option's jacobian calls
   the inherited Jacobian, whose first step is residual(values), and then
   residual(values) again as the reference for the r1, r2, r3 forward
   differences. Here the first result is kept and used as that reference.
3. Vertex factors. T(r)/T(0) of the K* and phi modes (vertex_factors) and of
   the s-channel modes (s_channel_factors) are recomputed for every energy-bin
   block. Here they are kept per (r, flavor matrices) key; the dictionary
   handed out is a fresh copy because the option updates it in place.

Equality is not assumed: the unit tests compare IEEE-754 bit patterns of the
residual and the Jacobian with IsospinBreakingVppLikelihood at the restart
endpoint, at r = 0 and at the scan point, at the fit quadrature (2,16,16).

Units and conventions: as in joint_vector_isospin_breaking_vpp_option.py (GeV,
GeV^2, r in GeV^-2; Alda:2026nrd).
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fit_babar_kkpi_global_vector_sector as joint_fit  # noqa: E402
import joint_vector_isospin_breaking_vpp_option as vpp_option  # noqa: E402
from fit_babar_kkpi_g_coefficients import CHANNELS  # noqa: E402

ZERO_EXCITATION = (0., 0.)
FACTOR_MEMORY_SIZE = 8
S_CHANNEL_FACTORS = vpp_option.s_channel_factors


class GroundPropagatorMemory:
    """running_propagator that reuses its ground pole for the following zero-excitation call."""

    def __init__(self):
        self.stored = None
        self.reused_calls = 0
        self.evaluated_calls = 0

    def __call__(self, name, invariant, beta, mass, shared_width, reference):
        stored = self.stored
        if (stored is not None and isinstance(beta, tuple) and beta == ZERO_EXCITATION
                and stored[0] == name and stored[1] is invariant and stored[2] == mass
                and stored[3] == shared_width and stored[4] is reference):
            self.stored = None
            self.reused_calls += 1
            breit_wigner, ground = stored[5], stored[6]
            return breit_wigner+(ground-ground)
        self.evaluated_calls += 1
        # Same operations, in the same order, as joint_fit.running_propagator.
        ground, reference_with_excitations = joint_fit.reference_propagator_pair(name, invariant, beta, reference)
        denominator = 1/ground
        on_shell = 1/joint_fit.reference_propagator(name, mass**2, reference, 0.)
        width = -denominator.imag/mass*shared_width/(-on_shell.imag/mass)
        excited = reference_with_excitations-ground
        breit_wigner = 1/(mass**2-invariant-1j*mass*width)
        self.stored = (name, invariant, mass, shared_width, reference, breit_wigner, ground)
        return breit_wigner+excited


def factor_key(coefficients, *matrices):
    """Hashable key of r and the flavor matrices (exact bytes)."""
    return (np.asarray(coefficients, dtype=float).tobytes(),
            *(None if matrix is None else np.asarray(matrix).tobytes() for matrix in matrices))


class IsospinBreakingVppLikelihoodWithReuse(vpp_option.IsospinBreakingVppLikelihood):
    """IsospinBreakingVppLikelihood with the three reuses of the module docstring.

    Same constructor, coordinates, bounds and outputs as the parent.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.factor_memory = {}
        self.jacobian_point = None
        self.jacobian_reference = None

    def remembered(self, key, compute):
        """Fresh copy of a remembered factor dictionary (computed on first use)."""
        if key not in self.factor_memory:
            if len(self.factor_memory) >= FACTOR_MEMORY_SIZE:
                self.factor_memory.clear()
            self.factor_memory[key] = compute()
        return dict(self.factor_memory[key])

    def factors(self, coefficients, phi_generator=None):
        """Remembered per-mode T(r)/T(0) (parent method on first use)."""
        return self.remembered(('modes',)+factor_key(coefficients, phi_generator),
                               lambda: super(IsospinBreakingVppLikelihoodWithReuse, self).factors(
                                   coefficients, phi_generator))

    def s_channel_factors(self, coefficients, spurion_diagonal, omega_generator, phi_generator):
        """Remembered vpp_option.s_channel_factors."""
        return self.remembered(
            ('s_channel',)+factor_key(coefficients, spurion_diagonal, omega_generator, phi_generator),
            lambda: S_CHANNEL_FACTORS(coefficients, spurion_diagonal, omega_generator, phi_generator))

    def _babar_blocks_full(self, values, state, derivatives=False, channels=CHANNELS):
        """Parent method with the propagator memory and remembered s-channel factors in place."""
        if not self.isospin_breaking_vpp:
            return super()._babar_blocks_full(values, state, derivatives, channels)
        original_running = vpp_option.RUNNING_PROPAGATOR
        original_s_channel = vpp_option.s_channel_factors
        vpp_option.RUNNING_PROPAGATOR = GroundPropagatorMemory()
        vpp_option.s_channel_factors = self.s_channel_factors
        try:
            return super()._babar_blocks_full(values, state, derivatives, channels)
        finally:
            vpp_option.RUNNING_PROPAGATOR = original_running
            vpp_option.s_channel_factors = original_s_channel
            joint_fit.running_propagator = original_running

    def residual(self, values):
        """Parent residual; kept when it is the first one at the current Jacobian point."""
        result = super().residual(values)
        if values is self.jacobian_point and self.jacobian_reference is None:
            self.jacobian_reference = result.copy()
        return result

    def jacobian(self, values, finite_difference_steps=None):
        """Parent Jacobian; the reference residual of the r columns is the inherited one."""
        self.check_length(values)
        self.jacobian_point, self.jacobian_reference = values, None
        try:
            baseline_jacobian = super(vpp_option.IsospinBreakingVppLikelihood, self).jacobian(
                values, finite_difference_steps)
            reference = self.jacobian_reference
        finally:
            self.jacobian_point, self.jacobian_reference = None, None
        if not self.added_coordinate_names:
            return baseline_jacobian
        if reference is None:
            reference = self.residual(values)
        added_columns = []
        for index in range(vpp_option.BASELINE_COORDINATE_COUNT, len(values)):
            step = (np.sqrt(np.finfo(float).eps)*max(1., abs(values[index]))
                    if finite_difference_steps is None else finite_difference_steps[index])
            shifted = np.array(values, dtype=float)
            shifted[index] += step
            added_columns.append((self.residual(shifted)-reference)/step)
        return np.column_stack([baseline_jacobian, *added_columns])


def bitwise_comparison(orders=(1, 3, 4), ledger='data/vector_sector/joint_vector_isospin_breaking_vpp_fit.json'):
    """Bit-pattern equality of residual and Jacobian with the parent class at three points."""
    options = dict(bin_block_size=4, isospin_breaking_vpp=True)
    plain = vpp_option.IsospinBreakingVppLikelihood(*orders, **options)
    reuse = IsospinBreakingVppLikelihoodWithReuse(*orders, **options)
    runs = json.loads((joint_fit.ROOT/ledger).read_text())['entries'][0]['runs']
    endpoint = np.asarray(min((run for run in runs if run.get('high_order', {}).get('valid')),
                              key=lambda run: run['high_order']['chi2'])['endpoint'])
    points = dict(restart=endpoint, r_zero=np.r_[endpoint[:32], 0., 0., 0.],
                  scan_point=np.r_[endpoint[:32], vpp_option.SCAN_POINT_GEV_MINUS2])
    result = {}
    for name, point in points.items():
        result[name] = {method: bool(np.array_equal(getattr(plain, method)(point).view(np.uint64),
                                                    getattr(reuse, method)(point).view(np.uint64)))
                        for method in ('residual', 'jacobian')}
    return result


if __name__ == '__main__':
    print(json.dumps(bitwise_comparison(), indent=1))
