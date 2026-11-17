#!/usr/bin/env python3
"""Inputs for the only-up phase model of DiLuzio:2026poh.

All masses and decay constants are in GeV and squared masses in GeV^2.
The scalar and pseudoscalar lists contain derivatives in x=a/f_a, not
Taylor coefficients. The radial field is set to zero. m1_squared is the
singlet anomaly mass parameter, fixed from the supplied eta-prime input.

"""
import math

import numpy as np

PION_MASS = 139.57061e-3
PION_DECAY_CONSTANT = 92.3e-3
UP_MASS = 2.20e-3
DOWN_MASS = 4.69e-3
STRANGE_MASS = 92.74e-3
B_ZERO = PION_MASS**2 / (UP_MASS + DOWN_MASS)
ETA_PRIME_MASS = 957.78e-3
M1_SQUARED = ETA_PRIME_MASS**2 - 2 * B_ZERO * (UP_MASS + DOWN_MASS + STRANGE_MASS) / 3


def only_up_inputs(kappa_u):
    """Return the phase sources and fixed singlet parameter for kappa_u > 0."""
    if not math.isfinite(kappa_u) or kappa_u <= 0:
        raise ValueError('kappa_u must be positive and finite; zero is a limiting case')
    zeros = lambda: np.zeros((3, 3), dtype=float)
    scalar = [np.diag([UP_MASS, DOWN_MASS, STRANGE_MASS]), zeros(),
              np.diag([-UP_MASS, 0., 0.]), zeros(),
              np.diag([UP_MASS, 0., 0.])]
    pseudoscalar = [zeros(), np.diag([-UP_MASS, 0., 0.]), zeros(),
                    np.diag([UP_MASS, 0., 0.])]
    return {'f_a': math.sqrt(2)*UP_MASS/kappa_u, 'm_a0': 0.,
            'hat_S': scalar, 'hat_P': pseudoscalar, 'hat_Theta': [0.]*4,
            'f_pi': PION_DECAY_CONSTANT, 'B_0': B_ZERO, 'm_0_squared': M1_SQUARED}

