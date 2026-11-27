import math
import numpy as np
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_up_only_eta_prime_heaviest as check
from only_up_inputs import ETA_PRIME_MASS, M1_SQUARED, B_ZERO, UP_MASS, DOWN_MASS, STRANGE_MASS


class EtaPrimeHeaviestTest(unittest.TestCase):
    def test_epsilon_kappa_roundtrip_and_grid_range(self):
        grid = check.default_epsilon_grid()
        self.assertAlmostEqual(grid[0], 0.1)
        self.assertAlmostEqual(grid[-1], 150.)
        for epsilon in grid:
            self.assertAlmostEqual(check.epsilon_from_kappa(check.kappa_from_epsilon(epsilon)) / epsilon, 1.)

    def test_singlet_parameter_reproduces_input_mass(self):
        diagonal = M1_SQUARED + 2 * B_ZERO * (UP_MASS + DOWN_MASS + STRANGE_MASS) / 3
        self.assertAlmostEqual(math.sqrt(diagonal), ETA_PRIME_MASS, places=12)
        self.assertAlmostEqual(M1_SQUARED, 0.7295546296, places=9)

    def test_small_mixing_limit_is_unmixed_eta_prime(self):
        # At epsilon -> 0 the pole is the largest eigenvalue of the 3x3 (pi0, eta, eta') block,
        # which differs from the input diagonal entry by eta-eta' mixing from m_d != m_u.
        row = check.scan([1e-3])[0]
        mass_squared = np.zeros((3, 3))
        generators = [np.diag([1., -1., 0.]), np.diag([1., 1., -2.]) / math.sqrt(3), np.eye(3) * math.sqrt(2 / 3)]
        chi = 2 * B_ZERO * np.diag([UP_MASS, DOWN_MASS, STRANGE_MASS])
        for i, first in enumerate(generators):
            for j, second in enumerate(generators):
                mass_squared[i, j] = (np.trace((first @ second + second @ first) @ chi) / 4
                                      + M1_SQUARED * np.trace(first) * np.trace(second) / 6)
        self.assertAlmostEqual(row["eta_prime_mass"], math.sqrt(np.linalg.eigvalsh(mass_squared)[-1]), places=6)
        self.assertTrue(row["eta_prime_is_heaviest"])

    def test_scan_is_deterministic_and_complete(self):
        grid = check.default_epsilon_grid(7)
        first, second = check.scan(grid), check.scan(grid)
        self.assertEqual(first, second)
        for row in first:
            self.assertGreater(row["eta_prime_mass"], 0.)
            self.assertIn(row["heaviest_label"], check.MESON_LABELS)


if __name__ == "__main__":
    unittest.main()
