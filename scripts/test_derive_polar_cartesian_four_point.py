"""Checks for the field expansion and mixed-propagator cancellation."""
import unittest
import sympy as sp
from scripts.derive_polar_cartesian_four_point import derive


class PolarCartesianEquivalence(unittest.TestCase):
    def test_channel_cancellation(self):
        self.assertEqual(derive()[-1], sp.ZeroMatrix(1, 1))

    def test_inverse_map(self):
        radial, phase, *_ = derive()
        symbols = {str(symbol): symbol for symbol in radial.free_symbols | phase.free_symbols}
        s, p, f = (symbols[name] for name in ('s', 'p', 'f_varphi'))
        t = sp.Symbol('t', real=True)
        radius = radial.subs({s: t*s, p: t*p})
        angle = phase.subs({s: t*s, p: t*p})
        real_part = (f + radius)*sp.cos(angle/f) - f
        imaginary_part = (f + radius)*sp.sin(angle/f)
        self.assertEqual(sp.simplify(sp.series(real_part, t, 0, 4).removeO() - t*s), 0)
        self.assertEqual(sp.simplify(sp.series(imaginary_part, t, 0, 4).removeO() - t*p), 0)


if __name__ == '__main__':
    unittest.main()
