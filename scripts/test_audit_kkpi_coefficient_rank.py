"""Focused independent algebraic checks for the KKpi rank audit.

Run: .venv/bin/python -m unittest scripts.test_audit_kkpi_coefficient_rank
"""

import unittest
import sympy as sp

from scripts.audit_kkpi_coefficient_rank import coefficient_matrix, derive_maps


class KKPiRankTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.results = derive_maps()

    def test_vppp_recovery_and_su3_limit(self):
        r = self.results
        g = r['g']
        matrix = r['maps']['VPPP both plus isoscalar']
        self.assertEqual(matrix.rank(),4)
        self.assertEqual(matrix.nullspace(),[sp.Matrix([1,0,1,0,0])])
        self.assertEqual(matrix.subs(r['pion'],r['kaon']).rank(),3)
        self.assertEqual(r['maps']['VPPP neutral plus isoscalar'].rank(),3)

    def test_local_energy_information(self):
        r = self.results
        q,s,t = r['invariants']
        fixed = coefficient_matrix([r['local_charged'].subs(q,sp.Rational(2))],r['g'],(s,t))
        self.assertEqual(fixed.rank(),2)
        self.assertEqual(r['maps']['VPPP charged'].rank(),3)

    def test_mass_direction_after_isoscalar_constraints(self):
        r = self.results
        g = r['g']
        # A fixed, g2 fixed and B fixed imply dg1=dg2=dg3=0,
        # dg4=1 and dg5=-2. Its D response must vanish in SU(3).
        response = sp.diff(r['local_charged'],g[3])-2*sp.diff(r['local_charged'],g[4])
        self.assertEqual(sp.expand(response),2*(r['pion']-r['kaon']))

    def test_full_amplitudes_with_preliminary_span(self):
        maps = self.results['maps']
        self.assertEqual(maps['preliminary coefficient span'].rank(),7)
        self.assertEqual(maps['full neutral plus preliminary'].rank(),9)
        self.assertEqual(maps['full charged plus preliminary'].rank(),11)
        self.assertEqual(maps['full both plus preliminary plus isoscalar'].rank(),12)
        # Full nulls are verified by multiplication, independently of rank.
        for label in ('full charged','full neutral','full both plus preliminary plus isoscalar'):
            for vector in maps[label].nullspace():
                self.assertEqual(sp.simplify(maps[label]*vector),sp.zeros(maps[label].rows,1))

    def test_cross_sector_compensation_for_symbolic_masses(self):
        self.assertEqual(self.results['compensation_responses'],[0,0])


if __name__ == '__main__':
    unittest.main()
