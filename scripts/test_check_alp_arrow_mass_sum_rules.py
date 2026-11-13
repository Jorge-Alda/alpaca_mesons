"""Run: .venv/bin/python -m unittest discover -s scripts -p test_check_alp_arrow_mass_sum_rules.py"""
import unittest
from check_alp_arrow_mass_sum_rules import check_sum_rules


class ArrowMassSumRulesTest(unittest.TestCase):
    def test_exact_pole_sum_and_product(self):
        for description, passed in check_sum_rules().items():
            with self.subTest(description=description):
                self.assertTrue(passed)


if __name__ == '__main__':
    unittest.main()
