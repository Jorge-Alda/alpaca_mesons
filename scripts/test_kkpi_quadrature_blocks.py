"""Energy-bin block equivalence; GeV and nb conventions.

Run: .venv/bin/python -m unittest scripts.test_kkpi_quadrature_blocks
"""
import unittest
from unittest.mock import patch

import numpy as np
from scripts import fit_babar_kkpi_global_vector_sector as fit


class EnergyBinBlockTests(unittest.TestCase):
    def test_block_boundaries_preserve_every_channel_bin(self):
        model=fit.GlobalLikelihood(1,2,2)
        blocks=list(fit.EnergyBinBlocks(model.babar,(2,16,16),7).blocks())
        for channel in fit.CHANNELS:
            edges=np.concatenate([block[channel]['edges'] for block in blocks])
            np.testing.assert_array_equal(edges,model.babar[channel]['edges'])
            self.assertTrue(all(len(block[channel]['edges'])<=7 for block in blocks))
            self.assertEqual(sum(len(block[channel]['errors']) for block in blocks),len(edges))

    def test_predictions_residuals_and_jacobians_match_full_grid(self):
        for orders in ((1,3,4),(2,16,16)):
            with self.subTest(orders=orders):
                full=fit.GlobalLikelihood(*orders)
                values=full.initial();values[20:26]=[.001,-.002,.003,-.001,.002,.001]
                expected_residual=full.residual(values)
                _,expected_predictions=full.babar_blocks(values,full.unpack(values))
                expected_jacobian=full.jacobian(values)
                for size in (1,7,100):
                    blocked=fit.GlobalLikelihood(*orders,bin_block_size=size)
                    residual=blocked.residual(values)
                    _,predictions=blocked.babar_blocks(values,blocked.unpack(values))
                    np.testing.assert_array_equal(residual,expected_residual)
                    for channel in fit.CHANNELS:
                        np.testing.assert_array_equal(predictions[channel],expected_predictions[channel])
                    np.testing.assert_array_equal(blocked.jacobian(values),expected_jacobian)
                    np.testing.assert_array_equal(blocked.residual(values),residual)
                    self.assertIsInstance(blocked.phase,fit.EnergyBinBlocks)
                    self.assertIsNone(blocked.babar_cache_key)

    def test_order_128_maps_are_built_one_bin_at_a_time(self):
        model=fit.GlobalLikelihood(1,2,2,bin_block_size=1)
        descriptor=fit.EnergyBinBlocks(model.babar,(2,128,128),1)
        first=next(descriptor.blocks())
        mapped=fit.phase_points(first,*descriptor.orders)
        for channel in fit.CHANNELS:
            phase=mapped[channel]
            self.assertEqual(phase['shape'],(1,2,384,256))
            self.assertTrue(np.all(np.isfinite(phase['weight'])))
            self.assertTrue(np.all(phase['weight']>=0))
            masses=(.493677,.493677,.1349768) if channel=='kkpi0' else (.497611,.493677,.13957039)
            first_mass,second_mass,third_mass=masses
            self.assertTrue(np.all(phase['s']>(first_mass+second_mass)**2))
            self.assertTrue(np.all(phase['s']<(np.sqrt(phase['q'])-third_mass)**2))
            outer=fit.triangle(phase['q'],phase['s'],third_mass**2)
            inner=fit.triangle(phase['s'],first_mass**2,second_mass**2)
            span=np.sqrt(outer*inner)/phase['s']
            other_center=(first_mass**2+third_mass**2+(phase['s']+first_mass**2-second_mass**2)
                *(phase['q']-phase['s']-third_mass**2)/(2*phase['s']))
            center=phase['q']+sum(mass**2 for mass in masses)-phase['s']-other_center
            self.assertTrue(np.all(phase['t']>center-span/2))
            self.assertTrue(np.all(phase['t']<center+span/2))

    def test_grid_builder_never_receives_more_than_one_block(self):
        model=fit.GlobalLikelihood(1,3,4,bin_block_size=7)
        original=fit.phase_points;counts=[]
        def build(data,*orders):
            counts.extend(len(data[channel]['edges']) for channel in fit.CHANNELS)
            return original(data,*orders)
        with patch.object(fit,'phase_points',side_effect=build):
            model.residual(model.initial())
        self.assertTrue(counts)
        self.assertLessEqual(max(counts),7)


if __name__=='__main__':unittest.main()
