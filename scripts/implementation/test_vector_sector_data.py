"""The calibrated sector reaches pinned observations without changing errors."""
import unittest
import numpy as np
from .vector_sector import VectorSector
from .vector_sector_amplitudes import VectorAmplitudeInputs
from .pseudoscalar_tff import TFFInputs
from .radiative_sm import Q
from .odd_parity_data import (load_cleo_tff_observations,load_babar_pi0_tff_observation,
    load_vector_omega_pi0_gamma_observation,load_omega_pi0_gamma_observation)
from .radiative_odd_parity import OmegaPiGammaInputs
from .alda_2026nrd_neutral_inputs import alda_2026nrd_neutral_inputs


class VectorSectorDataTests(unittest.TestCase):
    def test_calibrated_sector_reaches_all_tff_observation_blocks(self):
        sector=VectorSector.from_stage3(rho_input='rho_neutral',kstar_input='kstar_neutral',
            phi_input='phi_evaluation',source_hypothesis='double_trace')
        data=alda_2026nrd_neutral_inputs(mass_scheme='physical')
        selected={name:VectorAmplitudeInputs(p['F'],.77,p['mp'],p['TP'],p['chi'],sector)
                  for name,p in data.items()}
        observations=load_cleo_tff_observations({'pi0':selected['pi0'],'eta':selected['eta'],
                                               'eta_prime':selected['etap']})
        self.assertEqual(sum(o.values.size for o in observations),63)
        observations+= (load_babar_pi0_tff_observation(selected['pi0']),
                        load_vector_omega_pi0_gamma_observation(selected['pi0'],.3028,.00868))
        parameters={f'c{i}':.001*i for i in range(1,8)}
        parameters.update({f'd{i}':-.002*i for i in range(1,5)})
        for observation in observations:
            prediction=observation.predict(parameters,observation.config)
            self.assertEqual(prediction.shape,observation.values.shape)
            self.assertTrue(np.all(np.isfinite(prediction)))

    def test_generalized_factories_recover_baseline_likelihood_and_predictions(self):
        sector=VectorSector.universal(.77,.15)
        TP=np.diag([1.,-1.,0.])/np.sqrt(2);chi=np.diag([.135**2,.135**2,.47])
        x=VectorAmplitudeInputs(.093,.77,.135,TP,chi,sector)
        fq=np.trace((Q@Q+Q@Q)@TP).real
        fh=np.trace((Q@Q+Q@Q)@(TP@chi+chi@TP)).real
        old=TFFInputs(x.F,.15,x.MV,x.mp,fq,fh)
        parameters={f'c{i}':.001*i for i in range(1,8)}
        parameters.update({f'd{i}':-.002*i for i in range(1,5)})
        pairs=[(load_cleo_tff_observations({'pi0':old})[0],load_cleo_tff_observations({'pi0':x})[0]),
               (load_babar_pi0_tff_observation(old),load_babar_pi0_tff_observation(x)),
               (load_omega_pi0_gamma_observation(OmegaPiGammaInputs(.3,.093,.15,.77,.77,.135),.00868),
                load_vector_omega_pi0_gamma_observation(x,.3,.00868))]
        for baseline,current in pairs:
            np.testing.assert_array_equal(current.values,baseline.values)
            np.testing.assert_array_equal(current.covariance,baseline.covariance)
            self.assertEqual(current.nuisances,baseline.nuisances)
            np.testing.assert_allclose(current.predict(parameters,current.config),
                baseline.predict(parameters,baseline.config),rtol=1e-12,atol=1e-14)


if __name__=='__main__':unittest.main()
