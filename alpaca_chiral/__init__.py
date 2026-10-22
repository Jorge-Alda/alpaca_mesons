"""ALP chiral decay amplitudes and widths, with masses and widths in GeV."""
from .pseudoscalar_diag import ALPcouplings, cA_no_mixing, cP_no_mixing, cG_no_mixing
from .three_pseudoscalar import three_pseudoscalar_amplitude, three_pseudoscalar_width
from .vector_pair import neutral_vector_pair_amplitude, neutral_vector_pair_width
from .vector_photon import vector_photon_amplitude, vector_photon_width
from .vector_two_pseudoscalar import vector_two_pseudoscalar_amplitude, vector_two_pseudoscalar_width
from .photon_two_pseudoscalar import photon_two_pseudoscalar_amplitude, photon_two_pseudoscalar_width
from .two_photon import photon_coupling, two_photon_width

__all__ = [
    'ALPcouplings', 'cA_no_mixing', 'cP_no_mixing', 'cG_no_mixing',
    'three_pseudoscalar_amplitude', 'three_pseudoscalar_width',
    'neutral_vector_pair_amplitude', 'neutral_vector_pair_width',
    'vector_photon_amplitude', 'vector_photon_width',
    'vector_two_pseudoscalar_amplitude', 'vector_two_pseudoscalar_width',
    'photon_two_pseudoscalar_amplitude', 'photon_two_pseudoscalar_width',
    'photon_coupling', 'two_photon_width',
]
