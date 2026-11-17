"""Compatibility import for the finite-mixing three-pseudoscalar calculator."""
import numpy as np
from calculate_alp_meson_masses_finite import calculate_masses


def _matrix(value, name):
    result = np.asarray(value, dtype=complex)
    if result.shape != (3, 3) or not np.all(np.isfinite(result)):
        raise ValueError(f'{name} must be a finite 3 by 3 matrix')
    if not np.allclose(result, result.conj().T, rtol=1e-12, atol=1e-14):
        raise ValueError(f'{name} must be Hermitian')
    return result
