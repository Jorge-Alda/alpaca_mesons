"""Corrected TFF tensor derivation, including Fourier signs and both photons.

The previous one-assignment rank claim was invalid. The full calculation lives
in derive_radiative_tff_amplitude.py. Keep the historical invocation functional:
    .venv/bin/python scripts/derive_radiative_tff_kinematics.py
"""
from derive_radiative_tff_amplitude import main

if __name__ == '__main__':
    main()
