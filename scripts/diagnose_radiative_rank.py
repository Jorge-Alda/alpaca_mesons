"""Action-based radiative rank probe, with explicit mass and field assumptions.

Run: .venv/bin/python scripts/diagnose_radiative_rank.py
The previous provisional rows were replaced after reconstructing the action.
See memory/rcht/odd-parity/radiative-sm-amplitude.md for scope and provenance.
"""
from derive_radiative_sm_amplitude import rank_probe

if __name__=='__main__':
    rank_probe()
