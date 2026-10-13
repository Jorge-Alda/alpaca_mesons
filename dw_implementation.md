Implement numerically in Python the decay widths of ALPs into 3 pseudoscalar mesons that were derived in `/memory/alp-three-meson-amplitudes-widths.md`. The result will be saved to the folder `/scripts/implementation`.

The inputs of the function `decay_width` are:

- ALP mass, in GeV
- Coupling matrices $\hat{c}_q^A$ and $\hat{c}_q^P$ (the user already provides these in the basis where there is no meson mixing). These are $3\times 3$ hermitian numpy arrays.
- Three strings identifying the final mesons
- $f_a$, in GeV.

If the decay is not kinematically allowed, doesn't conserve charge or strangeness, the function must return 0.0 without doing any computation.

If two or three mesons are identical, use the formulas that require numerical integration only in one Mandelstam variable.

Use the conventions in `/memory/sources/Alda:2026nrd.pdf`.

Use the existing virtual environment in `.venv`