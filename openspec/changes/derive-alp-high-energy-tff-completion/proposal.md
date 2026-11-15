# Proposal

## Why

The doubly virtual ALP transition form factor (TFF) needed for $a\to\ell^+\ell^-$ is integrated over unbounded photon virtualities, so the low-energy hadronic representation must be completed by an ALP-specific short-distance behavior. Only the leading local hard-current coefficient for nearly equal, deeply Euclidean virtualities exists in `memory/alp/lepton-pair/alp-hard-current-ope.md`. Its renormalized ALP matrix elements, other virtuality limits, finite-$m_a$ corrections, flavor thresholds, and joining to the low-energy form are open (beads issue `alp-6os.2`, analytical part).

## What Changes

- Express the large-virtuality TFF in the no-mixing source basis of Alda:2026nrd through renormalized ALP matrix elements of quark axial currents and, where required, the gluonic operator, with the scale, scheme, and active-flavor thresholds stated.
- Derive the short-distance behavior in three regimes: both virtualities large with fixed ratio, one virtuality large with the other fixed (light-cone expansion), and the nearly symmetric limit already obtained; state which regime each result covers.
- Derive finite-$m_a$ corrections of relative order $m_a^2/Q^2$ and the heavy-quark threshold treatment ($c$, $b$) for the supported ALP mass range $2m_e<m_a<2m_\tau$.
- Derive the ALP-specific asymptotic coefficient without importing the $\eta$ or $\eta'$ value, and record what it depends on (source matrices, decay-constant-like matrix elements, anomaly/singlet treatment).
- Derive a matching prescription joining the asymptotic form to the low-energy representation (poles plus cuts) at a declared matching scale, with no double counting of strength already in the low-energy form, and compare alternative completions (effective-pole and interpolating forms) to quantify model dependence.
- Provide executable SymPy derivations under `scripts/` and paired indexed memory records. Numerical matching to the vector-fit-calibrated low-energy form is out of scope here and remains conditional on `alp-6os.1`.

## Capabilities

### New Capabilities
- `alp-high-energy-tff`: Source-resolved, ALP-specific short-distance behavior of the doubly virtual TFF, its scheme and threshold provenance, and its matching to the low-energy representation.

### Modified Capabilities
- None. The existing `low-energy-tff-spectra` change is an unarchived delta, and this change only consumes its low-energy representation.

## Impact

New records under `memory/alp/lepton-pair/` (with updated indexes), new derivation scripts under `scripts/` (written by the `coder` agent), and bibliography entries for consulted sources. Feeds `alp-6os.3` (library) and `alp-6os.4` ($a\to\ell^+\ell^-$). No existing script or record is changed except for index links and a status note in the input-gates record.
