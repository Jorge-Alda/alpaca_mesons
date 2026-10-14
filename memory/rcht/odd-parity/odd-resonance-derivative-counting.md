---
title: Hard vector momenta and residual derivative counting in the odd action
status: mixed
created: 2026-10-14
updated: 2026-10-14
bibtex_keys: [Alda:2026nrd, Kaiser:2000gs, Ruiz-Femenia:2003jdx]
provenance: analytical-derivation
---

# Inputs and selected prescription

Continue [the power-counting proposal](odd-action-power-counting-proposal.md)
and [its parent filtration](odd-hilbert-power-counting.md). Use the
[canonical two-boson contraction](odd-higher-two-boson-combinations.md) and
[explicit-vector convention](odd-covariant-closure.md). No vector EOM,
propagator expansion, quark-mass substitution, new coefficient prior or
loop calculation is introduced. Units are GeV; scalar virtualities have
units GeV squared. Applying relative Goldstone/OZI counting to explicit
resonances remains the proposal's assumption [@Kaiser:2000gs].

The reproducible calculation is
[derive_odd_power_counted_channels.py](../../../scripts/derive_odd_power_counted_channels.py).
Its companion [transcript](../../../data/vector_sector/odd_power_counted_channels.md)
derives the split below with SymPy, rather than assigning a derivative on an
external vector a soft momentum by definition.

# Exact hard/residual split

With $s=q_A^2$, $t=q_B^2$ and $h=q_A\cdot q_B$, every retained extra-two-
derivative cubic image is linear in these invariants. For
$P=k_ss+k_tt+k_hh$, substitute

$$
s=s_0+\delta_s,\quad t=t_0+\delta_t,\quad h=h_0+\delta_h.
$$

The script produces the exact decomposition

$$
P=P_{\rm hard}+P_{\rm residual},\qquad
P_{\rm hard}=k_ss_0+k_tt_0+k_hh_0,
\qquad P_{\rm residual}=k_s\delta_s+k_t\delta_t+k_h\delta_h.
$$

Only the residual term can be assigned suppression by small variations
around that reference. The hard term is still tied to the same parent
coefficient. It is not discarded, fitted independently or absorbed into
an unrelated leading source coefficient. The $1/M_V^2$ normalization of
a higher derivative operator does not suppress $M_i^2/M_V^2$ automatically.

# Channel prescriptions

For a genuinely soft low-virtuality photon/pseudoscalar domain, use the
existing formal cost $c=\ell+n_{\rm tr}-1$, where $\ell=1$ for higher seeds.
This assumes all contracted invariants are soft compared with the declared
$\Lambda^2$. It is conditional, not a numerical domain established by data.

For channels containing on-shell vectors, keep $s=M_1^2$ and, if appropriate,
$t=M_2^2$ exactly. For a real photon $t=0$, but
$h=(m_a^2-M_1^2)/2$ need not be soft. The $Ydd$ and $fmdd$ parents therefore
receive no automatic extra-two-derivative suppression. Conservatively promote
their parent cost to $n_{\rm tr}-1$, while retaining the ordinary additional
spurion penalty for the $X$-dressed derivative parents. Filter parents before
their images, including parents whose images subsequently vanish.

This promotion is a declared resonance-region bookkeeping prescription, not
a proof of a controlled heavy-resonance expansion. It retains all derivative
mass-source blocks at the relevant trace order. A more economical
reference-dependent choice would require evaluating the actual hard images
and preserving their parent correlations; that optimization is not assumed.

For $a\to\gamma\gamma$, both photon virtualities vanish but $h=m_a^2/2$.
The soft prescription applies only when that invariant is soft too. A heavy
ALP does not become a soft process because its photons are real. Vector
response functions are kept explicit in every channel.

# Interpretation and next step

[The joint source/channel map](odd-power-counted-channel-maps.md) executes
both prescriptions and keeps the first-correction dimensions separate.
Their difference is a truncation comparison, not a numerical error band.
Numerical adequacy still requires coefficient calibration and a specified
kinematic region. No cutoff or value of $\Lambda$ is chosen by this algebra.
