---
title: Odd-sector tree zero and conditional loop plan for ALP three-pseudoscalar decays
status: derived
created: 2026-10-13
updated: 2026-10-13
bibtex_keys: [Alda:2026nrd, Ruiz-Femenia:2003jdx, Ecker:1988te, Ecker:1989yg]
provenance: analytical-derivation
---

# Odd-sector audit of a to three pseudoscalars

## Existing calculation and inputs

The [existing vector-exchange amplitude](a3p-vector-exchange-amplitude.md)
and its executable/numerical kernels calculate the even-parity $G_V$ pole
and mandatory representation-induced contact. Importing generator utilities
from `odd_parity_rcht.py` does not make that interaction odd parity.
No prior $c_i,d_i$ odd contribution to the physical $a\to3P$ amplitude was found.
The existing even-vector calculation retains its older ideal-mixing and
line-shape assumptions; it has not been upgraded to the stage-4 calibrated
vector sector by this audit. That separate update does not introduce an
odd-coefficient dependence into the $G_V$ interaction.

The following result assumes a Lorentz-invariant vacuum, spin-zero external
ALP/pseudoscalars and one epsilon tensor from a single odd insertion.
It needs no CP-conservation assumption or mass-shell relation.

## Executable tree zero

Run [the script](../../../scripts/audit_alp_three_p_odd_parity.py):

```text
.venv/bin/python scripts/audit_alp_three_p_odd_parity.py
```

SymPy builds the determinant of the four external momenta and substitutes
$p_a=p_1+p_2+p_3$, yielding

$$
\epsilon_{\mu\nu\rho\sigma}
p_a^\mu p_1^\nu p_2^\rho p_3^\sigma=0.
$$

It also verifies zero when a fourth routed momentum is any linear combination
of the three independent momenta. Any vacuum four-scalar amplitude with one
epsilon reduces to such a contraction. Hence the single-odd contribution and
its linear sensitivity to every $c_i,d_i$, including the blind directions,
vanish.

In the retained leading tree action, bilinear mixings are absorbed in fields
and propagators and tadpoles are absent. For a connected tree,
$E=2+\sum_v(n_v-2)$. Four external scalar legs permit a single quartic vertex
or two cubic vertices. Two odd cubic VVP/VJP vertices leave vector/photon
external legs, not four scalars; they cannot evade the zero in this tree
truncation. This is a topology argument for the declared action, not a claim
about a general enlarged operator basis or all orders.

## Conditional higher-order calculation plan

No further nonzero odd tree amplitude should be invented for this channel.
If the target accuracy requires double-odd effects:

1. Specify chiral/resonance and $1/N_c$ orders and enumerate all diagrams with
   two odd insertions plus the needed even vertices. Track possible local
   counterterms at that order; the one-pseudoscalar VJP/VVP basis is not a
   complete multi-pseudoscalar basis.
2. Use independent Alda arXiv ALP source spurions and physical generator
   conventions; do not replace the existing even-parity contribution.
3. Derive full tensor propagator contractions, loop integrals and renormalization
   in executable scripts, retaining representation-induced local terms and
   avoiding resonance-saturated ChPT double counting.
4. Match or explicitly parameterize the required local terms. Audit the
   resulting amplitude along the current blind columns; two-epsilon products
   can give nonzero scalar functions, so this tree result cannot settle them.
5. Only then assemble the correction coherently with the existing even amplitude
   and integrate three-body phase space with the proper flavor multiplicities.

Deliver a script/record pair and an error/power-counting estimate before treating
such terms as a phenomenological uncertainty. These loop calculations have
not been performed in this audit.
