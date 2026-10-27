---
title: Conditional regular remainder of the retained ALP vector numerator
status: mixed
created: 2026-10-27
updated: 2026-10-27
bibtex_keys: [Alda:2026nrd, Ecker:1988te, Ecker:1989yg, Ruiz-Femenia:2003jdx]
provenance: analytical-derivation
---

# Scope and inputs

This record separates the strict vector pole from the regular polynomial
remainder of the retained ordered ALP production numerator, conditional on
using its existing affine off-shell continuation. The source matrices,
photon response $W_j$, numerators and normalization are those of the
[source-resolved tree TFF](alp-doubly-virtual-tff.md) and
[pole-factorization record](alp-three-pion-vector-pole-factorization.md),
with conventions from Alda:2026nrd and the retained RChT operators
[@Alda:2026nrd; @Ecker:1988te; @Ecker:1989yg; @Ruiz-Femenia:2003jdx].
Natural units are used. The result is not a representation-independent
physical off-shell amplitude.

# Executable division

Run `.venv/bin/python scripts/derive_alp_three_pion_regular_remainder.py`.
For either ordered photon placement, polynomial division gives
$$
\frac{N_i(s_1,s_2)}{M_i^2-s_1}
=\frac{N_i(M_i^2,s_2)}{M_i^2-s_1}
+\frac{N_i(s_1,s_2)-N_i(M_i^2,s_2)}{M_i^2-s_1}.
$$
The script reconstructs the retained VJP and VVP numerators from their
existing coefficient expressions, then performs this division in SymPy. The
regular quotient is independent of $s_1$ for these affine numerators. In the
ordered scalar coefficient convention of the pole-factorization record it
is
$$
\begin{aligned}
N_{i,\mathrm{reg}}(s_2)={}&
-\frac{4e}{M_V}(c_1-c_2-c_5+2c_6)
\operatorname{tr}(T_i Q U_a)\\
&-4\sqrt2\,e\,d_3\sum_j W_j(s_2)
\operatorname{tr}(T_iT_jU_a).
\end{aligned}
$$
The first term is the VJP remainder; the second is the VVP remainder with
the other photon response retained. The script verifies the decomposition by
recombining pole residue and quotient into each original numerator.

# Interpretation and limits

This quotient is the regular term associated with the declared continuation
of the retained numerator. It does not determine general regular production
terms, direct $V3\pi$ contact coefficients, rescattering contributions, or
the complementary local term of the full antisymmetric-tensor propagator.
The [off-shell sequential-exchange record](alp-isoscalar-vector-three-pion-offshell.md)
and [contact-basis record](../../rcht/odd-parity/vector-three-pion-contact.md)
show that independent endpoint terms remain. The script's focused test
passes; no numerical coupling point is selected.
