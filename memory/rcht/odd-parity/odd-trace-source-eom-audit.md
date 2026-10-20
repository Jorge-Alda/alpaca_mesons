---
title: Neutral U(3) source EOM and induced odd topological contacts
status: derived
created: 2026-10-14
updated: 2026-10-14
bibtex_keys: [Alda:2026nrd, Ruiz-Femenia:2003jdx, Roig:2013baa]
provenance: analytical-derivation
---

# Source-preserving EOM audit

## Inputs and assumptions

Use [Alda:2026nrd's leading action](../../extracts/Alda:2026nrd/lo-topological-action.md),
[Appendix-A sources](../../extracts/Alda:2026nrd/appendix-a-expansion.md),
[canonical field conversion](alda-neutral-amplitude-inputs.md), and
[the stage-A trace catalogue](odd-breaking-flavor-enumeration.md).
This calculation retains neutral linear fields, independent diagonal axial
and pseudoscalar sources and an independent topological phase source.
It uses the LO meson action with $F=f_\pi$ and $N_f=3$, natural units and
no fitted physical-mass substitution into its EOM.

This is an audit of a source/EOM subsector, not a complete covariant operator
quotient, a physical-state ALP calculation or an amplitude fit.

## Executable derivation

Run [the script](../../../scripts/audit_odd_trace_sources.py):

```text
.venv/bin/python scripts/audit_odd_trace_sources.py
```

The [transcript](../../../data/vector_sector/odd_trace_source_audit.md)
prints the executed expressions; the
[results JSON](../../../data/vector_sector/odd_trace_source_audit.json)
pins the source/script hashes.

1. Construct the actual canonical neutral generators
   $T_3=\operatorname{diag}(1,-1,0)/\sqrt2$,
   $T_8=\operatorname{diag}(1,1,-2)/\sqrt6$ and $T_0=I/\sqrt3$.
   Write $\Phi_R=\phi_3T_3+\phi_8T_8+\phi_0T_0$.
2. Expand $u=\exp(i\Phi_R/(\sqrt2F))$ through quadratic total field/source
   order, with $\chi=\chi_S+i\chi_P$ and independent $\chi_P$. Differentiate
   before setting any source to zero. SymPy obtains

$$
Y=i\chi_-\big|_{\rm linear}
=\frac{\sqrt2}{F}\{\Phi_R,\chi_S\}-2\chi_P.
$$

3. Introduce the independent dimensionless phase $\theta_{TC}$, which
   specializes to $c_Ga/f_a$ only after the source calculation. Alda:2026nrd's
   invariant becomes

$$
X_{TC}^{(1)}=-\frac{\sqrt2}{F}\operatorname{tr}(\Phi_R)-2\theta_{TC}.
$$

4. Extract the quadratic term of $F^2\operatorname{tr}(\chi_+)/4$.
   Add $-F^2m_0^2X_{TC}^2/12$ and the integrated axial-source cross term
   $\sqrt2F\operatorname{tr}(\Phi_R\partial\cdot a_{\rm ext})$ from
   $u_\mu=2a_{{\rm ext},\mu}-\sqrt2\partial_\mu\Phi_R/F$.
   Differentiate the action with respect to $\phi_3,\phi_8,\phi_0$ to
   obtain their Euler-Lagrange equations. The kinetic term is canonical.
5. Substitute the resulting $\Box\Phi_R$ into $\partial\cdot u$.
   SymPy gives the exact neutral-linear identity

$$
\partial\cdot u=\frac12Y-\frac{m_0^2}{3}X_{TC}^{(1)}I.
$$

The independent axial-source divergences cancel between the definition
of $u_\mu$ and the field EOM. The $\chi_P$ and topological source terms
remain. The script verifies the identity's full matrix residual is zero.

## Induced contacts and partial relations

Define the four VJP trace functions used in
[the normalized kernel record](odd-trace-scalar-kernels.md) as
$q_\ell(A,B;Z)=2F_\ell(A,B,Z)$, $\ell=0,1,2,3$.
An EOM replacement in a divergence operator $D_\ell$ gives

$$
D_\ell=\frac12Y_\ell-\frac{m_0^2}{3}X_{TC}^{(1)}q_\ell(A,B;I),
$$

with the common epsilon/tensor factors understood. Its contact traces are

$$
(q_0,q_1,q_2,q_3)(A,B;I)
=(2\operatorname{tr}AB,\;2\operatorname{tr}A\operatorname{tr}B,
\;2\operatorname{tr}A\operatorname{tr}B,\;6\operatorname{tr}AB).
$$

They span two directions. The script derives two neutral-linear VJP
relations whose induced topological terms cancel:

$$
D_2-D_1=\frac12(Y_2-Y_1),\qquad
D_3-3D_0=\frac12(Y_3-3Y_0).
$$

For the three symmetric VVP trace families, the corresponding contacts
also span two directions, and the script gives one relation,
$D_2-3D_0=(Y_2-3Y_0)/2$ in that family's ordering.

These are **three verified relations in this neutral-linear
divergence/mass subsector**. They are not a complete reduction of all 205
entries, and this audit does not subtract three from that catalogue and
report the result as a full-action independent count.

The remaining EOM replacements induce $X_{TC}\epsilon\operatorname{tr}(AB)$
and $X_{TC}\epsilon\operatorname{tr}(A)\operatorname{tr}(B)$ contacts.
For a photon source the latter may vanish because $\operatorname{tr}Q=0$;
that physical projection is not a general source identity. The contacts
are determined by the coefficient being redefined. Keeping them introduces
no new free parameter by itself; dropping them changes the theory.

## Consequence for the requested operator count

The exact EOM-minimal **new action-operator count is not established yet**.
It belongs after the full one-spurion covariant/source closure and all
permitted EOM reductions, before data calibration. The 205 catalogue entries
and the 29 additional leading trace tensors are pre-EOM source-resolved
counts. The 15 additional scalar EM kernel directions are a separate
projection count. None is silently promoted to the final action count.

To complete that count: retain the fixed induced topological contacts,
audit chiral Bianchi/covariant-derivative curvature remainders and resonance
EOM contacts across the one-spurion blocks, then row-reduce the closed action
space against the original action. Independent new $X_{TC}$ coefficient
functions remain excluded by the selected model; this exclusion does not
permit discarding contacts induced by a field redefinition.

The normalized SM kernels retain their pre-EOM coefficients and therefore
do not implement an invalid octet-only replacement. Applying a tree EOM
to Alda:2026nrd's phenomenological physical-mass scheme would also need its
truncation remainders, especially for eta.
