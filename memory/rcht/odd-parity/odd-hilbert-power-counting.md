---
title: Chiral grading and correction filtration for the bounded odd action
status: mixed
created: 2026-10-14
updated: 2026-10-14
bibtex_keys: [Alda:2026nrd, Graf:2020yxt, Kaiser:2000gs]
provenance: analytical-derivation
---

# Inputs and assumptions

Use [the bounded Hilbert quotient](odd-bounded-hilbert-series.md),
[the prior proposed hierarchy](odd-action-power-counting-proposal.md),
[covariant closure](odd-covariant-closure.md), and
[Kaiser--Leutwyler counting](../../extracts/Kaiser:2000gs/u1a-counting-and-topological-source.md).
The selected vector representation remains unconstrained and canonical.
Assigning it formal chiral order zero is a bookkeeping convention for this
explicit-resonance action, not an expansion theorem for on-shell vectors.
No physical mass substitutions precede source differentiation.

# Homogeneous grade

Use

$$
w(\nabla)=w(u_\mu)=1,\quad
w(f_\pm)=w(X)=w(Y)=w(Z_{TC})=2,\quad w(V_{\mu\nu})=0.
$$

The light-field weights agree with the source-character setup
[@Graf:2020yxt]. The assignment $w(Z_{TC})=2$ is consistent with retaining
the singlet mass at leading U(3) order. The covariant identities preserve
total weight, e.g.

$$
\nabla u\sim f_-,\quad [\nabla,\nabla]\sim [u,u]+f_+,
\quad Y=2\nabla\cdot u+2Z_{TC}I.
$$

The leading VJP seeds therefore have order $p^4$, leading VVP $p^2$;
their higher blocks have order $p^6$ and $p^4$, respectively. Equal mass
dimension across sectors is neither required nor a measure of relative
size. Higher VJP/VVP engineering dimensions are seven/six, as in the
normalized basis; $1/M_V^3$ and $1/M_V^2$ normalize dimensions without
establishing a breakdown scale.

Separate derivative and quark-spurion numbers are not exact quotient
grades: the source-preserving EOM relates $Y$ to a derivative divergence
and a topological contact. Separate trace number is not exact either:
two of the thirteen saved relations mix single/double traces. Accordingly
one should not attach independent mass/OZI fugacities to the quotient
representatives and assume that each relation remains homogeneous.

# Parent filtration and executed dimension

The equal-order correction assumption of the prior proposal is

$$
\epsilon_p=p_{\rm soft}^2/\Lambda^2\sim
\epsilon_\chi=B_0m_q/\Lambda^2\sim
\epsilon_N=1/N_c\sim\delta.
$$

Define the parent cost $c=\ell+n_{\rm tr}-1$, where $\ell=0$ for the five
leading seed families and $\ell=1$ for the five higher families. In particular,
$Vf_-X$ is leading and $Vf_-X\nabla\nabla$ is its derivative correction.
This treats leading $Y$ sources as leading, rather than assigning them an
additional relative mass penalty. Additional traces and the higher seed
order each cost one; their combination costs two. This is a nominal tree
model filtration, with all fixed source/curvature/EOM remainders retained.

Let $S_k\subset F$ be the span of parents with $c\le k$. The actual retained
space is its image in $Q=F/R$, not a subset of arbitrarily selected full-basis
labels. The script computes

$$
\dim\operatorname{im}(S_k\to Q)
=|S_k|-\operatorname{rank}R
+\operatorname{rank}R_{\overline S_k},
$$

where rows of $R$ are independent relations and $R_{\overline S_k}$ is
their restriction to omitted columns. This follows by computing
$\dim(R\cap S_k)=\operatorname{rank}R-\operatorname{rank}R_{\overline S_k}$.

| Ceiling | Selected raw parents | Independent image dimension |
|---|---:|---:|
| Leading single traces | 12 | 12 |
| Through first corrections | 107 | 101 |
| Through mixed corrections | 241 | 228 |

The associated graded **filtration** polynomial is

$$
H_{\mathrm{gr}}(\delta)=12+89\delta+127\delta^2.
$$

This explicitly resolves the previously uncomputed parent-action truncation
count. It does not establish an observable/neutral/ALP rank or change the
physical original $c_i,d_i$ coefficient definitions. The leading twelve
contain an extra covariant direction beyond the original eleven.

# Domain of validity

Kaiser--Leutwyler establish the simultaneous expansion for U(3) Goldstones
[@Kaiser:2000gs]; applying relative OZI counting to this explicit-resonance
seed model remains an assumption. Formal chiral weight one on a derivative
does not make a derivative acting on a vector of momentum $O(M_V)$ soft.
Near resonance poles or for $a\to VV$, an additional kinematic prescription
is required before numerical suppression can be claimed. Retain propagators
explicitly. A first-correction tree truncation with loops deferred is not a
complete next order in $1/N_c$.

The Hilbert series counts permitted coefficients; it cannot derive their
natural size, a value of $\Lambda$, priors, or convergence. Adequacy requires
normalized kernels, nested predictions and comparison with data in the stated
domain. If a suppressed parent gives a channel's first nonzero amplitude,
its contribution must be assessed at that channel's leading accuracy.

The next executable step has now been carried out in
[the resonance derivative prescription](odd-resonance-derivative-counting.md)
and [the joint channel map](odd-power-counted-channel-maps.md). Their
EOM-compatible soft primitive ranks are $(7,29,70)$; hard-derivative promotion
gives $(12,40,70)$. Those projections do not replace the parent-action counts
above or establish numerical convergence.
