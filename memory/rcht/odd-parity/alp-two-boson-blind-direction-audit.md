---
title: ALP two-boson tensor projections onto the three SM blind directions
status: derived
created: 2026-10-13
updated: 2026-10-15
bibtex_keys: [Alda:2026nrd, Ruiz-Femenia:2003jdx, Ecker:1988te, Ecker:1989yg]
provenance: analytical-derivation
---

Canonical tensor and odd-coupling conventions follow [the Levi-Civita audit](levi-civita-convention.md). The linked original scripts and saved artifacts use the explicitly identified reference backend; apply that record's signed coefficient and tensor adapter when interpreting them in the canonical convention.

# Two-boson tensor blindness

## Inputs and assumptions

Use the [physical-ALP source insertions](two-source-alp-rcht-odd-parity-extension.md),
[direct VJP insertions](odd-parity-coupling-prescriptions-and-direct-alp-vgamma.md),
and actual Fourier/index contractions from the
[SM tensor derivation](radiative-tff-amplitude.md).
The VJP/VVP basis comes from the indexed Ruiz-Femenia extracts.
Tree level, single-trace continuation, electromagnetic-compatible ALP sources,
$f_-=0$ for the one-ALP vacuum vertex and canonical spin-one external tensors
are the assumptions. Natural units, metric $(+---)$, $\epsilon^{0123}=+1, \epsilon_{0123}=-1$.

In Alda's invariant conventions let
$X=\mathbb M_a^{-1}(\mathcal C)$ and $Y=c_q^{AP}+X$. Then

$$
u_\mu|_a=-\frac{2\partial_\mu a}{f_a}X,\qquad
\chi_-|_a=-\frac{4ia}{f_a}\{Y,B_0m_q\}.
$$

The derivative and mass weights remain independent. The script uses formal
stripped symbols $D,M$ rather than pretending to normalize a full amplitude.
Neither an equal-spurion relation nor an OPE substitution is used.

## Executable derivation and result

Run [the script](../../../scripts/audit_alp_vector_blindness.py):

```text
.venv/bin/python scripts/audit_alp_vector_blindness.py
```

It calls the explicit indexed Lorentz contractions for both labeled boson
assignments, builds the coefficient polynomials with independent $D,M$, and
differentiates them with respect to the eleven odd coefficients. With $N$
the columns $e_{c_7}$, $e_{d_4}$ and the common $c_2,c_5,c_6$ shift, it obtains

$$
\begin{pmatrix}r_{VJP}\\r_{VVP}\end{pmatrix}N=0_{2\times3}.
$$

The script also repeats the raw contractions with arbitrary polarization
components before transversality or scalar-reference division. Each assignment
gives zero for $c_7$ and the common shift; the two VVP assignments give zero
for $d_4$.

For a wedge tensor $V^{\mu\nu}\propto q^\mu\varepsilon^\nu-q^\nu\varepsilon^\mu$,
the derivative in the seventh/fourth operators produces repeated $q$ factors
inside the epsilon contraction. The common-shift zero is obtained from the
actual sum of the second, fifth and sixth operator contractions.

Canonical external-vector LSZ tensors and photon-induced tensors have this
wedge form. Consequently these projected leading $aVV$ and $aV\gamma$
resonance vertices are insensitive to the three current structural blind
directions. No independent ALP flavor weight changes that algebraic result.

## Limits

This does not supply the missing complete ALP amplitudes, their normalization,
anomaly/counterterm assembly or daughter-channel line shapes. The zero is not
asserted for arbitrary off-shell tensor legs or new local odd operators.
Higher-order diagrams with more than one odd insertion require a separate
audit. Pair this result with the completion plans rather than inferring a
numerical prediction from stripped operator polynomials.
