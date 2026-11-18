---
title: Polar and Cartesian field coordinates in finite-mixing amplitudes
status: derived
created: 2026-11-18
updated: 2026-11-18
bibtex_keys: [Alda:2026nrd, DiLuzio:2026poh, Quevillon:2019zrd]
provenance: analytical-derivation
---

# Polar and Cartesian field coordinates in finite-mixing amplitudes

## Inputs and assumptions

Use natural units and the polar parametrization of [the radial and phase mass-source record](alp-murayama-radial-phase-mass-coefficients.md), following DiLuzio:2026poh [@DiLuzio:2026poh]. The tree-level pole vectors, all-incoming vertex convention and full mixed inverse propagator $D(k^2)$ follow [the pole-projected amplitude record](alp-lo-pole-projected-tree-amplitudes.md) and Alda:2026nrd [@Alda:2026nrd]. Work near a vacuum with $f_\varphi\ne0$, with stable, simple, isolated external poles and nonsingular quadratic kinetic terms. No expansion in the mixing is made. Quevillon:2019zrd compares linear and polar parametrizations in an axion model and demonstrates the equality of complete amplitudes after its additional fermion redefinition [@Quevillon:2019zrd].

## Field coordinates

Write the same dimensionful complex field in polar and Cartesian coordinates:
$$
\Phi=\frac{f_\varphi+\sigma}{\sqrt2}e^{ia/f_\varphi}
=\frac{f_\varphi+s+ip}{\sqrt2}.
$$
Thus
$$
s=(f_\varphi+\sigma)\cos(a/f_\varphi)-f_\varphi,\qquad
p=(f_\varphi+\sigma)\sin(a/f_\varphi),
$$
and locally around the vacuum the inverse begins
$$
\sigma=s+\frac{p^2}{2f_\varphi}+\cdots,\qquad
 a=p-\frac{sp}{f_\varphi}+\cdots.
$$
Collect the polar fields, including mesons and the radial field, in $\Phi^I$, and the corresponding Cartesian fields in $\hat{\Phi}^I$. This multiplet notation is distinct from the complex field $\Phi$ above. The local map is
$$
\Phi^I=\hat{\Phi}^I+\frac12 C^I{}_{JK}\hat{\Phi}^J\hat{\Phi}^K+\frac16 C^I{}_{JKL}\hat{\Phi}^J\hat{\Phi}^K\hat{\Phi}^L+\cdots.
$$
The map has identity Jacobian at the vacuum. Hence the quadratic operator and its normalized tree-level pole vectors agree in these coordinates. The nonlinear terms alter interaction vertices. In particular, the polar kinetic term contains $\frac12(1+\sigma/f_\varphi)^2(\partial a)^2$, while the Cartesian kinetic term is $\frac12[(\partial s)^2+(\partial p)^2]$.

## Four-point amplitude with finite mixing

Let $v_i$ be the pole vector for each of the four all-incoming external legs, including the physical ALP and the crossed final mesons, and impose $D(p_i^2)v_i=0$. Let $\mathcal V_{ij,I}$ denote the cubic vertex projected onto $v_i,v_j$ and $C_{ij}^I=C^I{}_{JK}v_i^Jv_j^K$. For the partition $12|34$, put $q=-(p_1+p_2)=p_3+p_4$. After transforming the fields and projecting external legs,
$$
\mathcal V'_{12,I}=\mathcal V_{12,I}+C_{12}^J D_{JI}(q^2),\qquad
\mathcal V'_{34,I}=\mathcal V_{34,I}+D_{IJ}(q^2)C_{34}^J.
$$
The omitted terms have an inverse propagator on an external pole vector and vanish. The projected four-point vertex becomes
$$
\mathcal V'_4=\mathcal V_4+\sum_{ij|kl}\left[
 C_{ij}^I\mathcal V_{I,kl}+\mathcal V_{ij,I}C_{kl}^I+C_{ij}^I D_{IJ}(q_{ij}^2)C_{kl}^J\right].
$$
The sum is over the three partitions $12|34$, $13|24$, $14|23$; terms containing $C^I{}_{JKL}$ vanish on the external poles. For each partition the transformed exchange contribution is
$$
\begin{aligned}
-\mathcal V'_{ij}D^{-1}\mathcal V'_{kl}={}&
-\mathcal V_{ij}D^{-1}\mathcal V_{kl}
-\mathcal V_{ij}C_{kl}-C_{ij}\mathcal V_{kl}-C_{ij}DC_{kl}.
\end{aligned}
$$
The last three terms cancel the corresponding changes in the four-point contact vertex. Therefore the full projected on-shell tree amplitude is equal in the two parametrizations for every external pole assignment, including $A\to P_1P_2P_3$. The inverse is the full mixed propagator, including radial exchange where present. The equality extends to an exchanged pole by comparison of residues or by taking a limit. Individual contact and exchange terms depend on the coordinates.

## General amplitudes

For a local invertible field change, a formal change of variables with sources coupled to the same interpolating operators gives
$$
Z[J]=\int\mathcal D\Phi\,e^{iS[\Phi]+i\int J_I\Phi^I}
=\int\mathcal D\hat\Phi\,\det\!\left(\frac{\delta\Phi}{\delta\hat\Phi}\right)
 e^{iS[\Phi(\hat\Phi)]+i\int J_I\Phi^I(\hat\Phi)}.
$$
The equality retains the transformed source and measure. For sources coupled linearly to either coordinate, the corresponding Green functions differ by composite-operator insertions. LSZ reduction with the appropriate one-particle residues gives the same on-shell S-matrix under the usual equivalence-theorem assumptions; at loop order the regulated measure, renormalization and residues require consistent treatment. The scalar field change alone is not a fermion chiral rotation. If a separate chiral rotation is made, include its anomalous contribution as well. For unstable states compare an appropriate pole residue or a process with stable external states. The tree-level four-point cancellation above does not by itself prove the quantum equivalence theorem.

## Verification

The displayed cancellation follows algebraically from the transformed vertices and $D^{-1}D=\mathbf1$ away from internal poles. [The SymPy calculation](../../../scripts/derive_polar_cartesian_four_point.py) expands the inverse polar map and derives the matrix exchange cancellation. It takes the contact shift from the chain rule as an input; it does not independently derive that shift or prove the general LSZ statement. It does not provide an independent numerical amplitude comparison.

```sh
.venv/bin/python scripts/derive_polar_cartesian_four_point.py
.venv/bin/python -m unittest scripts.test_derive_polar_cartesian_four_point
```
