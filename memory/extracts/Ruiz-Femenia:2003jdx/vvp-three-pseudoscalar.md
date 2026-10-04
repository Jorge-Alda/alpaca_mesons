---
title: Ruiz-Femenia:2003jdx VVP operators and three-pseudoscalar amplitude
status: literature
created: 2026-10-04
updated: 2026-10-04
bibtex_keys: [Ruiz-Femenia:2003jdx]
provenance: literature
---

# VVP Operators and Three-Pseudoscalar Amplitude in Ruiz-Femenia:2003jdx

## VVP Operator Basis (Eq. (11), source lines 384--397)

The paper defines four independent VVP operators in the antisymmetric-tensor formulation [@Ruiz-Femenia:2003jdx, Eq. (11)]:

$$
\begin{aligned}
\mathcal O^1_{VVP} &= \epsilon_{\mu\nu\rho\sigma}\,
\langle \{V^{\mu\nu},V^{\rho\alpha}\} \nabla_{\alpha}u^{\sigma}\rangle, \\
\mathcal O^2_{VVP} &= i\,\epsilon_{\mu\nu\rho\sigma}\,
\langle \{V^{\mu\nu},V^{\rho\sigma}\}\, \chi_{-}\rangle, \\
\mathcal O^3_{VVP} &= \epsilon_{\mu\nu\rho\sigma}\,
\langle \{\nabla_{\alpha}V^{\mu\nu},V^{\rho\alpha}\} u^{\sigma}\rangle, \\
\mathcal O^4_{VVP} &= \epsilon_{\mu\nu\rho\sigma}\,
\langle \{\nabla^{\sigma}V^{\mu\nu},V^{\rho\alpha}\} u_{\alpha}\rangle.
\end{aligned}
$$

The odd-parity VVP Lagrangian is [@Ruiz-Femenia:2003jdx, Eq. (12), lines 424--426]:

$$
\mathcal L_{VVP} = \sum_{a=1}^{4} d_a\, \mathcal O^a_{VVP}.
$$

**Crucial caveat (source lines 431--433):** "We stress that the set defined above is a complete basis for constructing vertices with only **one-pseudoscalar**; for a larger number of pseudoscalars additional operators may emerge."

## Which $d_i$ Enter a Genuine Three-Pseudoscalar Amplitude?

### $\omega\to 3\pi$ Amplitude (Eq. (30), source lines 913--924)

The $\rho$-mediated $\omega\to\pi^+\pi^-\pi^0$ amplitude (including cyclic permutations) is:

$$
i\mathcal M_{\omega\to3\pi}=i\epsilon_{\alpha\beta\rho\sigma}
k_1^\alpha k_2^\beta k_3^\rho\epsilon_\omega^\sigma
\frac{8G_V}{M_\omega F^3}\left[
\frac{m_\pi^2(d_1+8d_2-d_3)+(M_\omega^2+s_{12})d_3}{M_V^2-s_{12}}
+\{s_{12}\to s_{13}\}+\{s_{12}\to s_{23}\}\right].
$$

**The $\rho$-exchange amplitude depends on the combinations:**
- $d_1 + 8d_2 - d_3$ (multiplying $m_\pi^2$)
- $d_3$ (multiplying $M_\omega^2 + s_{ij}$)

**$d_4$ does not appear** in this $\rho$-exchange amplitude.

### Direct $\omega\to 3\pi$ Contact Term

The paper **explicitly states** that a direct $\omega\to 3\pi$ term is required phenomenologically (source lines 955--956, 974--985): "the intermediate meson exchange does not account entirely for the $\omega$ decay into three pions, and the direct terms must be considered... In our effective theory, these terms would be obtained by writing down the operators which give rise to local contributions to $\omega\to 3\pi$."

**However, the paper does not write down the explicit direct-contact operators for $3\pi$.** The VVP basis (Eq. 11) is stated to be complete only for **one-pseudoscalar** vertices. For three pseudoscalars, "additional operators may emerge" (line 433). The direct $\omega 3\pi$ vertex is **not provided** in the published basis.

### Short-Distance Constraints on $d_i$ (Eq. (19), source lines 683--686)

Matching the $\langle VVP\rangle$ Green's function to the QCD OPE yields:

$$
\begin{aligned}
d_1 + 8d_2 &= -\frac{N_C}{64\pi^2}\frac{M_V^2}{F_V^2} + \frac{F^2}{4F_V^2}, \\
d_3 &= -\frac{N_C}{64\pi^2}\frac{M_V^2}{F_V^2} + \frac{F^2}{8F_V^2}.
\end{aligned}
$$

**$d_4$ is not constrained** by the $\langle VVP\rangle$ short-distance matching. The constraints fix only the combinations $d_1+8d_2$ and $d_3$ individually.

### $\pi\to\gamma\gamma$ Correction (Eq. (34), source lines 1036--1043)

The VVP contribution to $\pi\to\gamma\gamma$ at $\mathcal O(m_\pi^2)$ is proportional to $(d_1+8d_2-d_3)$ [@Ruiz-Femenia:2003jdx, Eq. (34), line 1038]. After imposing short-distance constraints, this combination becomes fixed and yields a tiny $\sim 1\%$ correction.

## Summary: $d_i$ in Three-Pseudoscalar Amplitudes

| Coupling | Enters $\rho$-exchange $\omega\to 3\pi$? | Enters direct $V\to 3\pi$? | Constrained by $\langle VVP\rangle$ OPE? |
|----------|------------------------------------------|----------------------------|------------------------------------------|
| $d_1$    | Only in combination $d_1+8d_2-d_3$       | **Unknown** (basis incomplete for $3P$) | Only in combination $d_1+8d_2$ |
| $d_2$    | Only in combination $d_1+8d_2-d_3$       | **Unknown** (basis incomplete for $3P$) | Only in combination $d_1+8d_2$ |
| $d_3$    | Yes (two distinct kinematic structures)  | **Unknown** (basis incomplete for $3P$) | Yes, individually |
| $d_4$    | **No**                                   | **Unknown** (basis incomplete for $3P$) | **No** |

**Conclusion:** For the $\rho$-exchange mechanism (which is what the BaBar ISR fit assumes), the three-pion amplitude depends on **$d_1+8d_2-d_3$ and $d_3$**. The direct $V\to 3\pi$ contact term—which the paper argues is necessary and leading in $1/N_C$—**is not provided in the published operator basis**, which is explicitly complete only for one-pseudoscalar vertices. $d_4$ does not enter the $\rho$-exchange amplitude and is unconstrained by the $\langle VVP\rangle$ OPE matching.

## $e^+e^-\to 3\pi$ Amplitude Structure

The paper does not give an explicit $e^+e^-\to 3\pi$ amplitude. The $e^+e^-\to V$ vertex comes from $\mathcal L_2(V)$ (photon-vector mixing via $F_V$), and $V\to 3\pi$ proceeds via $\rho$ exchange as above. The full amplitude would be a coherent sum over $V=\omega,\phi,\omega',\omega''$ with propagators $D_V(m)$ and the $\rho$-exchange vertex. The relative phases are not predicted by the short-distance matching.