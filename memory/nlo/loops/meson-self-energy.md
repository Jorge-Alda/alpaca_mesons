---
title: One-loop meson self-energy matrix of the U(3) nonet
status: mixed
created: 2026-11-06
updated: 2026-11-10
bibtex_keys: [Gasser:1984gg, Gasser:1983yg, Kaiser:1998ds, Alda:2026nrd]
provenance: analytical-derivation
---

# Meson self-energy $\Sigma_{ij}(s)$

## Assumptions and inputs

- Convention: $Q(s)=sK-M^2+\Sigma(s)$, as in [loop-corrections-to-alp-meson-lsz-mixing.md](../loop-corrections-to-alp-meson-lsz-mixing.md). Units GeV.
- Vertex: the $\Phi^4$ term of [vertices-and-topologies.md](vertices-and-topologies.md). Propagator contractions: [loop-functions-and-subtraction.md](loop-functions-and-subtraction.md). LO masses and mixing: [neutral-meson-matrices.md](../neutral-meson-matrices.md), local terms: [quadratic-lagrangian.md](../quadratic-lagrangian.md).
- **Assumptions:** no resonance loops ([loop-classes.md](../counting/loop-classes.md)); $\Phi=\sum_a\lambda_aP_a$, $\operatorname{tr}(\lambda_a\lambda_b)=2\delta_{ab}$, $\chi_0=2B_0m_q$.

## Result

Contracting two fields of the $\Phi^4$ vertex gives the loop Lagrangian

$$
\mathcal L_{\rm loop}=\frac1{f_\pi^2}\sum_{c,d}\Big\{\frac{\langle\partial\phi_c\partial\phi_d\rangle}{24}\big[\langle\Phi\lambda_c\Phi\lambda_d\rangle-\langle\Phi^2\lambda_c\lambda_d\rangle\big]
+\frac{\langle\phi_c\phi_d\rangle}{24}\big[\langle\lambda_c\partial_\mu\Phi\lambda_d\partial^\mu\Phi\rangle-\langle\lambda_c\lambda_d\partial_\mu\Phi\partial^\mu\Phi\rangle\big]
+\frac{\langle\phi_c\phi_d\rangle}{48}\sum_{\rm 6\ orderings}\langle\lambda_c\lambda_d\Phi^2\chi_0\rangle\Big\},
$$

where the last sum runs over the six orderings of $\lambda_c,\lambda_d,\Phi,\Phi$ in the trace that are distinct cyclically (explicit list in the script). $\Sigma_{ij}$ is the Hessian: the $\partial P\partial P$ part gives $s\,\partial_s\Sigma$ and the $PP$ part gives $\Sigma(0)$, with $\Sigma_{ij}(s)=\Sigma_{ij}(0)+s\,\partial_s\Sigma_{ij}$.

Pole parts (coefficient of $\lambda(\mu)$), for example:

$$
\partial_s\Sigma_{\pi^\pm}\big|_\lambda=-\frac{B_0(5m_u+5m_d+2m_s)}{3f_\pi^2},\qquad
\partial_s\Sigma_{88}\big|_\lambda=-\frac{B_0(m_u+m_d+2m_s)}{f_\pi^2},
$$

$$
\Sigma_{\pi^\pm}(0)\big|_\lambda=\frac{B_0^2[5(m_u+m_d)^2+2m_s(m_u+m_d)]}{3f_\pi^2}+\frac{2m_0^2B_0(m_u+m_d)}{3f_\pi^2}.
$$

The singlet row and column of $\partial_s\Sigma$ vanish ($\partial_s\Sigma_{00}=\partial_s\Sigma_{30}=\partial_s\Sigma_{80}=0$), because the singlet commutes with every flavour matrix. The full entries in the $(3,8,0)$ basis, in the LO mass eigenbasis, and the numbers at $\mu=0.5,0.77,1.0\,$GeV are printed by `--meson`; they are not repeated here. The test point (illustrative, not fitted) is $m_u=2.16$, $m_d=4.67$, $m_s=93.4\,$MeV, $B_0=2.67\,$GeV, $f_\pi=93\,$MeV, $m_0^2=0.7\,$GeV$^2$.

## Checks

```bash
.venv/bin/python scripts/derive_nlo_loop_self_energies.py --meson
```

- $\Sigma$ is symmetric and vanishes in the chiral limit even for $m_0\neq0$.
- Octet theory (traceless $\Phi$), isospin limit: $M_\pi^2$, $M_K^2$, $M_\eta^2$, the pion, kaon and eta decay constants agree with Gasser:1984gg [@Gasser:1984gg]. Degenerate limit: $M^2(1+\tfrac23\mu_P)$ and $f_\pi(1-3\mu_P)$ with $\mu_P=\frac{M^2}{32\pi^2f_\pi^2}\log\frac{M^2}{\mu^2}$.
- No thresholds: see [vertices-and-topologies.md](vertices-and-topologies.md); $\Sigma_{ij}$ is real for all $s$.

## Pole cancellation against the local terms (`--cancel`)

The normalization can be checked directly from the topological terms in the source conventions: Alda:2026nrd has $-f_\pi^2m_0^2X_{TC}^2/(4N_f)$, while Kaiser:1998ds has $-\tau X_{TC}^2/2$. For $N_f=3$, equating coefficients gives $\tau=f_\pi^2m_0^2/6$ and $4\tau/f_\pi^4=2m_0^2/(3f_\pi^2)$. The unit test `test_topological_normalization` derives and checks this relation. The difference with the scheme record is therefore a factor in that record, not a different normalization of the pion decay constant. Sources: [Alda:2026nrd topological term](../../extracts/Alda:2026nrd/lo-topological-action.md) and [Kaiser:1998ds singlet term](../../extracts/Kaiser:1998ds/qcd-scale-dependence-za.md).

Adding the local terms of [quadratic-lagrangian.md](../quadratic-lagrangian.md) with $L_i\to\Gamma_i\lambda$ and $B\to B^r(1+k_Bm_0^2\lambda/f_\pi^2)$, the poles of $Q(s)$ cancel entry by entry only after a divergent redefinition $P\to(1+\lambda A)P$, with $A$ built from $\{\chi_0,\Phi\}$, $\langle\chi_0
angle\Phi$, $\langle\chi_0\Phi
angle\mathbb 1$ and $\langle\Phi
angle\chi_0$ (coefficients $-\tfrac12,-\tfrac13,\tfrac13,\tfrac13$ with the full $U(3)$ table). The redefinition shifts $Q$ off shell but not the pole positions, so $\Sigma_{ij}$ alone does not fix $\Gamma_4,\Gamma_5,\Gamma_7$ or the $\Lambda_2$ pole; the axial-current mixed function of task 3.1 is needed. The same holds for the $SU(3)$ table of Gasser:1984gg.

Combinations fixed by $\Sigma_{ij}$ and compared with the table of [counterterm-basis.md](../counting/counterterm-basis.md):

| Combination | Derived | Table |
|---|---|---|
| $\Gamma_6-\Gamma_4/2$ | $0$ | $0$ |
| $\Gamma_8-\Gamma_5/2$ | $0$ | $0$ |
| $\Gamma_{25}-\Gamma_4/2$ | $-1/16$ | $-1/16$ |
| $\Gamma_{18}+\Gamma_4-2\Gamma_7$ | $-1/8$ | $-1/8$ |

No mismatch for $\Gamma_6$ or $\Gamma_{25}$ is found. After the redefinition, the pole positions of the nine mesons are independent of $\mu\in\{0.5,0.77,1.0\}\,$GeV to numerical precision.

**Resolved normalization:** The former $1/3$ coefficient arose from an incorrect matching of the topological terms. The source-matched coefficient is $2/3$, and the scheme record and executable check now use it. The resulting running is $\mu\,dB^r/d\mu=-m_0^2B^r/(24\pi^2f_\pi^2)$.
