---
title: Combined delta counting for the U(3) NLO Lagrangian and explicit resonances
status: mixed
created: 2026-11-05
updated: 2026-11-05
bibtex_keys: [Kaiser:1998ds, Kaiser:2000gs, Ecker:1988te, Alda:2026nrd]
provenance: analytical-derivation
---

# Combined $\delta$ counting

## Assumptions and inputs

- The fields are those of [field-content.md](field-content.md).
- The local Goldstone operators are $\mathcal L^{(0)}$, $\mathcal L^{(1)}$ and $\mathcal L^{(2)}$ of [lagrangian.md](../lagrangian.md) [@Kaiser:1998ds; @Alda:2026nrd].
- The counting rules are those of Kaiser:2000gs (see [the extract](../../extracts/Kaiser:2000gs/u1a-counting-and-topological-source.md)) [@Kaiser:2000gs]. This record extends them to the resonance couplings of Ecker:1988te.
- **Assumption:** the resonance couplings scale like $F\sim\sqrt{N_c}$ and the masses like $N_c^0$. This is the usual large-$N_c$ scaling of tree-level meson couplings, and it is not derived here.
- **Assumption:** one power of $1/N_c$ suppresses OZI-violating differences between singlet and octet resonance couplings.
- The ALP is kept at first order in $1/f_a$. $1/f_a$ is not part of $\delta$, and $m_a^2$ counts like $p^2$.

The executable derivation is [`scripts/derive_nlo_u3_resonance_counting.py`](../../../scripts/derive_nlo_u3_resonance_counting.py):

```bash
.venv/bin/python scripts/derive_nlo_u3_resonance_counting.py --orders
```

## Step 1: the expansion parameter

$$
p=\mathcal O(\sqrt\delta),\qquad m_q=\mathcal O(\delta),\qquad 1/N_c=\mathcal O(\delta).
$$

Each derivative, $\bm{u_\mu}$ and $f_\pm^{\mu\nu}$ counts as $p$, and each $\bm{\chi_\pm}$ counts as $p^2$. $X_{TC}$, including the ALP through $c_G a/f_a$, counts as $p^0$.

## Step 2: the $N_c$ power of a Goldstone operator

A term of chiral dimension $p^{2n}$ with coefficient $\propto N_c^k$ is of order $\delta^{\,n-k}$. The leading Lagrangian is $\delta^0$. For an operator built from Goldstone fields and sources, the $N_c$ power follows from the large-$N_c$ expansion of the effective Lagrangian of Kaiser:2000gs:

- The pure-gluon vacuum energy is $\propto N_c^2$.
- A single flavor trace is $\propto N_c$.
- Each further flavor trace costs $1/N_c$.
- Each factor of $X_{TC}=-(\psi+\theta)$ costs $1/N_c$, because $\theta$ enters only as $\theta/N_c$ and $\psi=-i\langle\log\bm U\rangle$ is itself a flavor trace.

Hence

$$
k = 2-(\text{number of flavor traces})-(\text{number of factors of }X_{TC}).
$$

The script checks this rule against the explicit prefactor of each term, using $F^2\sim N_c$, $\Lambda_{1,2}\sim1/N_c$ and $\tau\sim N_c^0$.

## Step 3: the Goldstone operators

| Operator | $k$ | $p^{2n}$ | Order |
|---|---|---|---|
| $\frac{f_\pi^2}{4}\langle\bm{u_\mu}\bm{u^\mu}\rangle$, $\frac{f_\pi^2}{4}\langle\bm{\chi_+}\rangle$ | 1 | $p^2$ | $\delta^0$ |
| $-\frac{f_\pi^2m_0^2}{4N_f}X_{TC}^2$ | 0 | $p^0$ | $\delta^0$ |
| $L_5\langle\bm{u_\mu}\bm{u^\mu}\bm{\chi_+}\rangle$, $\frac{L_8}{2}\langle\bm{\chi_+}^2+\bm{\chi_-}^2\rangle$ | 1 | $p^4$ | $\delta$ |
| $\frac{f_\pi^2}{12}\Lambda_1\langle\bm{u_\mu}\rangle\langle\bm{u^\mu}\rangle$ | 0 | $p^2$ | $\delta$ |
| $i\frac{f_\pi^2}{12}\Lambda_2X_{TC}\langle\bm{\chi_-}\rangle$ | 0 | $p^2$ | $\delta$ |
| $L_4\langle\bm{u_\mu}\bm{u^\mu}\rangle\langle\bm{\chi_+}\rangle$, $L_6\langle\bm{\chi_+}\rangle^2$, $L_7\langle\bm{\chi_-}\rangle^2$, $L_{18}\langle\bm{u_\mu}\rangle\langle\bm{u^\mu}\bm{\chi_+}\rangle$ | 0 | $p^4$ | $\delta^2$ |
| $-iL_{25}X_{TC}\langle\bm{\chi_+}\bm{\chi_-}\rangle$ | 0 | $p^4$ | $\delta^2$ |
| single-trace $\mathcal O(N_cp^6)$ | 1 | $p^6$ | $\delta^2$ |
| one-loop graph with LO vertices | 0 | $p^4$ | $\delta^2$ |

This reproduces the split of [lagrangian.md](../lagrangian.md) into $\mathcal L^{(1)}$ and $\mathcal L^{(2)}$. From the topological term, $m_0^2=4N_f\tau/f_\pi^2=\mathcal O(\delta)$, so the $\eta'$ mass counts like a Goldstone mass. $L_7$ is $\mathcal O(N_c^0)$ in the $U(3)$ theory. Its $\mathcal O(N_c^2)$ size in $SU(3)$ ChPT comes from integrating out $\eta'$, which is a dynamical field here. The last row shows that one-loop Goldstone graphs and $\mathcal L^{(2)}$ are of the same order. This is why the loop divergences are absorbed by the $\delta^2$ constants (see the counterterm records of this folder).

## Step 4: resonance couplings and masses

| Quantity | $k$ | $p^{2n}$ | Order |
|---|---|---|---|
| $M_R^2$ | 0 | $p^0$ | $\delta^0$ |
| $c_d,c_m,d_m,F_V,G_V,F_A$ (vertex with one resonance) | $\tfrac12$ | $p^2$ | $\delta^{1/2}$ |
| singlet $\widetilde c_d,\widetilde c_m,\widetilde d_m$ | $\tfrac12$ | $p^2$ | $\delta^{1/2}$ |
| OZI-violating difference, e.g. $\widetilde c-c/\sqrt3$ | $-\tfrac12$ | $p^2$ | $\delta^{3/2}$ |

A single vertex has a fractional order. Only graphs as a whole are assigned to the declared order.

## Step 5: single-resonance exchange

Tree-level exchange between two Goldstone/source vertices gives

$$
\frac{(\sqrt{N_c}\,p^2)^2}{M_R^2-s}
=\frac{N_cp^4}{M_R^2}\left(1+\frac{s}{M_R^2}+\cdots\right).
$$

| Term | $k$ | $p^{2n}$ | Order |
|---|---|---|---|
| $1/M_R^2$ | 1 | $p^4$ | $\delta$ |
| $s/M_R^4$ | 1 | $p^6$ | $\delta^2$ |
| one OZI-violating vertex $\times1/M_R^2$ | 0 | $p^4$ | $\delta^2$ |

The leading term has the order of $\mathcal L^{(1)}$, as resonance saturation of $L_5$ and $L_8$ requires. The $s/M_R^4$ term has the order of the $\mathcal O(N_cp^6)$ terms of $\mathcal L^{(2)}$.

**Unexpanded propagators.** This expansion is used only to assign orders. For $\sqrt s\lesssim1\,$GeV, $s/M_R^2$ is not small: $s/M_\rho^2\approx1.7$ at $1\,$GeV. Every resonance-exchange graph is therefore kept with its full propagator $1/(M_R^2-s)$. A graph is included at the declared order when its low-energy limit is of order $\delta$ or $\delta^2$. This is a **choice of the framework**. It resums a series of local $N_cp^{2n}$ terms that is not convergent near $1\,$GeV.

## Step 6: Goldstone–pseudoscalar-resonance mixing

The coupling $id_m\langle P\bm{\chi_-}\rangle$ contains a bilinear $R^TBP$ with $B\sim d_m\chi_0/F$ (see [resonance-saturation-and-dynamical-masses.md](../resonance-saturation-and-dynamical-masses.md)). For canonically normalized fields, $B\sim N_c^0p^2$ and the induced mass shift is $B^2/M_R^2\sim N_c^0p^4$. Inserting canonical-field coefficients into $\delta^{\,n-k}$ gives an order that is off by one power of $\delta$, because the canonical LO mass term itself reads $N_c^0p^2$. The order must be read from the coefficient normalized with $F^2$. That coefficient is $F^2B^2/M_R^2\sim N_cp^4$, which is $\delta$, the order of the $L_8$-type term it saturates. The ratio of the shift to the LO mass, $\delta$, is independent of normalization.

## Result

The declared calculation is $\mathcal O(\delta^2)$ in the Goldstone/ALP sector. It contains:

- $\mathcal L^{(0)}$, $\mathcal L^{(1)}$ and $\mathcal L^{(2)}$ of [lagrangian.md](../lagrangian.md);
- tree-level exchange of the resonance nonets of [field-content.md](field-content.md), with unexpanded propagators, including Goldstone–resonance bilinear mixing;
- one-loop graphs at $\mathcal O(N_c^0p^4)$, which are classified in the loop-class record of this folder.

## Verification

The script asserts every order in the tables above and exits with code 0. It also runs its unit tests with `--test`. All orders in this record are copied from its `--orders` output.
