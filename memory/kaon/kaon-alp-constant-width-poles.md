---
title: Physical neutral-kaon and ALP poles with constant widths
status: derived
created: 2026-11-16
updated: 2026-11-16
bibtex_keys: [Alda:2026nrd, ParticleDataGroup:2024cfk]
provenance: analytical-derivation
---

# Physical neutral-kaon and ALP poles with constant widths

## Inputs, conventions and assumptions

The unmixed physical-kaon basis follows [the kaon projection record](neutral-kaon-amplitude-projection.md) and [the neutral-kaon mixing extract](../extracts/ParticleDataGroup:2024cfk/neutral-kaon-mixing.md) [@ParticleDataGroup:2024cfk]. The distinction between an exact quadratic mixing calculation and its weak-mixing limit follows [the finite ALP-meson mixing record](../alp/meson-mixing/finite-mixing/spectrum-and-poles/quadratic-mixing/finite-alp-meson-mixing.md) and [the near-degeneracy record](../alp/meson-mixing/finite-mixing/pole-limits/alp-lo-near-degeneracy-two-level.md) [@Alda:2026nrd]. The latter treats a real mass-squared matrix; the present complex-energy Hamiltonian includes the physical weak kaon widths and is an effective extension, not a direct identification of its entries with the chiral mass-squared entries. The exact degenerate strong-sector kaon result in [the $K^0$ doublet record](../alp/meson-mixing/finite-mixing/spectrum-and-poles/pole-equations/alp-lo-k0-doublet-exact-pole.md) does not include the observed weak $K_S$--$K_L$ splitting.

Use $\hbar=c=1$, with all energies, masses, widths and Hamiltonian entries in GeV. Let $\lambda_i=m_i-i\Gamma_i/2$ for $i=S,L,a$ be the constant complex energies before ALP--kaon mixing, and let $\mu_{ia}$ and $\mu_{ai}$ be the directed entries from $a$ to $K_i$ and from $K_i$ to $a$, respectively. These notation choices indicate the direction of mixing with subscripts. The baseline kaon Hamiltonian preserves CPT, and its two eigenvalues are distinct. The three-state Hamiltonian is an assumption for a low-energy, constant-width treatment; it does not calculate the entries from ALP-quark couplings. The observed $K_S,K_L$ are the baseline poles; after ALP mixing the eigenstates are identified by continuous connection where this is unambiguous. Diagonal ALP-induced corrections, if present, are included in the specified diagonal inputs. A physical decay matrix must satisfy the constraints appropriate to the underlying open system, including positivity; independent directed entries below are a parametrization and do not assert that every choice is physically allowed.

## Step 1: physical kaons and directed entries

In the CPT-preserving convention, $|K_S\rangle=p|K^0\rangle+q|\bar K^0\rangle$ and $|K_L\rangle=p|K^0\rangle-q|\bar K^0\rangle$, with $|p|^2+|q|^2=1$ and nonzero $p,q$. When $|p|\ne|q|$, these right eigenstates are not orthogonal. The coordinate transformation has these states as columns, so its inverse supplies the dual left basis. If $H_{0a}$, $H_{\bar0a}$, $H_{a0}$ and $H_{a\bar0}$ are entries of the flavor-basis Hamiltonian, SymPy gives

$$
\begin{aligned}
\mu_{Sa}&=\frac{H_{0a}}{2p}+\frac{H_{\bar0a}}{2q},&
\mu_{La}&=\frac{H_{0a}}{2p}-\frac{H_{\bar0a}}{2q},\\
\mu_{aS}&=pH_{a0}+qH_{a\bar0},&
\mu_{aL}&=pH_{a0}-qH_{a\bar0}.
\end{aligned}
$$

The inverse basis transformation is essential: ordinary conjugate bras are not the dual basis in general. For the exact complex convention

$$
p=\frac{1+\epsilon_K}{\sqrt{2(1+|\epsilon_K|^2)}},\qquad q=\frac{1-\epsilon_K}{\sqrt{2(1+|\epsilon_K|^2)}},
$$

the right eigenvectors are the columns $(p,q,0)^T$, $(p,-q,0)^T$, and $(0,0,1)^T$. The dual left basis consists of the rows of the inverse transformation, not ordinary conjugate bras. The general directed entries, written directly in terms of the four flavor-basis matrix elements, are

$$
\begin{aligned}
\mu_{Sa}&=\frac{\sqrt{2(1+|\epsilon_K|^2)}}{2(1+\epsilon_K)}H_{0a}+\frac{\sqrt{2(1+|\epsilon_K|^2)}}{2(1-\epsilon_K)}H_{\bar 0a},&
\mu_{La}&=\frac{\sqrt{2(1+|\epsilon_K|^2)}}{2(1+\epsilon_K)}H_{0a}-\frac{\sqrt{2(1+|\epsilon_K|^2)}}{2(1-\epsilon_K)}H_{\bar 0a},\\
\mu_{aS}&=\frac{(1+\epsilon_K)H_{a0}+(1-\epsilon_K)H_{a\bar 0}}{\sqrt{2(1+|\epsilon_K|^2)}},&
\mu_{aL}&=\frac{(1+\epsilon_K)H_{a0}-(1-\epsilon_K)H_{a\bar 0}}{\sqrt{2(1+|\epsilon_K|^2)}}.
\end{aligned}
$$

These expressions allow general CP-violating ALP couplings, without relating the four flavor entries. Now assume CP conservation in the ALP--kaon mixing, with $CP|K^0\rangle=|\bar K^0\rangle$ and $CP|a\rangle=-|a\rangle$. CP is unitary, so $H_{\bar 0a}=-H_{0a}$ and $H_{a\bar 0}=-H_{a0}$, without complex conjugation. Keeping the baseline kaon parameter $\epsilon_K$ finite, the entries become

$$
\mu_{Sa}=-\frac{\sqrt{2(1+|\epsilon_K|^2)}\,\epsilon_K H_{0a}}{1-\epsilon_K^2},\qquad
\mu_{La}=\frac{\sqrt{2(1+|\epsilon_K|^2)}\,H_{0a}}{1-\epsilon_K^2},\qquad
\mu_{aS}=\frac{2\epsilon_K H_{a0}}{\sqrt{2(1+|\epsilon_K|^2)}},\qquad
\mu_{aL}=\frac{2H_{a0}}{\sqrt{2(1+|\epsilon_K|^2)}}.
$$

The denominators here are $1-\epsilon_K^2$, not $1-|\epsilon_K|^2$. A non-Hermitian effective Hamiltonian does not generally imply $\mu_{ai}=\mu_{ia}^*$ or $\mu_{ai}=\mu_{ia}$. For a Hermitian dispersive mixing matrix in the orthonormal flavor basis, $H^{\rm disp}_{a0}=(H^{\rm disp}_{0a})^*$ and $H^{\rm disp}_{a\bar0}=(H^{\rm disp}_{\bar0a})^*$; this relation need not hold for the transformed entries in the nonorthogonal kaon basis. Absorptive contributions make the full effective Hamiltonian non-Hermitian. Rephasing or rescaling the basis changes the individual entries, while their products and the poles remain unchanged. These entries describe propagation and are distinct from final-state decay-amplitude projections: the latter use the physical bras with $p^*,q^*$, as in the kaon projection record.

## Step 2: three-pole equation

In the ordered right-eigenstate basis $(K_S,K_L,a)$ the assumed effective Hamiltonian and its characteristic equation are

$$
H=\begin{pmatrix}
\lambda_S&0&\mu_{Sa}\\
0&\lambda_L&\mu_{La}\\
\mu_{aS}&\mu_{aL}&\lambda_a
\end{pmatrix},
$$

$$
\det(zI-H)=(z-\lambda_S)(z-\lambda_L)(z-\lambda_a)
-(z-\lambda_L)\mu_{Sa}\mu_{aS}
-(z-\lambda_S)\mu_{La}\mu_{aL}=0.
$$

This polynomial remains valid at the unmixed kaon poles. Away from $z=\lambda_S,\lambda_L$, elimination of the kaon components gives the equivalent equation

$$
z-\lambda_a-\frac{\mu_{Sa}\mu_{aS}}{z-\lambda_S}
-\frac{\mu_{La}\mu_{aL}}{z-\lambda_L}=0.
$$

The denominator form alone can miss decoupled poles and must not be used at its singular points.

## Step 3: weak mixing and right eigenvectors

Scale every directed entry by a common formal parameter and expand about three distinct baseline poles. For a fixed basis normalization, assume both directed entries are small compared with the corresponding complex kaon--ALP separation, and the induced shifts are small compared with all relevant separations. The pole positions through second order are

$$
\begin{aligned}
\lambda_S^{\rm phys}&=\lambda_S+\frac{\mu_{Sa}\mu_{aS}}{\lambda_S-\lambda_a}+O(\mu^4),\\
\lambda_L^{\rm phys}&=\lambda_L+\frac{\mu_{La}\mu_{aL}}{\lambda_L-\lambda_a}+O(\mu^4),\\
\lambda_a^{\rm phys}&=\lambda_a+\frac{\mu_{Sa}\mu_{aS}}{\lambda_a-\lambda_S}
+\frac{\mu_{La}\mu_{aL}}{\lambda_a-\lambda_L}+O(\mu^4).
\end{aligned}
$$

Here $O(\mu^n)$ denotes order $n$ in that common scaling, not a new physical coupling. The corresponding right eigenvectors, with the dominant coordinate fixed to one, are

$$
\begin{aligned}
|K_S^{\rm phys}\rangle&=|K_S\rangle+\frac{\mu_{aS}}{\lambda_S-\lambda_a}|a\rangle+O(\mu^2),\\
|K_L^{\rm phys}\rangle&=|K_L\rangle+\frac{\mu_{aL}}{\lambda_L-\lambda_a}|a\rangle+O(\mu^2),\\
|a^{\rm phys}\rangle&=|a\rangle+\frac{\mu_{Sa}}{\lambda_a-\lambda_S}|K_S\rangle
+\frac{\mu_{La}}{\lambda_a-\lambda_L}|K_L\rangle+O(\mu^2).
\end{aligned}
$$

In particular, near a kaon--ALP degeneracy these expansions fail; near a $K_S$--$K_L$ degeneracy induced kaon-to-kaon mixing also requires care at second order. Left eigenvectors are independent and enter pole residues; these right-coordinate coefficients alone are not probabilities. For a complex pole shift $\delta\lambda$, the mass and width shifts are $\operatorname{Re}\delta\lambda$ and $-2\operatorname{Im}\delta\lambda$.

## Step 4: CP-conserving two-level limit

Assume the full effective Hamiltonian conserves CP and $a$ has definite CP. Only the kaon with matching CP can couple to it. In the conventional phase choice, $K_S$ and $K_L$ tend to the CP-even and CP-odd kaon states, respectively; a CP-odd ALP therefore couples to $K_L$ in this limit. The opposite-CP pole is unchanged. With the matching state denoted $K_i$, $i=S$ or $L$, setting both directed entries of the other kaon to zero gives

$$
\lambda_\pm=\frac{\lambda_i+\lambda_a}{2}\pm\frac12
\sqrt{(\lambda_i-\lambda_a)^2+4\mu_{ia}\mu_{ai}},
\qquad
(\lambda_+-\lambda_-)^2=(\lambda_i-\lambda_a)^2+4\mu_{ia}\mu_{ai}.
$$

The square-root branch exchanges the two labels; continuity fixes a label away from degeneracies. A right eigenvector in the two-level basis can be chosen as $(\mu_{ia},\lambda_\pm-\lambda_i)^T$ when nonzero. If the two unmixed complex energies coincide, the complex pole separation is $2\sqrt{\mu_{ia}\mu_{ai}}$ up to its branch, rather than necessarily a real mass splitting.

As an additional explicit assumption, for real symmetric dispersive mixing $\mu_{ia}=\mu_{ai}=\mu$ and equal unmixed masses, the squared complex separation becomes

$$
(\lambda_+-\lambda_-)^2=4\mu^2-\frac{(\Gamma_i-\Gamma_a)^2}{4}.
$$

Thus widths can convert a real mass splitting into a width splitting. At $|\Gamma_i-\Gamma_a|=4|\mu|$ with nonzero $\mu$, the two poles and their right eigenvectors coalesce: this exceptional point requires the full matrix and does not admit two independent pole eigenvectors. If both ALP couplings vanish simultaneously at degeneracy, the matrix is instead diagonal and the degeneracy is not exceptional. With CP violation, either in the kaon baseline or in the ALP mixing, both kaon entries can be nonzero; use the three-pole equation. Their values depend on the phase convention and on the underlying dynamics, so their mere presence is not by itself a convention-independent measure of CP violation.

## Verification and limits

The executable derivation is [the SymPy script](../../scripts/derive_kaon_alp_nonhermitian_poles.py); its results above are from its output. It constructs the similarity transformation, general flavor-directed entries, finite-$\epsilon_K$ CP-odd expressions, characteristic polynomial, two-level roots, perturbative pole shifts and right eigenvectors. Individual unit checks verify the general transformation and all four finite-$\epsilon_K$ CP-odd formulas. It prints the steps of the derivation, including the Schur complement obtained from the inverse kaon block.

```bash
.venv/bin/python scripts/derive_kaon_alp_nonhermitian_poles.py
```

The calculation uses constant widths and a finite-dimensional effective Hamiltonian. It does not include energy-dependent self-energies, thresholds, residues from a full propagator, ultraviolet matching, or a numerical estimate of the $|\Delta S|=1$ entries. It does not determine decay-channel rates or experimental time acceptance.
