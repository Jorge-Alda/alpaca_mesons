---
title: ALP-meson inverse propagator under a linear field change
status: derived
created: 2026-11-11
updated: 2026-11-11
bibtex_keys: [Alda:2026nrd]
provenance: analytical-derivation
---

# ALP-meson inverse propagator under a linear field change

## Assumptions and conventions

Follow `Alda:2026nrd` [@Alda:2026nrd] and the [finite-mixing record](finite-alp-meson-mixing.md). The quadratic Lagrangian of the neutral mesons and the ALP is
$$
\mathcal L_2=\frac12\partial_\mu\Phi^{\mathsf T}K\,\partial^\mu\Phi-\frac12\Phi^{\mathsf T}M^2\Phi,\qquad
\Phi=(P_i,a),\qquad D(s)=sK-M^2 .
$$
Assumptions: $K$ and $M^2$ are real symmetric constant matrices, $K$ is positive definite (true at LO, see the finite-mixing record), and the selected poles are simple. Only bilinears with at most two derivatives are kept, as in the restricted Lagrangian of the [partonic record](alp-general-chiral-redefinition-invariants.md); loops are not included. $K$ is dimensionless; $s$ and $M^2$ are squared masses.

The field change is $\Phi=T\Phi'$ with a constant invertible matrix $T$. The amputated vertex linear in the external leg is written $\mathcal V_I\Phi_I$, as in the finite-mixing record.

Notation (approved by the user): $v^{(k)}$, with components $v^{(k)}_I$, denotes a null vector of $D(m_k^2)$ in $(P_i,a)$ space, i.e. the field-space direction of pole $k$. When its ALP component is nonzero it is normalized as in `Alda:2026nrd`, $v^{(k)}=(w_i^{(k)}(m_k^2),1)$; the general symbol is needed when its ALP component vanishes or is rescaled. It is not the axial overlap $(f_\pi/f_a)\operatorname{tr}(c_q^A\lambda_i)$ of the finite-mixing record.

## Derivation

1. Substituting $\Phi=T\Phi'$ gives $K'=T^{\mathsf T}KT$ and $M^{2\prime}=T^{\mathsf T}M^2T$, hence
$$
D'(s)=T^{\mathsf T}D(s)\,T .
$$
The kinetic and mass matrices change separately; only their combination $D(s)$ transforms as one object. This is why removing kinetic mixing or mass mixing alone has no physical meaning.

2. Poles. $\det D'(s)=(\det T)^2\det D(s)$, so the zeros $s=m_k^2$ are the same.

3. Pole vectors. $D(m_k^2)v^{(k)}=0$ is equivalent to $D'(m_k^2)\,T^{-1}v^{(k)}=0$, so $v'^{(k)}=T^{-1}v^{(k)}$, i.e. $v^{(k)}=Tv'^{(k)}$.

4. Residues. $D'(s)^{-1}=T^{-1}D(s)^{-1}T^{-\mathsf T}$. Near a simple pole, $D(s)=D(m_k^2)+(s-m_k^2)K$ gives the residue
$$
\operatorname*{Res}_{s=m_k^2}D(s)^{-1}=\frac{v^{(k)}v^{(k)\mathsf T}}{v^{(k)\mathsf T}Kv^{(k)}},
$$
and the residue in the new coordinates is $T^{-1}(\cdots)T^{-\mathsf T}$. The pole normalization $v^{(k)\mathsf T}Kv^{(k)}$, which equals $K_{\hat a\hat a}+2K_{i\hat a}w_i^{(k)}+K_{ij}w_i^{(k)}w_j^{(k)}$ for $v^{(k)}=(w^{(k)},1)$, is unchanged.

5. Amplitudes. The vertex transforms as $\mathcal V'=T^{\mathsf T}\mathcal V$, so the LSZ one-leg amplitude of the finite-mixing record,
$$
\mathcal M=\frac{\mathcal V^{\mathsf T}v^{(k)}}{\sqrt{v^{(k)\mathsf T}Kv^{(k)}}},
$$
is unchanged when $v'^{(k)}=T^{-1}v^{(k)}$. If the primed pole vector is instead rescaled to unit ALP component, the amplitude can change sign only; this is the phase convention of the external state. The contraction $\mathcal V^{\mathsf T}\operatorname{Res}D^{-1}\,\mathcal V$ is invariant without any sign choice.

## Verification

[`derive_alp_inverse_propagator_congruence.py`](../../../scripts/derive_alp_inverse_propagator_congruence.py) uses exact rationals for $N_f=3$ neutral mesons on $(P_1,P_2,P_3,a)$. $K$ has the LO form $K_{ij}=\delta_{ij}$, $K_{i\hat a}=-(f_\pi/f_a)\operatorname{tr}(c_q^A\lambda_i)$, $K_{\hat a\hat a}=1+\sum_i[(f_\pi/f_a)\operatorname{tr}(c_q^A\lambda_i)]^2$, with rational values chosen directly for $(f_\pi/f_a)\operatorname{tr}(c_q^A\lambda_i)$; $M^2$ is an arbitrary symmetric matrix with rational poles, built from a rational orthogonal matrix (Cayley transform). The test data are not physical values. For a generic rational $T$ (three seeds in the test) it checks the congruence, $\det D'=(\det T)^2\det D$, the pole positions, the null vectors before and after, the residue formula and its transformation, the invariance of the pole normalization, of $\mathcal V^{\mathsf T}\operatorname{Res}D^{-1}\mathcal V$ and of the amplitude. All checks pass.

```sh
.venv/bin/python scripts/derive_alp_inverse_propagator_congruence.py
.venv/bin/python -m unittest discover -s scripts -p test_derive_alp_inverse_propagator_congruence.py
```

Related records: [single-pole alignment](alp-single-pole-alignment.md), [several poles](alp-multi-pole-alignment.md).
