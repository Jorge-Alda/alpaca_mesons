---
title: Double-trace spurion extension and mass-compatible vector branches
status: derived
created: 2026-10-13
updated: 2026-10-13
bibtex_keys: [Ecker:1988te, Alda:2026nrd]
provenance: analytical-derivation
---

# Double-trace spurion extension and mass-compatible vector branches

## Inputs and assumptions

Use the canonical tensor normalization, positive-strange ideal phi phase, vacuum spurion $\chi_+=\operatorname{diag}(x_l,x_l,x_s)$ and strong-isospin assumptions of [the minimal mass action](vector-nonideal-mass-action.md), with its cited kinetic conventions [@Ecker:1988te]. Numerical spurion inputs must follow the downloaded arXiv version of Alda:2026nrd [@Alda:2026nrd], rather than a local manuscript. This derivation requires no numerical values. It retains tree-level canonical kinetic terms; loops, electromagnetic shifts and higher spurion orders remain omitted.

**Assumed extension:** add one symmetry-allowed double-trace spurion mass operator,

$$
\mathcal L_{\rm extra}=\frac{\kappa_V}{4}\langle V_{\mu\nu}\rangle\langle V^{\mu\nu}\chi_+\rangle.
$$

This equation defines the dimensionless $\kappa_V$ normalization exactly. Its coefficient is a new independent parameter, not a value inferred from Alda:2026nrd or a claim that the operator exhausts all breaking terms. The Lorentz contraction is suppressed only during flavor mass extraction. Every mass and spurion entry has dimension two in natural units.

## Executable calculation

Run [the script](../../../scripts/derive_vector_mass_extension.py):

```text
.venv/bin/python scripts/derive_vector_mass_extension.py
```

SymPy constructs the established explicit ideal generators, evaluates the two flavor traces, extracts the Hessian of the quadratic action, calibrates the charged masses, solves the neutral trace equation, and derives the determinant polynomial and discriminant. Its `derive()` function exposes these intermediate symbolic expressions for reuse. The calculation gives

$$
\mathcal M_V^2=\begin{pmatrix}
M_V^2-4e_m^Vx_l+2\Delta_V/3+2\kappa_Vx_l & \sqrt2\Delta_V/3+\kappa_V(x_l+x_s)/\sqrt2\\
\sqrt2\Delta_V/3+\kappa_V(x_l+x_s)/\sqrt2 & M_V^2-4e_m^Vx_s+\Delta_V/3+\kappa_Vx_s
\end{pmatrix}.
$$

The charged traceless generators have $\langle T\rangle=0$, so their masses remain

$$
r\equiv m_\rho^2=M_V^2-4e_m^Vx_l,\qquad
k\equiv m_{K^*}^2=M_V^2-2e_m^V(x_l+x_s).
$$

For $x_s\ne x_l$, these fix $e_m^V=(r-k)/[2(x_s-x_l)]$ and $M_V^2=r+4e_m^Vx_l$. Let $W=m_\omega^2$, $P=m_\phi^2$, $T=W+P-2k$, and $\delta x=x_s-x_l$. The trace fixes

$$
\Delta_V=T-\kappa_V(2x_l+x_s).
$$

After that substitution the neutral entries are

$$
A=r+\frac{2T}{3}-\frac{2\kappa_V\delta x}{3},\quad
D=2k-r+\frac{T}{3}+\frac{2\kappa_V\delta x}{3},\quad
B=\frac{\sqrt2T}{3}+\frac{\kappa_V\delta x}{3\sqrt2}.
$$

The remaining equation is $AD-B^2-WP=0$, or

$$
-\frac{(\delta x)^2}{2}\kappa_V^2-\frac{4(k-r)\delta x}{3}\kappa_V+C=0,
\qquad
C=-WP+\frac{(W+P)(4k-r)}3-\frac{8k^2}3+\frac{8kr}3-r^2.
$$

SymPy obtains the discriminant

$$
\mathcal D=-\frac{2(\delta x)^2}{9}(3P-4k+r)(3W-4k+r).
$$

## Branching and identifiability

For $W<P$ and $\delta x\ne0$, real mass-compatible branches exist exactly when

$$
W\le\frac{4k-r}{3}\le P.
$$

A strict inequality gives two distinct parameter solutions. Equality gives one repeated solution. A negative discriminant rules out an exact real fit within this extended truncation, without ruling out more general vector actions. The two solutions are

$$
\kappa_V^{(\pm)}=\frac{-4(k-r)/3\pm\sqrt{-2(3P-4k+r)(3W-4k+r)/9}}{\delta x},
\qquad
\Delta_V^{(\pm)}=T-\kappa_V^{(\pm)}(2x_l+x_s).
$$

Their angles must be calculated separately with the established convention $\theta_V=\tfrac12\operatorname{atan2}(2B,D-A)$, assigning the lower eigenvalue to omega. Masses alone generally do not select a unique physical generator. Electromagnetic residues and a declared interaction-breaking prescription can discriminate branches; preference for a small angle or a small coefficient is an additional assumption, not a mass measurement. Positivity of $W,P$ ensures neutral mass stability at an exact solution, but gives no control over omitted interactions or truncation errors.

The symmetric limit $x_l=x_s$ requires separate treatment: the charged calibration is singular unless $k=r$, and the two double-trace structures become proportional. Then only a combination of $\Delta_V$ and $\kappa_V$ is identifiable. A degenerate neutral spectrum similarly leaves the physical basis undetermined.

## Verification

The executable script passed. Assertions independently verify the neutral trace, both determinant roots and the quadratic coefficient. No measured mass or decay input is fitted in this record; numerical calibration and uncertainty choices belong to the stage-3 calibration record.
