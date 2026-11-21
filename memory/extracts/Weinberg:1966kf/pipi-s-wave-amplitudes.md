---
title: Leading-order pi pi amplitude, isospin channels and S-wave normalization (Weinberg:1966kf)
status: mixed
created: 2026-11-21
updated: 2026-11-21
bibtex_keys: [Weinberg:1966kf, Chivukula:1992gi]
provenance: literature
---

# Leading-order $\pi\pi$ S-wave amplitudes

## Access and verification status

- Full text of Weinberg:1966kf (Phys. Rev. Lett. 17 (1966) 616) was not accessible (no arXiv version; APS and OSTI pages give no text). No source archive exists under `memory/sources/Weinberg:1966kf/`.
- Verified from the INSPIRE-HEP abstract: current commutation relations and the PCAC assumption determine the matrix elements for emission and absorption of soft pions, hence the scattering length of a pion on any target; the paper gives a formula for pion scattering on any particle except a pion, and then extends the result to $\pi\pi$ scattering.
- The equations below are NOT quoted from Weinberg:1966kf. They are the leading-order chiral amplitudes written in the conventions of Chivukula:1992gi (source read directly: arXiv hep-ph/9206222, Sections 1-2), specialised to $N=2$. They are listed here because they are the standard content of the soft-pion $\pi\pi$ result. Any statement that Weinberg:1966kf contains these exact equations or numerical values is unverified.

## Leading-order amplitude (Chivukula:1992gi, general $N$, $N=2$ case)

With $\Sigma=\exp(2i\pi^aT^a/f)$, $\mathrm{tr}\,T^aT^b=\delta^{ab}/2$, and $\mathcal{L}_2$ the leading Lagrangian with a common Goldstone mass $m$, Chivukula:1992gi gives (text after Eq. (2.3)):

$$A(s,t,u)=\frac{2}{N}\frac{s-m^2}{f^2},\qquad B(s,t,u)=\frac{s-m^2}{f^2}.$$

For $N=2$ there are no $d^{abc}$ symbols, so $B$ drops out and

$$A(s,t,u)=\frac{s-m^2}{f^2},\qquad f\simeq 93\ \mathrm{MeV}\ \text{(pion decay-constant normalization)}.$$

## Isospin channels ($N=2$)

From the channel decomposition in Chivukula:1992gi (Eq. (2.3), channels $\Delta,F,X$ for isospin $0,1,2$) with $N=2$:

$$a_{\Delta}=3A(s,t,u)+A(t,s,u)+A(u,t,s),\quad a_{F}=A(t,s,u)-A(u,t,s),\quad a_{X}=A(t,s,u)+A(u,t,s).$$

Inserting $A=(s-m^2)/f^2$ and $s+t+u=4m^2$:

$$a_\Delta=\frac{2s-m^2}{f^2},\qquad a_F=\frac{t-u}{f^2},\qquad a_X=\frac{2m^2-s}{f^2}.$$

(The last three forms are an elementary substitution by the present memory, labelled derived; the preceding projections are from the source.)

## Partial-wave normalization (Chivukula:1992gi, Eq. (2.4))

$$a_{I\ell}(s)=\frac{1}{64\pi}\int_{-1}^{1}a_I(s,\cos\theta)\,P_\ell(\cos\theta)\,d\cos\theta.$$

Resulting S-wave and P-wave amplitudes (Chivukula:1992gi, Eq. (2.5), $N=2$):

$$a_{\Delta0}=\frac{2s-m^2}{32\pi f^2},\qquad a_{X0}=\frac{2m^2-s}{32\pi f^2},\qquad a_{F1}=\frac{s-4m^2}{96\pi f^2}.$$

These are the source's general-$N$ formulas, e.g. $a_{\Delta0}=\frac{Ns}{32\pi f^2}-\frac{m^2}{16\pi N f^2}$, evaluated at $N=2$; I checked the $N=2$ reduction of $a_{X0}$ and $a_{F1}$ by hand against the projection above.

## Threshold values (derived, not quoted from Weinberg:1966kf)

At $s=4m^2$:

$$a_{\Delta0}\big|_{\rm thr}=\frac{7m^2}{32\pi f^2},\qquad a_{X0}\big|_{\rm thr}=-\frac{2m^2}{32\pi f^2}=-\frac{m^2}{16\pi f^2}.$$

These coincide with the widely quoted soft-pion scattering lengths, but the attribution of these numbers to the text of Weinberg:1966kf is unverified here.

## Normalization remark

With the $1/(64\pi)$ convention above, $a_{I\ell}$ equals the standard partial wave $\frac{1}{32\pi}\int_{0}^{1}\ldots$ for $\pi\pi$ scattering, so that the elastic unitarity circle has $|\mathrm{Re}\,a|\le 1/2$ (derived remark; the source speaks only of the Argand circle and of amplitudes "exceeding 1", see the Chivukula:1992gi extract).
