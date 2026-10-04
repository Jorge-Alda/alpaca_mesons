---
title: Timelike propagator and omega-phi mixing conventions in RChT and dispersive analyses
status: literature
created: 2026-10-04
updated: 2026-10-04
bibtex_keys: [Ruiz-Femenia:2003jdx, Ecker:1989yg, Niecknig:2012sj, BaBar:2004ytv]
provenance: literature
---

# Timelike Propagator and $\omega$-$\phi$ Mixing Conventions

## Ruiz-Femenia:2003jdx — Ideal Mixing and Breit-Wigner Propagators

**Vector nonet and mixing (source lines 784--797):**
The physical $\omega$ and $\phi$ are superpositions of octet ($\omega_8$) and singlet ($\omega_1$) components. Assuming **ideal mixing**:

$$
|\omega\rangle = \sqrt{\frac{2}{3}}|\omega_1\rangle + \sqrt{\frac{1}{3}}|\omega_8\rangle,
\qquad
|\phi\rangle = -\sqrt{\frac{1}{3}}|\omega_1\rangle + \sqrt{\frac{2}{3}}|\omega_8\rangle.
$$

**Propagators in $\omega\to\pi\gamma$ and $\omega\to 3\pi$:**
The amplitudes use a relativistic Breit-Wigner propagator for the intermediate $\rho$ [@Ruiz-Femenia:2003jdx, Eq. (29), line 827--830]:

$$
\frac{1}{M_V^2 - s_{ij}},
$$

where $M_V$ is the **octet mass** (not the physical $\rho$ mass). For the external $\omega/\phi$, the mass $M_\omega$ or $M_\phi$ is used in kinematic factors.

**No dynamical $\omega$-$\phi$ mixing in propagators.** The mixing is fixed to ideal at the level of the external states. The relative phase between $\omega$ and $\phi$ amplitudes is not predicted by the short-distance matching.

## Ecker:1989yg — Antisymmetric Tensor Normalization

Ecker:1989yg establishes the antisymmetric tensor representation for massive spin-1 fields [@Ecker:1989yg, Sec. 1]. The kinetic term is [@Ecker:1989yg, Eq. (2.3) in spin-one-representations.md]:

$$
\mathcal L_{\mathrm{Kin}}(V) = -\frac{1}{2}\langle \nabla^\lambda V_{\lambda\mu} \nabla_\nu V^{\nu\mu} - \frac{M_V^2}{2} V_{\mu\nu}V^{\mu\nu}\rangle.
$$

**External state normalization:** The antisymmetric tensor field $V_{\mu\nu}$ creates/annihilates vector mesons with the normalization [@Ecker:1989yg, spin-one-representations.md]:

$$
V_{\mu\nu} = \frac{p_\mu \varepsilon_\nu - p_\nu \varepsilon_\mu}{M_V},
$$

where $\varepsilon_\mu$ is the polarization vector. This is the convention used in Ruiz-Femenia:2003jdx (source lines 823--824, 914--915) and in the local numerical kernel (`odd-parity-rcht-numerical-kernel.md` cites Ecker:1989yg for external-state normalization).

**Internal vs. external propagator:** The Ecker normalization applies to **external** vector states (on-shell). For **internal** (timelike) propagators in loops or exchange diagrams, the propagator is derived from the kinetic term:

$$
\Delta_{\mu\nu,\rho\sigma}(p) = \frac{1}{M_V^2 - p^2} \left( g_{\mu\rho}g_{\nu\sigma} - g_{\mu\sigma}g_{\nu\rho} \right) + \text{gauge terms}.
$$

The local kernel uses this for internal vector exchange. The Ecker normalization does **not** directly apply to the internal timelike propagator; the internal propagator is the standard Feynman propagator from the kinetic term.

## Niecknig:2012sj — Dispersive Propagator (Omnès Function)

The dispersive approach does not use a Breit-Wigner propagator. Instead, the $\pi\pi$ P-wave rescattering is resummed into the **Omnès function** [@Niecknig:2012sj, Eq. (12)]:

$$
\Omega(s) = \exp\left\{ \frac{s}{\pi} \int_{4m_\pi^2}^\infty ds' \frac{\delta(s')}{s'(s'-s)} \right\},
$$

where $\delta(s)$ is the $\pi\pi$ P-wave phase shift. The $\rho$ resonance emerges dynamically from the phase shift. The amplitude is [@Niecknig:2012sj, Eq. (14)]:

$$
\mathcal F(s) = \Omega(s) \left\{ a + \frac{s}{\pi} \int_{4m_\pi^2}^\infty \frac{ds'}{s'} \frac{\sin\delta(s')\hat{\mathcal F}(s')}{|\Omega(s')|(s'-s)} \right\}.
$$

**No explicit $\omega$-$\phi$ mixing.** The $\omega$ and $\phi$ decays are treated in parallel with the same formalism; only the mass $M_V$ and normalization (partial width) differ [@Niecknig:2012sj, Sec. 7, lines 1352--1360]. The $e^+e^-\to 3\pi$ cross section would use the same formalism with $s_{e^+e^-}$-dependent normalization.

## BaBar:2004ytv — Fit Model with Fixed Relative Phase

The BaBar fit uses a coherent sum of four Breit-Wigner resonances [@BaBar:2004ytv, Eq. (14)]:

$$
\sigma_{3\pi}(m) = \frac{12\pi}{m^3}F_{\rho\pi}(m)
\left| \sum_{V=\omega,\phi,\omega',\omega''}
\frac{\Gamma_V m_V^{3/2}\sqrt{\mathcal B(V\to e^+e^-)\mathcal B(V\to 3\pi)}}{D_V(m)}
\frac{e^{i\phi_V}}{\sqrt{F_{\rho\pi}(m_V)}} \right|^2,
$$

with $D_V(m) = m_V^2 - m^2 - i m \Gamma_V(m)$ and mass-dependent widths $\Gamma_V(m)$ from SND2003 for $\omega,\phi$ [@BaBar:2004ytv, Sec. 6, lines 1684--1690].

**$\omega$-$\phi$ mixing:** The relative phase $\phi_\phi = (163\pm7)^\circ$ is **fixed from SND2003** [@BaBar:2004ytv, Sec. 6, lines 1767--1768]. Phases of $\omega,\omega',\omega''$ are fixed at $0^\circ, 180^\circ, 0^\circ$ [@BaBar:2004ytv, Sec. 6, line 1769--1770, citing Clegg]. No dynamical mixing matrix is used.

## Summary of Conventions

| Aspect | Ruiz-Femenia (RChT) | Ecker:1989yg | Niecknig (Dispersive) | BaBar (Fit) |
|--------|---------------------|--------------|----------------------|-------------|
| **Vector representation** | Antisymmetric tensor $V_{\mu\nu}$ | Antisymmetric tensor $V_{\mu\nu}$ | Not used (dispersive) | Not used (phenomenological) |
| **External state norm.** | $V_{\mu\nu} = (p_\mu\varepsilon_\nu-p_\nu\varepsilon_\mu)/M_V$ | Same (source) | N/A | N/A |
| **Internal propagator** | $1/(M_V^2-s)$ from kinetic term | From kinetic term | Omnès function $\Omega(s)$ | Breit-Wigner $D_V(m)$ |
| **$\omega$-$\phi$ mixing** | Ideal mixing (fixed) | Not discussed | Not in formalism | Fixed phase $\phi_\phi=163^\circ$ |
| **Mass-dependent width** | $\Gamma_V(m)$ for $\omega,\phi$ from SND2003 | Not discussed | Implicit in phase shift | $\Gamma_V(m)$ from SND2003 |

## Implications for RChT Predictor

1. **External vector normalization:** The local kernel's use of Ecker antisymmetric-tensor normalization for external states is consistent with Ruiz-Femenia. This normalization must be used when computing $V\to 3\pi$ amplitudes from the RChT Lagrangian.

2. **Internal timelike propagator:** For $\rho$ exchange in $V\to 3\pi$, the internal propagator is $1/(M_V^2-s_{ij})$ with the octet mass $M_V$. This is what the local kernel implements.

3. **$\omega$-$\phi$ mixing in $e^+e^-\to 3\pi$:** The BaBar fit uses a fixed relative phase. In RChT, the $\omega$-$\phi$ mixing should be implemented at the level of the vector nonet diagonalization (photon-vector mixing), not as an ad hoc phase. The local kernel's `photon_vector_transform` handles this.

4. **Dispersive vs. Breit-Wigner:** The dispersive Omnès function resums $\pi\pi$ rescattering to all orders, while the Breit-Wigner is a single-pole approximation. For precision line shapes (especially $\phi\to 3\pi$ where the $\rho$ peak is in the physical region), the dispersive approach is superior. The RChT predictor should use the dispersive $\pi\pi$ phase shift input for the $\rho$ propagator, not a fixed-width Breit-Wigner.