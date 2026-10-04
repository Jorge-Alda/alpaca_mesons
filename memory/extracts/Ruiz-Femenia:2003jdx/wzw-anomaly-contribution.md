---
title: Anomalous WZW contribution to V -> 3pi and e+e- -> 3pi
status: literature
created: 2026-10-04
updated: 2026-10-04
bibtex_keys: [Ruiz-Femenia:2003jdx, Niecknig:2012sj]
provenance: literature
---

# Anomalous / WZW Contribution to $V\to 3\pi$ and $e^+e^-\to 3\pi$

## WZW Term in Ruiz-Femenia:2003jdx

The Wess-Zumino action $Z_{\mathrm{WZ}}[v,a]$ is the $\mathcal O(p^4)$ term that fulfills the chiral anomaly [@Ruiz-Femenia:2003jdx, Sec. 2, lines 317--329]. The relevant piece for the $\langle VVP\rangle$ Green's function is the pseudoscalar decay into two photons [@Ruiz-Femenia:2003jdx, Eq. (17), lines 662--668]:

$$
\mathcal L^{(4)}_{\mathrm{WZ}}=-\frac{\sqrt{2} N_C}{8\pi^2 F}
\epsilon_{\mu\nu\alpha\beta}
\langle \Phi \,\partial^{\mu}v^{\nu}\,
\partial^{\alpha}v^{\beta}\,
\rangle.
$$

This term contributes to the $\langle VVP\rangle$ three-point function as the last term in Eq. (18) [@Ruiz-Femenia:2003jdx, line 652]:

$$
-\frac{N_C}{8\pi^2 r^2}.
$$

**For $V\to 3\pi$:** The WZW term does **not** directly generate a $V\to 3\pi$ vertex. The WZW term has one pseudoscalar and two vector currents; to get $V\to 3\pi$, one would need to replace one vector current by two pseudoscalars (via the chiral Lagrangian), which is a higher-order effect. The paper does not write a WZW-induced $V\to 3\pi$ amplitude.

**For $e^+e^-\to 3\pi$:** The process $e^+e^-\to 3\pi$ proceeds via $e^+e^-\to \gamma^* \to V \to 3\pi$ (with $V=\omega,\phi,\omega',\omega''$). The WZW term contributes to the $\gamma^*\to 3\pi$ vertex only through the anomaly-induced $\pi^0\gamma^*\gamma$ coupling (which is not the $3\pi$ final state) or through the $\langle VVP\rangle$ function with one vector current being the photon. The paper does not give an explicit $e^+e^-\to 3\pi$ amplitude from WZW.

## $\phi\to 3\pi$ and the Anomaly in Niecknig:2012sj

Niecknig:2012sj uses a dispersive framework based on $\pi\pi$ P-wave phase shifts [@Niecknig:2012sj, Sec. 1]. The amplitude is decomposed as [@Niecknig:2012sj, Eq. (1)]:

$$
\mathcal M(s,t,u)=i\epsilon_{\mu\nu\alpha\beta}n^\mu p_+^\nu p_-^\alpha p_0^\beta \mathcal F(s,t,u),
$$

where $n^\mu$ is the polarization vector of the decaying vector meson. The dispersive representation **does not include an explicit WZW/anomaly term**; the amplitude is fully determined by the $\pi\pi$ P-wave phase shift (elastic unitarity) and subtraction constants fixed by the partial width.

The paper notes that the $K\bar K$ intermediate state in $\phi\to 3\pi$ factorizes and can be absorbed into a wave-function renormalization constant (an unobservable overall phase) [@Niecknig:2012sj, Sec. 2, lines 253--264]. The $\rho\pi$ state is **not** a distinct two-body intermediate state in dispersion theory; the $\rho$ resonance appears as the $\pi\pi$ P-wave interaction [@Niecknig:2012sj, lines 267--272].

## Is the Anomaly Required for a Quantitative $\phi\to 3\pi$ Line Shape?

**No explicit WZW/anomaly term appears in the dispersive $\phi\to 3\pi$ amplitude.** The line shape is determined by:
1. The $\pi\pi$ P-wave phase shift (input from Roy equations)
2. Crossed-channel rescattering (three-particle cuts)
3. Subtraction constants (normalization from $\Gamma_{\phi\to 3\pi}$)

The dispersive representation achieves $\chi^2/\mathrm{ndof} \approx 1.17$ (KLOE) and $\approx 1.0$ (twice-subtracted) **without an explicit anomaly term** [@Niecknig:2012sj, Tables 1 and 2]. The anomaly is implicitly contained in the $\pi\pi$ phase shift input (which comes from data/Roy equations that include all QCD effects).

**Conclusion:** For a quantitative $\phi\to 3\pi$ line shape, the dispersive approach uses the empirical $\pi\pi$ phase shift as input, which already encodes the full QCD dynamics including anomaly effects. An explicit WZW term is not added separately. In RChT, the anomaly is matched onto the resonance couplings via the $\langle VVP\rangle$ short-distance constraints (which fix combinations of $c_i,d_i$), but the direct WZW contribution to $V\to 3\pi$ is not a separate piece in the amplitude.