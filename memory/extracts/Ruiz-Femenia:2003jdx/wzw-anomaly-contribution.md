---
title: Anomalous WZW contribution to V -> 3pi and e+e- -> 3pi
status: literature
created: 2026-10-04
updated: 2026-10-26
bibtex_keys: [Wess:1971yu, Ruiz-Femenia:2003jdx, Niecknig:2012sj]
provenance: literature
---

Displayed tensors and odd couplings use $\epsilon^{0123}=+1$ and
$\epsilon_{0123}=-1$ for $(+---)$. The canonical coefficients are
$c_i=-c_i^R$, $d_i=-d_i^R$ relative to the opposite source convention.
The [conversion audit](../../rcht/odd-parity/levi-civita-convention.md)
derives the action-preserving map; archived source text remains literal.


# Anomalous / WZW Contribution to $V\to 3\pi$ and $e^+e^-\to 3\pi$

## WZW Term in Ruiz-Femenia:2003jdx

The Wess-Zumino action $Z_{\mathrm{WZ}}[v,a]$ is the $\mathcal O(p^4)$ term that fulfills the chiral anomaly [@Ruiz-Femenia:2003jdx, Sec. 2, lines 317--329]. The relevant piece for the $\langle VVP\rangle$ Green's function is the pseudoscalar decay into two photons [@Ruiz-Femenia:2003jdx, Eq. (17), lines 662--668]:

$$
\mathcal L^{(4)}_{\mathrm{WZ}}=\frac{\sqrt{2} N_C}{8\pi^2 F}
\epsilon_{\mu\nu\alpha\beta}
\langle \Phi \,\partial^{\mu}v^{\nu}\,
\partial^{\alpha}v^{\beta}\,
\rangle.
$$

This term contributes to the $\langle VVP\rangle$ three-point function as the last term in Eq. (18) [@Ruiz-Femenia:2003jdx, line 652]:

$$
\frac{N_C}{8\pi^2 r^2}.
$$

**Correction:** The gauged WZW functional does contain a direct one-photon, three-pion anomalous vertex. In the convention $U=\exp(i\Phi/f_\pi)$ and $v_\mu=-eQA_\mu$, its magnitude is
$$
\mathcal L_{\gamma3\pi}
=\frac{N_C e}{12\pi^2 f_\pi^3}
\epsilon^{\mu\nu\alpha\beta}A_\mu
\partial_\nu\pi^+\partial_\alpha\pi^-\partial_\beta\pi^0,
$$
up to the common WZW orientation and field-phase convention. This is the low-energy anomaly vertex [@Wess:1971yu]. The earlier statement that the WZW term does not generate $\gamma3\pi$ was incorrect.

In antisymmetric-tensor RChT, however, $V_{\mu\nu}$ is a separate resonance field and $Z_{\rm WZ}[v,a]$ depends on external sources. Diagonalizing photon--vector mixing inserts a vector-field component into the photon appearing in this vertex. The resulting local $V_i3\pi$ term is proportional to $e(g_i/M_i)$ and hence is $O(e^2)$ when $g_i=O(e)$; it is an electromagnetic, field-basis-induced contribution, not the leading strong decay vertex. Keep it distinct from any independent chirally invariant odd-parity resonance operator and do not double count it against explicit photon exchange. See the [derived WZW contribution record](../../rcht/odd-parity/wzw-gamma3pi-contribution.md).

**For $e^+e^-\to 3\pi$:** The process is initiated by $e^+e^-\to\gamma^*$. Its low-energy current amplitude contains the direct WZW $\gamma^*\to3\pi$ term above, together with resonance and rescattering contributions at higher energies. Ruiz-Femenia:2003jdx does not assemble a full $e^+e^-\to3\pi$ amplitude; this does not remove the WZW contact fixed by the anomaly.

## $\phi\to 3\pi$ and the Anomaly in Niecknig:2012sj

Niecknig:2012sj uses a dispersive framework based on $\pi\pi$ P-wave phase shifts [@Niecknig:2012sj, Sec. 1]. The amplitude is decomposed as [@Niecknig:2012sj, Eq. (1)]:

$$
\mathcal M(s,t,u)=i\epsilon_{\mu\nu\alpha\beta}n^\mu p_+^\nu p_-^\alpha p_0^\beta \mathcal F(s,t,u),
$$

where $n^\mu$ is the polarization vector of the decaying vector meson. The dispersive representation **does not include an explicit WZW/anomaly term**; the amplitude is fully determined by the $\pi\pi$ P-wave phase shift (elastic unitarity) and subtraction constants fixed by the partial width.

The paper notes that the $K\bar K$ intermediate state in $\phi\to 3\pi$ factorizes and can be absorbed into a wave-function renormalization constant (an unobservable overall phase) [@Niecknig:2012sj, Sec. 2, lines 253--264]. The $\rho\pi$ state is **not** a distinct two-body intermediate state in dispersion theory; the $\rho$ resonance appears as the $\pi\pi$ P-wave interaction [@Niecknig:2012sj, lines 267--272].

## Is the Anomaly Required for a Quantitative $\phi\to 3\pi$ Line Shape?

**No separate WZW term is added to the strong dispersive $\phi\to 3\pi$ decay amplitude in Niecknig:2012sj.** That analysis is a strong decay representation; its line shape is determined by:
1. The $\pi\pi$ P-wave phase shift (input from Roy equations)
2. Crossed-channel rescattering (three-particle cuts)
3. Subtraction constants (normalization from $\Gamma_{\phi\to 3\pi}$)

The dispersive representation achieves $\chi^2/\mathrm{ndof} \approx 1.17$ (KLOE) and $\approx 1.0$ (twice-subtracted) without a separately displayed anomaly term [@Niecknig:2012sj, Tables 1 and 2]. This does not imply that the $\gamma3\pi$ anomaly is absent or encoded in the $\pi\pi$ phase shift. The strong-decay subtraction constants are fixed from the decay normalization; the electromagnetic $\gamma^*\to3\pi$ dispersive amplitude requires its own low-energy normalization, conventionally matched to the anomaly.

**Conclusion:** Keep three statements separate: the WZW functional fixes the low-energy $\gamma3\pi$ anomalous vertex; photon--vector diagonalization induces an electromagnetic $V_i3\pi$ term in the chosen field basis; and a strong $V_i\to3\pi$ resonance amplitude may also contain independent odd-parity local operators and rescattering. The $\langle VVP\rangle$ short-distance constraints on $c_i,d_i$ do not by themselves complete that one-vector, three-pion operator basis.
