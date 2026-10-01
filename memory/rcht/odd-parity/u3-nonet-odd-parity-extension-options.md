---
title: Decision record for a U(3) odd-parity nonet and ALP-source extension
status: mixed
created: 2026-09-28
updated: 2026-09-30
bibtex_keys: [Alda:2026bar, Bai:2025fvl, Fujiwara:1984mp, Ruiz-Femenia:2003jdx, Roig:2013baa, Kaiser:2000gs, Herrera-Siklody:1996tqr, Ecker:1988te, Ecker:1989yg]
provenance: literature
---

# Scope and current point of departure

The repository currently implements the seven $VJP$ and four $VVP$ single-trace monomials of Ruiz--Femenia et al. in antisymmetric-tensor notation [@Ruiz-Femenia:2003jdx]. The leading-$N_c$ statement applies to the displayed trace topology, not to an arbitrary $U(3)$ completion: the general anomalous RChT action is organized by resonance-field number, and the one-vector/two-vector sector is only one truncation [@Roig:2013baa].

## Concrete choices

**A. Minimal nonet continuation (recommended now).** Define a leading-OZI model by promoting $V$ and the Goldstone matrix to $U(3)$ and evaluating the existing seven/four single traces with vector singlet generator $T_0=\mathbf1/\sqrt3$, $\operatorname{tr}(T_0^2)=1$, and Alda pseudoscalar singlet $\lambda_0=\sqrt{2/3}\,\mathbf1$, $\operatorname{tr}(\lambda_0^2)=2$. Take the same $c_i,d_i$ for octet and singlet components, set independent multi-trace coefficients to zero, and omit additional $X_{TC}$ coefficient functions and derivative-$X_{TC}$ operators. Ideal $\omega$--$\phi$ mixing and physical masses are declared baseline approximations, not consequences of $U(3)$. Kaiser--Leutwyler's generating-functional organization associates one trace generically with one quark loop and motivates OZI suppression for additional trace topologies [@Kaiser:2000gs]. This is nominal coefficient/topology counting, not a proof of uniform numerical suppression after field normalization, singlet insertions, or channel projections. The common singlet/octet $c_i,d_i$ are a separate ansatz, not a consequence of trace counting, and this model is neither a symmetry theorem nor a complete $U(3)$ basis.

**B. Minimal plus controlled singlet/OZI fit.** Add the independent invariants obtained by replacing one trace with products of traces, including singlet-vector projections such as $\langle V_{\mu\nu}\rangle\langle\{V^{\mu\nu},\ldots\}\rangle$ and $\langle V_{\mu\nu}V^{\mu\nu}\rangle\langle\chi_-\rangle$ when the Lorentz/parity structure permits it. Use one nominal $1/N_c$ suppression per extra trace only as generic OZI prior bookkeeping, then account separately for field/source normalization and fit the combinations visible in the target channels. This rule is separate from anomalous-source counting in powers or derivatives of $X_{TC}$. This option is appropriate if eta-prime/singlet-vector data are part of the accuracy goal, but it expands the parameter space and loses the current short-distance closure unless new matching constraints are supplied [@Kaiser:2000gs; @Roig:2013baa].

**C. Full $U(3)$ source-function treatment.** Denote Herrera-Siklody et al.'s invariant by $X_{\rm HS}=\langle\log U\rangle+\hat\theta$. Alda Eq. (III.14) instead defines

$$
X_{TC}=2i\left(\langle\log u\rangle+i c_G\frac{a}{f_a}\right).
$$

On a common logarithm branch, $U=u^2$ relates the mesonic logarithms, but the recorded conventions do not derive the normalization or sign map between $\hat\theta$ and $c_Ga/f_a$. No exact factor relation between $X_{\rm HS}$ and $X_{TC}$ is therefore claimed. Translating a Herrera-Siklody coefficient function requires an independently stated source map [@Herrera-Siklody:1996tqr; @Alda:2026bar]. A complete treatment allows parity-appropriate coefficient functions of Alda's $X_{TC}$ and independent derivative-$X_{TC}$ structures. Such terms can contribute directly to dynamical-ALP amplitudes and are omitted in option A as an explicit large-$N_c$/operator truncation, not because topological observables are absent.

## Mandatory versus optional at the current accuracy

Mandatory for a consistent leading nonet extension: (i) a normalized $U(3)$ singlet in $V$ and $U$, (ii) a declared singlet/octet mass and mixing prescription, (iii) the leading anomaly/topological mass term if eta0 is dynamical, and (iv) the same WZW plus resonance matching convention used for the octet [@Kaiser:2000gs; @Herrera-Siklody:1996tqr]. The existing $\chi_-$ and $\chi_+$ insertions already provide the leading explicit-breaking structures in the Ruiz basis; no new breaking coefficient is forced if physical breaking is not claimed beyond that truncation [@Ruiz-Femenia:2003jdx].

Optional, and formally beyond the current leading-OZI model: independent singlet-vector couplings, all double/multi-trace odd operators, $X_{TC}$-dependent coefficient functions beyond their retained constant terms, derivative-$X_{TC}$ structures, additional $O(p^6)$ odd operators, nonideal mass mixing fitted from data, and resonance-loop/width effects. Roig--Sanz-Cillero's explicit $i\langle VV\rangle\langle\chi_-\rangle$ versus $i\langle\{V,V\}\chi_-\rangle$ pair is the recommended template for adding the first breaking/double-trace nuisance parameter [@Roig:2013baa]. A complete upgrade requires an operator enumeration followed by integration-by-parts, equation-of-motion, and $U(3)$ trace-identity reduction; it is not obtained merely by appending selected double traces.

## ALP/topological-source prescription

For the baseline, use the current $u_\mu$ and $\chi_-$ insertions and do not invent a separate $aVV$ anomaly coefficient: the $c_i,d_i$ are independent homogeneous resonance couplings [@Ruiz-Femenia:2003jdx]. With external sources and current assignment fixed, WZW plus the declared Bardeen-type counterterm convention reproduces the corresponding anomalous variation; it does not determine finite invariant local odd-parity terms or resonance matching [@Bai:2025fvl; @Fujiwara:1984mp]. $X_{TC}$-dependent coefficient functions can generate dynamical ALP interactions even without a separately targeted topological observable; omitting them is the explicit option-A large-$N_c$/operator truncation [@Alda:2026bar; @Herrera-Siklody:1996tqr; @Kaiser:2000gs].

## Recommendation and unresolved issues

Implement A first, explicitly documenting “defined leading-OZI single-trace nonet model; common singlet/octet $c_i,d_i$ assumed; independent double traces and nonconstant/derivative $X_{TC}$ structures set to zero.” This is a controlled model truncation at the accuracy of the current Ruiz basis and is preferable to adding isolated unconstrained singlet terms. Add B when singlet-vector/eta-prime observables demand subleading OZI accuracy, with one-trace and double-trace coefficients separated. Add C when the target accuracy requires anomalous-source dependence. Do not impose universal singlet couplings as a QCD theorem.

The singlet normalization and first-prediction mixing convention are fixed above. Nonideal vector mixing and new double-trace amplitudes become fit parameters only when the target accuracy or dataset requires them. The cited literature does not fix those additional singlet odd couplings.
