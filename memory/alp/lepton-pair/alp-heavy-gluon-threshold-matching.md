---
title: Leading heavy-flavor matching of the ALP gluonic coupling
status: derived
created: 2026-11-15
updated: 2026-11-15
bibtex_keys: [Alda:2026nrd, Messerli:2025rnv]
provenance: analytical-derivation
---

# Scope and inputs

Use the conventions and natural units of the [high-energy conventions](alp-high-energy-tff-conventions.md). The input is the invariant $c_{GP}=c_G+\sum_f c_{ff}^P$ and the heavy-flavor triangle's zero-momentum limit, as recorded in the [triangle calculation](alp-heavy-quark-triangle.md) following Alda:2026nrd. The matching is at leading order in $\alpha_s$ and $1/f_a$, for external gluon momenta and $m_a$ small compared with the heavy-quark mass $m_h$. It does not establish local matching in the charmonium region.

# Derivation

The heavy-quark triangle's pseudoscalar term approaches the constant proportional to $c_{hh}^P$, while its derivative term vanishes, as $m_h$ tends to infinity. Equivalently, matching the invariant $c_{GP}$ above and below the threshold gives
$$
c_G^{(n_f-1)}+\sum_{f\ne h}c_{ff}^P=c_G^{(n_f)}+\sum_f c_{ff}^P,
\qquad c_G^{(n_f-1)}=c_G^{(n_f)}+c_{hh}^P.
$$
For $c_{GA}=c_G+\sum_f c_{ff}^A$ this implies
$$
c_{GA}^{(n_f-1)}=c_{GA}^{(n_f)}-c_{hh}^A+c_{hh}^P.
$$
This is a source matching relation, not a derivation of the renormalized singlet-current matching coefficient. Its corrections include powers of external momenta over $m_h$ and QCD corrections. When $m_a$ approaches or exceeds $2m_c$, local charm matching is not justified by this expansion.

# Verification and remaining limitations

The [script](../../../scripts/derive_alp_symmetric_tff.py) solves the matching equation, checks invariance under the stated heavy-flavor rotation, and checks the heavy triangle's limits and the two Ward contractions and Bose symmetry of the leading tensor. Run:

```text
.venv/bin/python scripts/derive_alp_symmetric_tff.py
.venv/bin/python scripts/derive_alp_symmetric_tff.py --test
```

The script also prints a conditional sum of the contact term, heavy triangle, and light axial-current term from the [hard-current OPE](alp-hard-current-ope.md). This sum is not established as the complete coefficient: the overlap between the explicit hard light-quark contribution and the renormalized axial matrix element must be fixed, and the [singlet-current record](alp-singlet-axial-current-evolution.md) leaves the gluonic Wilson coefficient and threshold matching open. The [heavy-flavor record](alp-heavy-flavor-thresholds.md) states the threshold-continuity assumption. No numerical uncertainty is assigned to these omissions.

# Perturbative overlap and the symmetric local OPE

The [triangle calculation](alp-heavy-quark-triangle.md) also provides a direct test of how the hard triangle and the local axial-current term overlap. In its Euclidean momentum-cutoff convention, retaining the leading $1/|q^2|$ term and the first term in $m_a^2/m_h^2$, the script subtracts the axial-current term from the full triangle's mass-dependent term and obtains
$$
-\frac{N_c Q_h^2m_h^2(c_{hh}^A-c_{hh}^P)}{\pi^2|q^2|f_a}
\left[\log\frac{|q^2|}{\Lambda^2}+1\right].
$$
This is generally nonzero at the same power in $1/|q^2|$ as the local axial-current term. It is a hard remainder in the stated cutoff calculation, not an $\overline{\rm MS}$ matching coefficient: both the finite constant and the separation depend on the prescription. Adding a complete free triangle and its own perturbative axial-current matrix element would count the latter twice; conversely, adding only the axial-current term omits this remainder. The comparison is limited to $m_a^2\ll m_h^2$ and does not by itself supply the corresponding light-flavor coefficient or the matching of the ChPT current to the renormalized QCD current.

For the symmetric local expansion in massless QCD, Ludtke:2020moa gives the axial current as the only dimension-three operator through $O(\alpha_s)$, with Wilson coefficient $1-\alpha_s/\pi$; see the [extract](../../extracts/Ludtke:2020moa/ope-nlo-wilson-coefficient.md). This statement does not calculate the mass-dependent source remainder above, nor does it determine a finite singlet-current threshold matching constant. At strictly leading order in $\alpha_s$, tree-level current matching is continuous; corrections and the renormalized identification of the low-energy matrix element remain to be specified. Thus the printed conditional sum is not promoted to a complete source-linear result.

The script also derives the same remainder without expanding in $m_a^2/m_h^2$, for $0\le m_a^2<4m_h^2$. In the triangle Feynman parameters put $z$ on the line between the two photon vertices and $x$ on the remaining edge. The denominator is
$$
m_h^2-x(1-x)(1-z)^2m_a^2+z(1-z)|q^2|.
$$
Split the $z$ integral at $\delta$ with $m_h^2/|q^2|\ll\delta\ll1$. The small-$z$ integral supplies $\log[|q^2|\delta/(m_h^2-x(1-x)m_a^2)]/|q^2|$; the complementary integral supplies $-\log\delta/|q^2|$. Their sum reproduces the full leading logarithm. Subtracting the cutoff axial-current matrix element, whose integrand is $[\log(\Lambda^2/(m_h^2-x(1-x)m_a^2))-1]/|q^2|$, gives the same $[\log(|q^2|/\Lambda^2)+1]/|q^2|$ for every $x$. Thus the displayed remainder holds for finite $m_a$ below threshold at leading $1/|q^2|$, within the same cutoff prescription. The result does not determine how a nonperturbative light-quark matrix element should be matched in $\overline{\rm MS}$.

## One-loop $\overline{\rm MS}$ subtraction

For the free-quark comparison, use dimensional regularization in $d=4-2\epsilon$ and subtract the $\overline{\rm MS}$ pole of the scalar two-point integral. With $D=m_h^2-x(1-x)m_a^2>0$, the dimensionless factor after extracting $i/(16\pi^2)$ is
$$
e^{\gamma_E\epsilon}\Gamma(\epsilon)\left(\frac{\mu^2}{D}\right)^\epsilon
=\frac1\epsilon+\log\frac{\mu^2}{D}+O(\epsilon).
$$
The script obtains this expansion and subtracts the finite two-point contribution from the triangle's leading hard logarithm. In the one-loop current prescription that subtracts this pole without an additional finite source-current counterterm, the resulting hard remainder is
$$
-\frac{N_cQ_h^2m_h^2(c_{hh}^A-c_{hh}^P)}{\pi^2|q^2|f_a}
\log\frac{|q^2|}{\mu^2}.
$$
This reproduces the cutoff result under $\Lambda^2=e\mu^2$ at this order. The comparison uses anticommuting $\gamma_5$ in the two-point graph, where two $\gamma_5$ matrices occur, and the triangle normalization already fixed by the Pauli--Villars calculation in the [triangle record](alp-heavy-quark-triangle.md). It is a perturbative one-loop matching calculation, not a nonperturbative determination of the $\overline{\rm MS}$ current or of finite source-current contact terms. Changing the finite source-current prescription shifts the current term and the hard remainder oppositely. The singlet finite renormalization and QCD corrections are beyond this free-quark calculation. As with the cutoff result, this expression is obtained at leading $1/|q^2|$ for $m_a^2<4m_h^2$ and is not by itself a charm-threshold hadronic calculation.
