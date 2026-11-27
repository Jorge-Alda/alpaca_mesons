# Proposal

## Why

The change `estimate-sigma-pipi-one-loop-murayama` sizes the one-loop correction to $\sigma\to\pi^+\pi^-$ in the generalized Murayama model with the polar fields $(\sigma,a)$ of $\varphi=(f_a+\sigma)e^{ia/f_a}$. The tree amplitude does not depend on this choice (`radial-decay-polar-cartesian-general-sources.md`), but a size estimate that replaces each scalar moment by $(m^2)^{2-r}/(16\pi^2)$ and each loop-momentum product by $m^2$ gives only an order of magnitude for each graph class, and can miss cancellations between graphs with the same $\epsilon$ scaling. Such cancellations can differ between parametrizations. In the Cartesian fields $\varphi=f_a+s+ip$ the kinetic term is free and all interactions sit in the potential, so the same procedure applied there may keep or miss different cancellations. A second, independent estimate is a diagnostic of this limitation.

## What Changes

- Carry out the estimate of the one-loop contribution to $\sigma\to\pi^+\pi^-$ **from scratch in the Cartesian fields** $(s,p)$. The source is expanded directly in the Cartesian fields, $\widehat M(s,p)=\sum_{r,s}(\widehat S^C_{r,s}+i\widehat P^C_{r,s})(s/f_a)^r(p/f_a)^s$, and the estimate is written in the invariants $\widehat S^C_{r,s}$, $\widehat P^C_{r,s}$. The Feynman rules, the graph list and the traces are derived from the Cartesian Lagrangian, not translated from the polar records or through the field map.
- Keep the assumptions and approximations of the polar estimate: size rule, loop scales (D3), routing rule (D12), exclusion of the meson-only chiral background (D4), CP and flavour-diagonality (D11), large-$\epsilon$ regime (D10). Each is restated in the Cartesian invariants (for example no tadpoles becomes $\widehat P^C_{0,0}=0$, $\langle\widehat S^C_{0,1}\rangle=0$).
- Derive the vertices up to five legs, graphs, leg factors, generator-sum treatment, external residues, traces, $M_1^{\rm est}$, $R=|M_1^{\rm est}|/|M_0|$, large-$\epsilon$ powers and the up-only limit.
- The only comparison with the polar results used in the derivation is the equivalence of the tree amplitude $\sigma\to\pi^+\pi^-$ (and consistency with the existing notes). The final comparison of the two estimates is a diagnostic. **No agreement is expected**, including with the polar up-only prefactor $-21/(64\pi^2)$, and neither estimate is authoritative.
- **No loop integral, scalar loop function, loop coefficient or loop amplitude is computed.** All tests are exact symbolic identities, with compact trace notation only (no component notation, no $E_{ij}$, no random or numerical matrices).

## Capabilities

### New Capabilities

- `radial-mode-pipi-loop-estimate-cartesian`: independent analytical size estimate of the one-loop $\sigma\to\pi^+\pi^-$ amplitude in the Cartesian fields, in the invariants $\widehat S^C_{r,s}$, $\widehat P^C_{r,s}$, under the assumptions of the polar estimate.

### Modified Capabilities

None. The polar requirements belong to `radial-mode-pipi-loop-estimate`.

## Impact

- New records under `memory/alp/radial-mode/finite-mixing/two-pion-loop/` (a Cartesian subfolder with its own index) and indexes updated.
- New symbolic scripts under `scripts/`, written by the coder agent only; generic trace-algebra tools may be reused, the vertices and graphs are not.
- No change to the polar records and scripts, to `latex/`, or to the main specs. No new bibliography entries expected (Alda:2026nrd, DiLuzio:2026poh).
