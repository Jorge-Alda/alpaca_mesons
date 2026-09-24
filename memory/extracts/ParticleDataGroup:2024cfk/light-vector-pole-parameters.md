---
title: PDG 2024 light-vector masses and total widths for Breit-Wigner inputs
status: literature
created: 2026-09-24
updated: 2026-09-24
bibtex_keys: [ParticleDataGroup:2024cfk]
provenance: literature
---

# Light-vector mass and width inputs

The following are central values and quoted one-standard-deviation uncertainties transcribed from the 2024 PDG particle listings. All values are in MeV; ``mass'' and ``width'' mean the listing's fitted quantities, not S-matrix pole parameters. [@ParticleDataGroup:2024cfk]

| state and listing category | mass | total width | source location |
|---|---:|---:|---|
| $\rho(770)^0$, neutral only, $e^+e^-$ | $775.26\pm0.23$ | $147.4\pm0.8$ | rho listing pp. 1, 7 |
| $\rho(770)^\pm$, charged only, $\tau$ decays and $e^+e^-$ | $775.11\pm0.34$ | $149.1\pm0.8$ | rho listing pp. 2, 8 |
| $K^*(892)^\pm$, charged only, hadroproduced | $891.67\pm0.26$ | $51.4\pm0.8$ | K-star listing pp. 1, 5 |
| $K^*(892)^0$, neutral only | $895.55\pm0.20$ | $47.3\pm0.5$ | K-star listing pp. 2--3, 6 |
| $\omega(782)$ | $782.66\pm0.13$ | $8.68\pm0.13$ | omega listing pp. 1--2 |
| $\phi(1020)$ | $1019.461\pm0.016$ | $4.249\pm0.013$ | phi listing pp. 1--2 |

The rho charged and neutral entries are kept as separate PDG category averages. The K-star charged mass/width entries likewise use the charged-only hadroproduction averages, while the neutral entries use the neutral-only averages. These category choices are literature transcription choices, not an isospin-equality assumption. [@ParticleDataGroup:2024cfk]

## Modeling convention

Using these literature values as $m_V$ and $\Gamma_V$ in a constant-width Breit--Wigner denominator, for example $D_V(s)=m_V^2-s-i m_V\Gamma_V$, is a modeling assumption. It is not a claim that the PDG listing itself endorses a constant-width line shape: the listings contain fit- and channel-dependent parameterizations, and the rho listing explicitly labels the displayed quantities as fit averages. The table therefore records experimental inputs separately from the subsequent constant-width approximation.

## Auditable source material and URLs

Cached PDFs and text conversions are under `memory/sources/ParticleDataGroup:2024cfk/`:

- `rpp2024-list-rho-770.pdf` / `.txt`, https://pdg.lbl.gov/2024/listings/rpp2024-list-rho-770.pdf (rho mass p. 1--2; width p. 7--8).
- `rpp2024-list-K-star-892.pdf` / `.txt`, https://pdg.lbl.gov/2024/listings/rpp2024-list-K-star-892.pdf (mass p. 1--3; width p. 5--6).
- `rpp2024-list-omega-782.pdf` / `.txt`, https://pdg.lbl.gov/2024/listings/rpp2024-list-omega-782.pdf (mass p. 1; width p. 2).
- `rpp2024-list-phi-1020.pdf` / `.txt`, https://pdg.lbl.gov/2024/listings/rpp2024-list-phi-1020.pdf (mass p. 1; width p. 2).
- Listing index: https://pdg.lbl.gov/2024/html/contents_listings.html.

The PDG listings are dated 2024-05-31 and cite S. Navas et al. (Particle Data Group), *Phys. Rev. D* 110 (2024) 030001, keyed in this repository as `ParticleDataGroup:2024cfk`. [@ParticleDataGroup:2024cfk]
