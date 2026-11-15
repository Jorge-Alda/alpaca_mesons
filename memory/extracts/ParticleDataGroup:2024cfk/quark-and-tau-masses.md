---
title: PDG 2024 heavy-quark, strange-quark, tau and electron masses and alpha_s(M_Z)
status: literature
created: 2026-11-15
updated: 2026-11-15
bibtex_keys: [ParticleDataGroup:2024cfk]
provenance: literature
---

# PDG 2024 masses and strong coupling

All values are read from the PDG 2024 edition (S. Navas et al., Phys. Rev. D 110, 030001; files dated 31 May 2024) [@ParticleDataGroup:2024cfk]. The downloaded PDFs and their `pdftotext -layout` text are under `memory/sources/ParticleDataGroup:2024cfk/` (ignored by Git; downloaded 2026-11-15 from `https://pdg.lbl.gov/2024/listings/` and `https://pdg.lbl.gov/2024/reviews/`).

| Quantity | Value | Listing / review (2024) | Source file | Definition and remarks |
|---|---|---|---|---|
| $m_c(m_c)$ | $1.2730\pm0.0046$ GeV | c-quark listing, "$\overline{\rm MS}$ MASS", "OUR EVALUATION" (CL = 90%) | `rpp2024-list-c-quark` | running mass at $\mu=m_c$; other schemes converted at two loops with $\alpha_s(m_c)=0.38\pm0.03$; corresponds to pole mass $1.67\pm0.07$ GeV |
| $m_b(m_b)$ | $4.183\pm0.007$ GeV | b-quark listing, "$\overline{\rm MS}$ MASS", "OUR EVALUATION" (CL = 90%) | `rpp2024-list-b-quark` | running mass at $\mu=m_b$; the listing header text quotes $4.18^{+0.04}_{-0.03}$ GeV, corresponding to pole mass $4.78\pm0.06$ GeV (two loops), with $\alpha_s(m_b)=0.223\pm0.008$ |
| $m_s$ (2 GeV) | $93.5\pm0.8$ MeV | light-quarks listing, "s-QUARK MASS", $\overline{\rm MS}$ MASS, "OUR EVALUATION" (CL = 90%) | `rpp2024-list-light-quarks` | $\overline{\rm MS}$ at $\mu=2$ GeV; the number of flavors is not stated in the listing |
| $m_t$ | $172.57\pm0.29$ GeV | t-quark listing, "t-Quark Mass (Direct Measurements)", "OUR AVERAGE" | `rpp2024-list-t-quark` | error includes scale factor 1.5; average of LHC and Tevatron kinematic reconstructions; this is the Monte Carlo generator mass, usually interpreted as the pole mass, with an interpretation uncertainty the listing says is hard to quantify |
| $m_\tau$ | $1776.93\pm0.09$ MeV | tau listing, "$\tau$ MASS", "OUR AVERAGE" | `rpp2024-list-tau` | |
| $m_e$ | $0.51099895000\pm0.00000000015$ MeV | electron listing, "e MASS" (MeV) | `rpp2024-list-electron` | 2018 CODATA value (TIESINGA 21); in $u$: $548.579909065(16)\times10^{-6}$ u |
| $\alpha_s(M_Z^2)$ | $0.1180\pm0.0009$ | review 9, "Quantum Chromodynamics", Sec. 9.4.8, Eq. (9.25) "PDG 2023 average"; revised August 2023 | `rpp2024-rev-qcd` | combination of a non-lattice average $0.1175\pm0.0010$ (Eq. 9.24) and the FLAG lattice estimate |

## Additional values read from the same pages

- Top-quark pole mass from cross-section measurements: $172.4\pm0.7$ GeV ("OUR AVERAGE", t-quark listing, "t-Quark Pole Mass from Cross-Section Measurements").
- Review 60, "Quark Masses" (revised August 2023), Eq. (60.5): lattice $m_s=(92.74\pm0.22\pm0.49[\pm0.54])$ MeV in $\overline{\rm MS}$ at 2 GeV in the $N_L=4$ theory. This differs from the listing value $93.5\pm0.8$ MeV above; the listing value is the one requested.

## Unresolved items

- Number of flavors for $\alpha_s(M_Z)$: the review text and Fig. 9.5 caption searched in `rpp2024-rev-qcd.txt` quote $\alpha_s(m_Z^2)=0.1180\pm0.0009$ without stating $n_f$. It is not listed; the standard five-flavor interpretation is an assumption, not a PDG statement.
- Number of flavors for the $m_s$ listing value: not stated in the listing.
- The review is labelled "PDG 2023 average" in Eq. (9.25) but is part of the 2024 edition files.
