# Tasks

## 1. Diagnostics without new operators

- [x] 1.1 (coder) Add a script that recomputes the on-shell $\Gamma_{K^*}$ at `basin_probes/even_shared_widths`, compares it with the tree-level RChT width at physical masses, and computes the irreducible $\chi^2$ of the two $K^*$ width rows for a common width; verify with a unittest on the irreducible-$\chi^2$ formula.
- [ ] 1.2 (coder) Add a script that rotates the whitened $\phi\to K\bar K$ residuals to $K^+K^-$ and $K^0\bar K^0$; verify that the rotation reproduces the block $\chi^2$.
- [ ] 1.3 (coder) Run refits without the $\eta'$ radiative rows and without the mass block, reporting which rows move and the odd-block and width-block $p$-values; verify optimizer status is recorded.
- [ ] 1.4 Record the results in an indexed memory file, state which candidate causes they favour, and update `memory/rcht/odd-parity/INDEX.md`; verify links.

## 2. Extension 1: separate $K^*$ charged and neutral masses and widths

- [ ] 2.1 (coder) Add options for separate $K^{*+}$ and $K^{*0}$ masses and widths and for the revised $K\pi$ treatment, defaulting to baseline; verify bitwise-equal residuals at defaults.
- [ ] 2.2 Refit from the accepted endpoints, re-evaluate at $(4,128,128)$, apply the adequacy criterion; store the ledger entry.

## 3. Extension 2: SU(3)-breaking VPP coupling

- [ ] 3.1 Propose and get approval for the notation of the $K^*K\pi$ to $\rho\pi\pi$ ratio; derive its effect on $\Gamma_{K^*}$, $\Gamma_\phi$ and $\phi\to K\bar K$ with a SymPy script and a memory record.
- [ ] 3.2 (coder) Implement the option (default baseline), refit, apply the criterion, store the ledger entry. If the criterion is met, skip to section 7.

## 4. Extension 3: double-trace source term

- [ ] 4.1 Derive whether $y\,\mathrm{tr}(\chi_+Q)\mathbb{1}$ enters the vector mass spectrum and the $KK\pi$ current (SymPy script, memory record); record "not applicable" if it does not.
- [ ] 4.2 (coder) If applicable, implement the $y$ coordinate (notation approved first), refit, apply the criterion, store the entry.

## 5. Extension 4: isospin breaking in the even sector

- [ ] 5.1 Derive $\rho$–$\omega$ mixing and the $\phi\to K^+K^-$ versus $K^0\bar K^0$ difference with a SymPy script and a memory record, using the channel diagnostic of 1.2 to choose the minimal parameterization.
- [ ] 5.2 (coder) Implement, refit, apply the criterion, store the entry.

## 6. Extension 5: $\eta'$ and odd-block tension

- [ ] 6.1 Using the $\eta'$-exclusion result of 1.3, decide and record what $\eta'$ treatment is justified; implement it as an option if so, refit, apply the criterion.

## 7. Outcome

- [ ] 7.1 Write the final memory record: adequate extended model with its accepted endpoint, or the documented reason the criterion cannot be met, with every extension and combination tried; link it from the indexes.
- [ ] 7.2 Update `bd` for alp-uzh.10 and alp-uzh.7 according to the outcome; propose a commit message.
