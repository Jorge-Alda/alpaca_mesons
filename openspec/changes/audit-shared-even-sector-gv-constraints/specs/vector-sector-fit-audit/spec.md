# Spec Delta

## Purpose

This capability produces reproducible diagnoses of vector-sector likelihood blocks whose even-sector and conditional decay constraints share measurements or model inputs. It makes tensions and retained assumptions auditable before downstream joint-fit minima are interpreted.

## ADDED Requirements

### Requirement: Audit inputs and shared information

The audit SHALL identify the provenance, observable, uncertainty treatment, and model role of every included even-sector and conditional $G_V$ input, and SHALL identify overlaps that prevent likelihood contributions from being treated as independent.

#### Scenario: Shared width or branching input is reused

- **WHEN** a measurement summary contributes to both a conditional diagnostic and a joint-fit constraint
- **THEN** the audit identifies the overlap and reports the conditional diagnostic separately without adding it as independent likelihood information

#### Scenario: Input covariance is incomplete

- **WHEN** source-level cross-correlations are unavailable
- **THEN** the audit states which covariance entries are retained and which are approximated as zero

### Requirement: Quantify blockwise tension

The audit SHALL report the stored joint-fit block contributions and a reproducible decomposition of the shared even-sector and conditional $G_V$ tension, including the change in likelihood for every justified comparison it evaluates.

#### Scenario: Baseline joint fit is not converged

- **WHEN** the stored fit reports optimizer failure or an evaluation limit
- **THEN** the audit labels all comparisons as diagnostics against that stored snapshot and does not describe it as a converged minimum

#### Scenario: A physically justified alternative is evaluated

- **WHEN** an alternative mass, width, saturation, or covariance prescription is supported by the available inputs
- **THEN** the audit records its assumptions, compared likelihood blocks, and numerical effect without conflating it with the baseline result

### Requirement: Preserve mass and width conventions

The audit SHALL distinguish physical masses used in external phase space from fitted tensor masses and any averaged masses used in vertices or propagators, and SHALL state the charged/neutral prescriptions for shared widths and $K^*\to K\pi$ saturation.

#### Scenario: Charged and neutral decay rows are combined

- **WHEN** charged and neutral $K^*$ observables enter a shared model constraint
- **THEN** the audit reports the row-level mass inputs, combination rule, and saturation assumption used to predict them
