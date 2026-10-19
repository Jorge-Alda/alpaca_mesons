# Instructions

This folder will contain a research project in the field of hep-ph. Your task is to create all the structure necessary for an AI research assistant based on a codex harness.

## Priorities

### 1. Scientific consistency

- Every source used must be included in the `bibliography.bib` file, using the bibtex information from inspire-hep
- Every concept, definition, equation, or piece of information learnt from a bibliographical source must be documented as a markdown file in the `/memory` folder, with a reference to the bibtex key
- Analytical derivations must be reproducible, using two files: a python script, in the `/scripts` folder, containing the derivation in `sympy`, and a markdown file in `/memory` detailing step-by-step solutions
- Numerical derivations must be reproducible, using python scripts in `/scripts` and explanations in `/memory`
- Derivations, both analytical and numerical, must indicate the provenance of every input by linking to the corresponding files in `/memory`
- Any approximation or assumption must be cleared stated as such
- All files in `/memory` must include metadata regarding their evidenciary status: obtained from literature, derived, conjectured,...

### 2. Context management and token efficiency

- Keep all files, including the AI scaffolding, small and modular: information and instructions must be loaded only as needed
- Save all new information to memory as markdown files, and keep to a minimum direct output to the terminal: only list of modified files, and follow-up questions to the user, if any
- Use a concise and precise language. Do not waste tokens in embellishments, connectors or unnecessary fluff
- **Never ever** praise the user "Good question!"

## Bibliographical management

- Store the bibliographical sources in `/memory/sources`, indexed by their bibtex key. This folder is not synched to the remote repository, so you might need to re-download sources between sessions or users.
- If a paper is in arXiv, fetch the source files instead of the rendered pdf.
- The results of digesting a paper will be saved in `/memory/extracts`, with a subdirectory for each paper corresponding to their bibtex key. Inside the folder, multiple extracts might exist, according to the modularity guidelines.

## Subagents

Create files for the following subagents, with specific instructions and configurations for codex:

- Orchestrator: GPT 5.6-Sol Medium, with default temperature
- Brainstormer: GPT 5.6-Sol Medium, higher temperature
- Coder: GPT 5.6-Sol Medium, low temperature
- Bibliography manager: GPT 5.6-Luna High, low temperature

## Skills

Create a skill file for bibliography management, that will be used by the Bibliography manager subagent.