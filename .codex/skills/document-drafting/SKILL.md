---
name: document-drafting
description: Draft HEP-ph research documents based on memory excerpts from other agents.
---

# Document Drafting Skill

## Inputs

- Memory excerpts from other agents containing relevant HEP-ph knowledge. They are hierarchically organized, with higher-level summaries providing context for more detailed records. Load to the context only the selected memory records.
- Specific instructions or prompts for drafting the document.
- Any relevant bibliography entries or references that should be cited in the draft. Copy the bibtex entries from the repository-level bibliography files as needed.

## Outputs

- A first draft of a coherent and well-structured HEP-ph research document, in LaTeX format.
- Properly cited references and bibliography entries.
- Clear distinction between established knowledge and conjecture.
- Supplementary material or appendices for implementation details, if applicable.

## Guidelines

- Ensure clarity, coherence, and adherence to HEP-ph conventions.
- Maintain a clear distinction between established knowledge and conjecture.
- Highlight any assumptions made in the reasoning process.
- Separate theoretical discussions from practical implementation details.
- Include supplementary material or appendices for implementation details, if necessary.
- Ensure the main text is self-contained and comprehensible without requiring the reader to refer to supplementary material.
- Keep minimal prose. Focus on conveying complex ideas clearly and concisely. Do not use unnecessary words or overly complex sentence structures. Do not use connectors between sentences. Use `itemize` profusely.
- Use consistent notation and terminology throughout the document.
- Ensure that all mathematical expressions and equations are clearly defined and properly formatted in LaTeX.
- Use LaTeX for internal references and citations within the document. For bibliography management, use BibTeX and ensure that all cited references are included in the bibliography section. Never introduce a reference by the name of the authors, always use the `\cite{}` command.
- The structure of the LaTeX project must be modular, with separate files for the sections and appendices. Ensure that the main file includes these modular files correctly using `\input{}` or `\include{}` commands.