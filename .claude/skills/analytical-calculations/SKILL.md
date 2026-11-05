---
name: analytical-calculations
description: Derive analytical results reproducibly with SymPy scripts, step-by-step LaTeX output, and project-consistent notation.
---

# Analytical Calculations

Use this skill for every analytical derivation, symbolic manipulation, series expansion, matrix calculation, equation solution, or algebraic simplification.

## Mandatory Environment

1. Use the repository's existing `.venv` for every Python operation.
2. Run scripts with `.venv/bin/python`; never use the system `python`, `python3`, or another environment.
3. Install a missing dependency only into `.venv`, using `uv add <package>` followed by `uv pip install -r pyproject.toml`.
4. Do not recreate, replace, or modify the environment configuration unless explicitly requested.
5. If `.venv` is absent or unusable, stop and ask the user how to proceed. Do not silently fall back to another interpreter.

## SymPy Requirement

1. Perform every analytical calculation in an executable Python script under `scripts/` using SymPy.
2. Use SymPy's standard objects, transformations, solvers, matrix expressions, tensor tools, quantum operators, series methods, simplifiers, printers, and other built-in functionality whenever applicable.
3. Search SymPy's documented API before writing new symbolic manipulation logic.
4. Create custom logic only when the required operation does not exist in SymPy. Keep such logic minimal and document why no SymPy operation is suitable.
5. Do not implement replacement symbolic classes, expression trees, commutator algebras, polynomial systems, matrix systems, simplifiers, or renderers when SymPy already provides them.

## Derivation Workflow

1. Read the relevant `memory/*.md` and `memory/extracts/**/*.md` files before writing the script.
2. Reuse every established variable, field, coefficient, operator, and convention name exactly as it appears in memory. Do not invent aliases or replacement names for quantities that already have names.
3. Encode the actual step-by-step derivation in SymPy. The script must be the calculation, not a transcription of a hand-derived result followed by a numerical or symbolic check.
4. Let each result be produced from earlier SymPy expressions through explicit operations such as expansion, substitution, differentiation, integration, coefficient extraction, solving, simplification, or series expansion.
5. Execute the script with `.venv/bin/python scripts/<script>.py`.
6. Take the partial and final analytical results from the executed script's output. Do not substitute independently derived expressions for the script output.
7. Record the derivation and invocation in the corresponding `memory/` file according to the repository research protocol.

## Script Output

1. Print the important input expressions, intermediate steps, substitutions, partial results, and final result.
2. Render every displayed analytical expression with SymPy's standard LaTeX renderer, `sympy.latex`.
3. Give each printed expression a concise textual label so the derivation order is unambiguous.
4. Keep diagnostic assertions separate from the derivation output. Assertions may verify the result, but they must not replace the step-by-step calculation.
5. Make output deterministic and state assumptions, truncation orders, conventions, and units where relevant.

A minimal output helper is:

```python
import sympy as sp


def show(label: str, expression: sp.Expr) -> None:
    print(f"{label}:\n{sp.latex(expression)}\n")
```

## Matrices And Noncommutative Algebra

1. Represent abstract matrices with `sympy.MatrixSymbol` and abstract noncommuting quantities with `sympy.Symbol(..., commutative=False)` unless the user explicitly requests component matrices or supplies explicit entries.
2. Use SymPy's standard matrix-expression and noncommutative operations, including `MatAdd`, `MatMul`, `MatPow`, `Trace`, `Transpose`, `Adjoint`, `Inverse`, `Commutator`, and `AntiCommutator`, as appropriate.
3. Do not replace an abstract matrix calculation with generic component matrices merely to make simplification easier unless the user explicitly asks for components.
4. Do not assume that matrices, operators, fields, or couplings commute unless that assumption is stated in memory or explicitly requested.
5. State every additional matrix property used, including dimensions, Hermiticity, diagonality, symmetry, invertibility, or trace cyclicity.

## Completion Checks

- The script runs successfully with `.venv/bin/python`.
- SymPy performs the derivation itself and produces all reported expressions.
- Partial and final expressions are printed through `sympy.latex`.
- Existing memory notation is preserved exactly.
- Any custom symbolic logic is justified by a missing SymPy capability.
- The memory derivation cites every literature input and links the executable script.
