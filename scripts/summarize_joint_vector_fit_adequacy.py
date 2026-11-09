"""Markdown summary of the goodness-of-fit and sensitivity ledgers of the joint vector fit.

Run: .venv/bin/python scripts/summarize_joint_vector_fit_adequacy.py

Reads (read only, no fit is run):
  data/vector_sector/joint_vector_fit_adequacy.json          (four accepted minima at (4,128,128))
  data/vector_sector/joint_vector_block_exclusion.json       (one-block-removed refits)
  data/vector_sector/joint_vector_sensitivity_variations.json (predefined variations)
and prints Markdown tables to standard output. The output is deterministic.

Conventions: chi^2 is the sum of squared standardized residuals (pulls);
p-values are the upper tail of the chi^2 distribution, quoted in scientific
notation. Per-block p-values use the block row count as degrees of freedom
(assumption stated in the ledgers; parameters are not subtracted, so these
p-values only rank blocks). Parameter shifts are in units of the baseline
(4,128,128) local errors of the best accepted endpoint. A refit is flagged
"max_nfev" when the optimizer stopped at the 5000 function-evaluation limit.
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
ADEQUACY = 'data/vector_sector/joint_vector_fit_adequacy.json'
EXCLUSION = 'data/vector_sector/joint_vector_block_exclusion.json'
VARIATIONS = 'data/vector_sector/joint_vector_sensitivity_variations.json'
MAXIMUM_FUNCTION_EVALUATIONS = 5000
WIDTH_AND_PHI_KK_BLOCKS = ('shared_widths', 'phi_KK')


def format_p_value(value):
    """p-value in scientific notation; None (non-positive dof) is shown as n/a."""
    return 'n/a' if value is None else f'{value:.3e}'


def table(header, rows):
    """Markdown table from a header list and a list of row lists."""
    lines = ['| ' + ' | '.join(header) + ' |', '|' + '---|'*len(header)]
    lines += ['| ' + ' | '.join(str(cell) for cell in row) + ' |' for row in rows]
    return '\n'.join(lines)


def hit_function_evaluation_limit(refit_record):
    """True when the optimizer stopped at the function-evaluation limit."""
    evaluations = refit_record.get('cumulative_nfev', refit_record.get('nfev', 0))
    return (evaluations >= MAXIMUM_FUNCTION_EVALUATIONS
            or 'maximum number of function evaluations' in refit_record.get('message', ''))


def worst_block(blocks):
    """Block with the smallest row-count p-value."""
    name = min(blocks, key=lambda block: blocks[block]['nominal_p_value_row_count_dof'])
    return name, blocks[name]


def degrees_of_freedom_text(evaluation):
    """Default dof, with the count that ignores zero Jacobian columns when different."""
    default = evaluation['degrees_of_freedom']
    alternative = evaluation.get('degrees_of_freedom_without_zero_columns', default)
    return str(default) if alternative == default else f'{default} ({alternative} w/o zero columns)'


def p_value_text(evaluation):
    """Default p-value, with the zero-column alternative when the dof differ."""
    text = format_p_value(evaluation['nominal_p_value'])
    if evaluation.get('degrees_of_freedom_without_zero_columns', evaluation['degrees_of_freedom']) \
            != evaluation['degrees_of_freedom']:
        text += f" ({format_p_value(evaluation['nominal_p_value_without_zero_columns'])})"
    return text


def block_table(adequacy):
    """Per-block chi^2, rows and nominal p for each accepted minimum, plus totals."""
    endpoints = adequacy['endpoints']
    header = ['block', 'rows'] + [f'{e["id"]} chi2 / p' for e in endpoints]
    rows = []
    for name in endpoints[0]['blocks']:
        cells = [f"{e['blocks'][name]['chi2']:.4f} / "
                 f"{format_p_value(e['blocks'][name]['nominal_p_value_row_count_dof'])}" for e in endpoints]
        rows.append([name, endpoints[0]['blocks'][name]['count']] + cells)
    rows.append(['**total**', endpoints[0]['residual_count']]
                + [f"{e['chi2']:.6f} / {format_p_value(e['nominal_p_value'])}" for e in endpoints])
    rows.append(['dof', ''] + [e['degrees_of_freedom'] for e in endpoints])
    rows.append(['chi2/dof', ''] + [f"{e['chi2_per_degree_of_freedom']:.4f}" for e in endpoints])
    return table(header, rows)


def largest_pull_table(adequacy):
    """Fifteen largest pulls at the best minimum, with the same rows at the other minima."""
    endpoints = adequacy['endpoints']
    pulls = [{row['index']: row['pull'] for row in e['standardized_residuals']} for e in endpoints]
    header = ['rank', 'row', 'block', 'chi2 contribution'] + [f'pull {e["id"]}' for e in endpoints]
    rows = [[item['rank'], item['label'], item['block'], f"{item['chi2_contribution']:.4f}"]
            + [f'{pull[item["index"]]:+.4f}' for pull in pulls] for item in endpoints[0]['largest_pulls']]
    best_set = {item['index'] for item in endpoints[0]['largest_pulls']}
    differing = [e['id'] for e in endpoints[1:] if {item['index'] for item in e['largest_pulls']} != best_set]
    note = ('Same top-15 row set at all four minima.' if not differing
            else f'Top-15 row set differs at: {", ".join(differing)}.')
    return table(header, rows) + '\n\n' + note


def sign_pattern_by_range(pattern):
    """Residual sign sub-sequence of each energy subrange (bins assigned by centre)."""
    signs = []
    for item in pattern['energy_ranges']:
        low, high = item['energy_low_gev'], item['energy_high_gev']
        signs.append(''.join('+' if pull > 0 else '-'
                             for centre, pull in zip(pattern['energy_center_gev'], pattern['pulls'])
                             if low <= centre < high))
    return signs


def babar_table(adequacy):
    """BaBar chi^2 and residual sign pattern per channel and energy subrange."""
    ranges = adequacy['energy_ranges_gev']
    header = (['minimum', 'channel', 'bins', 'chi2', 'runs (expected, z)']
              + [f'{low:.2f}-{high:.2f} GeV: chi2 (n) [+/-] signs' for low, high in ranges])
    rows = []
    for endpoint in adequacy['endpoints']:
        for channel, pattern in endpoint['babar_energy_pattern'].items():
            runs = pattern['runs']
            cells = [f"{item['chi2']:.3f} ({item['count']}) [{item['positive']}/{item['negative']}] {signs}"
                     for item, signs in zip(pattern['energy_ranges'], sign_pattern_by_range(pattern))]
            rows.append([endpoint['id'], channel, pattern['count'], f"{pattern['chi2']:.4f}",
                         f"{runs['runs']} ({runs['expected_runs']:.2f}, {runs['standardized_runs']:+.2f})"]
                        + cells)
    full_sequences = '\n'.join(f"- {e['id']} {channel}: `{pattern['sign_sequence']}`"
                               for e in adequacy['endpoints'][:1]
                               for channel, pattern in e['babar_energy_pattern'].items())
    return table(header, rows) + '\n\nFull sign sequences (increasing energy) at the best minimum:\n' \
        + full_sequences


def exclusion_table(exclusion):
    """One-block-removed refits: optimizer status and chi^2 of remaining and removed rows."""
    header = ['excluded block', 'optimizer status', 'remaining chi2', 'dof', 'p',
              'excluded-block chi2 (rows)', 'baseline chi2 of that block', 'max |shift| (parameter)',
              'conditional']
    rows = []
    for entry in exclusion['entries']:
        refit, evaluation = entry['refit'], entry['high_order']
        status = ('success' if refit['success'] else 'no success') + f" (status {refit['status']}"
        status += ', max_nfev)' if hit_function_evaluation_limit(refit) else ')'
        name = entry['excluded_blocks'][0]
        rows.append([name, status, f"{evaluation['chi2']:.4f}", degrees_of_freedom_text(evaluation),
                     p_value_text(evaluation),
                     f"{evaluation['removed_rows_chi2']:.4f} ({len(evaluation['removed_row_pulls'])})",
                     f"{entry['baseline_blocks'][name]['chi2']:.4f}",
                     f"{entry['maximum_absolute_shift']:.3f} ({entry['maximum_shift_parameter']})",
                     entry['conditional']])
    return table(header, rows)


def basin_rerun_text(entry):
    """Outcome of the reruns from the other accepted endpoints."""
    if not entry['basin_switch_reruns']:
        return 'not triggered'
    runs = entry['runs']
    accepted = sum(run['accepted'] for run in runs)
    return (f"{len(runs)} starts, {accepted} accepted; selected run {entry['selected_run']} "
            f"({runs[entry['selected_run']]['start_endpoint']})")


def variation_table(variations):
    """Predefined variations: acceptance, chi^2, worst block, largest shift, reruns."""
    header = ['variation', 'outcome', 'chi2', 'dof', 'p', 'worst block: chi2 (rows) p',
              'max |shift| (parameter)', 'basin rerun']
    rows = []
    for entry in variations['entries']:
        run = entry['runs'][entry['selected_run']]
        evaluation = run['high_order']
        if run['accepted']:
            outcome = 'accepted'
        else:
            reason = 'max_nfev' if hit_function_evaluation_limit(run) else \
                run['quadrature_refinement'].get('status', 'audit failed')
            outcome = f'unaccepted ({reason})'
        limit_runs = sum(hit_function_evaluation_limit(item) for item in entry['runs'])
        if limit_runs and len(entry['runs']) > 1:
            outcome += f'; {limit_runs}/{len(entry["runs"])} runs at max_nfev'
        name, block = worst_block(evaluation['blocks'])
        rows.append([entry['name'], outcome, f"{evaluation['chi2']:.4f}", degrees_of_freedom_text(evaluation),
                     p_value_text(evaluation),
                     f"{name}: {block['chi2']:.3f} ({block['count']}) "
                     f"{format_p_value(block['nominal_p_value_row_count_dof'])}",
                     f"{run['maximum_absolute_shift']:.3f} ({run['maximum_shift_parameter']})",
                     basin_rerun_text(entry)])
    deferred = '\n'.join(f"- deferred: {item['name']} ({item['reason']})"
                         for item in variations.get('not_implemented', []))
    accepted = sum(entry['runs'][entry['selected_run']]['accepted'] for entry in variations['entries'])
    counts = (f"{len(variations['entries'])} variations run, {accepted} with an accepted selected run, "
              f"{len(variations.get('not_implemented', []))} deferred.")
    return table(header, rows) + '\n\n' + counts + '\n' + deferred


def width_and_phi_table(adequacy):
    """chi^2 contribution of each shared-width row and of the phi -> K Kbar rows."""
    endpoints = adequacy['endpoints']
    pulls = [{row['index']: row for row in e['standardized_residuals']} for e in endpoints]
    header = ['row', 'block'] + [f'chi2 (pull) {e["id"]}' for e in endpoints]
    rows = []
    for index, label in enumerate(adequacy['row_labels']):
        if label['block'] not in WIDTH_AND_PHI_KK_BLOCKS:
            continue
        name = label['label'] + (' [Gamma_etap]' if label['label'] == 'total width etap' else '')
        rows.append([name, label['block']]
                    + [f"{pull[index]['pull']**2:.4f} ({pull[index]['pull']:+.4f})" for pull in pulls])
    return table(header, rows)


def summary(adequacy, exclusion, variations):
    """Full Markdown summary."""
    sections = [
        ('Per-block chi2 at the four accepted minima (4,128,128)', block_table(adequacy)),
        ('Fifteen largest pulls (ranked at the best minimum)', largest_pull_table(adequacy)),
        ('BaBar residual pattern by channel and energy subrange', babar_table(adequacy)),
        ('Block exclusions (refit from the best minimum, evaluated at (4,128,128))', exclusion_table(exclusion)),
        ('Predefined variations (selected run)', variation_table(variations)),
        ('Shared-width and phi -> K Kbar rows', width_and_phi_table(adequacy)),
    ]
    return '\n\n'.join(f'## {title}\n\n{body}' for title, body in sections) + '\n'


def main():
    """Load the three ledgers and print the summary."""
    ledgers = [json.loads((ROOT/path).read_text()) for path in (ADEQUACY, EXCLUSION, VARIATIONS)]
    sys.stdout.write(summary(*ledgers))


if __name__ == '__main__':
    main()
