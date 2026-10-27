"""Extract BABAR:2021cde Tables IV/V and published low-energy covariances.

Run: .venv/bin/python scripts/extract_babar_2021_three_pion.py
Tests: .venv/bin/python -m unittest discover -s scripts -p extract_babar_2021_three_pion.py

Inputs are locally downloaded original LaTeX archive and manuscript-linked
covariance files. Mass units GeV, cross sections nb, covariances nb^2.
Low/high energy tables deliberately retain different error schemas. No
cross-section fit, vacuum-polarization removal or covariance combination is
performed. The tabulated correlated systematic error is preserved as a
vector; treating it as one Gaussian nuisance would be a further fit choice.
"""

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import re
import tarfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT/'memory/sources/BABAR:2021cde'
OUTPUT = ROOT/'data/odd_parity'
NUMBER = r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?'


def table_rows(source_text, label, error_count):
    """Read both printed table columns and order bins by increasing mass."""
    start = source_text.index('\\label{'+label+'}')
    content = source_text[start:source_text.index('\\end{tabular}', start)]
    pattern = (rf'({NUMBER})\s*--\s*({NUMBER})\s*&\s*\$?\s*({NUMBER})'
               + rf'\s*\$?\\pm\$?\s*({NUMBER})'*error_count)
    rows = sorted(tuple(float(value) for value in match.groups())
                  for match in re.finditer(pattern, content))
    if not rows:
        raise ValueError('No numerical rows parsed from '+label)
    values = np.array(rows)
    if np.any(values[:, 1] <= values[:, 0]) or not np.allclose(values[:-1, 1], values[1:, 0], atol=1e-12, rtol=0):
        raise ValueError('Noncontiguous or invalid bins in '+label)
    return values


def read_covariance(path, bin_count):
    """Read scalar-per-line source matrix after its explanatory header."""
    entries = [float(line.strip()) for line in path.read_text().splitlines()
               if re.fullmatch(NUMBER, line.strip())]
    if len(entries) != bin_count**2:
        raise ValueError(f'{path.name}: expected {bin_count**2} entries, found {len(entries)}')
    matrix = np.array(entries).reshape(bin_count, bin_count)
    if not np.isfinite(matrix).all() or not np.allclose(matrix, matrix.T, rtol=1e-9, atol=1e-12):
        raise ValueError(path.name+': covariance is nonfinite or asymmetric')
    eigenvalues = np.linalg.eigvalsh(matrix)
    if eigenvalues[0] < -1e-10*max(1., eigenvalues[-1]):
        raise ValueError(path.name+': covariance is not positive semidefinite')
    return matrix


def validate_covariance_diagonal(matrix, quoted_error):
    """Compare exact matrix errors to the table's 0.001 nb rounding."""
    difference = np.max(np.abs(np.sqrt(np.diag(matrix))-quoted_error))
    if difference > 0.000501:
        raise ValueError(f'Covariance diagonal differs from printed uncertainty by {difference} nb')
    return float(difference)


def write_table(path, header, rows):
    """Write round-trip-safe numeric CSV with a fixed column order."""
    stream = io.StringIO(newline='')
    writer = csv.writer(stream, lineterminator='\n')
    writer.writerow(header)
    writer.writerows([[format(value, '.17g') for value in row] for row in rows])
    path.write_text(stream.getvalue())


def extract(source_directory=SOURCE, output_directory=OUTPUT):
    """Validate sources, write separate tables/matrices and provenance manifest."""
    source_directory, output_directory = Path(source_directory), Path(output_directory)
    archive_path = source_directory/'arxiv-source.tar'
    with tarfile.open(archive_path) as archive:
        latex_bytes = archive.extractfile('babar-pub-2104.tex').read()
    source_text = latex_bytes.decode('utf-8')
    low = table_rows(source_text, 'table_cross', 3)
    high = table_rows(source_text, 'tab_cshigh', 2)
    if low.shape != (114, 6) or high.shape != (72, 5):
        raise ValueError(f'Unexpected table dimensions: {low.shape}, {high.shape}')
    statistical = read_covariance(source_directory/'rec_cov_cr.dat', len(low))
    systematic = read_covariance(source_directory/'matrix_cov_syst_sum.dat', len(low))
    diagonal_checks = dict(statistical_nb=validate_covariance_diagonal(statistical, low[:, 3]),
                           residual_systematic_nb=validate_covariance_diagonal(systematic, low[:, 5]))
    output_directory.mkdir(parents=True, exist_ok=True)
    names = {
        'low': 'babar_2021_3pi_below_1p1.csv',
        'high': 'babar_2021_3pi_above_1p1.csv',
        'statistical': 'babar_2021_3pi_below_1p1_statistical_covariance.csv',
        'residual_systematic': 'babar_2021_3pi_below_1p1_residual_systematic_covariance.csv'}
    write_table(output_directory/names['low'], ['mass_low_GeV', 'mass_high_GeV', 'cross_section_nb',
        'statistical_error_nb', 'correlated_systematic_error_nb', 'residual_systematic_error_nb'], low)
    write_table(output_directory/names['high'], ['mass_low_GeV', 'mass_high_GeV', 'cross_section_nb',
        'statistical_error_nb', 'systematic_error_nb'], high)
    for name, matrix in [('statistical', statistical), ('residual_systematic', systematic)]:
        np.savetxt(output_directory/names[name], matrix, delimiter=',', fmt='%.17g')
    checksum = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    manifest = {
        'bibtex_key': 'BABAR:2021cde',
        'method': 'Parse original source Tables IV/V; retain all published error components separately.',
        'source_files_sha256': {path.name: checksum(path) for path in [archive_path,
            source_directory/'rec_cov_cr.dat', source_directory/'matrix_cov_syst_sum.dat']},
        'latex_member': 'babar-pub-2104.tex',
        'latex_member_sha256': hashlib.sha256(latex_bytes).hexdigest(),
        'source_urls': {'archive': 'https://arxiv.org/src/2110.00520',
            'statistical_covariance': 'https://raw.githubusercontent.com/eakozyrev/3pi_unfold/main/rec_cov_cr.dat',
            'residual_systematic_covariance': 'https://raw.githubusercontent.com/eakozyrev/3pi_unfold/main/matrix_cov_syst_sum.dat'},
        'units': {'mass': 'GeV', 'cross_section': 'nb', 'covariance': 'nb^2'},
        'low_energy': {'table': 'IV', 'bins': len(low), 'range_GeV': [float(low[0, 0]), float(low[-1, 1])],
            'covariance_order': 'rows/columns ordered by increasing mass, identical to below_1p1.csv',
            'statistical_covariance': names['statistical'],
            'residual_systematic_covariance': names['residual_systematic'],
            'correlated_systematic': 'Second printed error retained as a vector; includes luminosity, radiative correction, efficiency and unfolding. No rank-one covariance or independent-bin approximation imposed.',
            'residual_systematic': 'Third printed error: background subtraction and Lorentzian/Gaussian resolution uncertainties; use supplied covariance.',
            'maximum_diagonal_error_rounding_difference': diagonal_checks},
        'high_energy': {'table': 'V', 'bins': len(high), 'range_GeV': [float(high[0, 0]), float(high[-1, 1])],
            'errors': 'First statistical, second total systematic. No high-energy covariance is supplied here. For its a_mu integral over 1.1-2.0 GeV, the paper conservatively treats systematics as 100% correlated; this is not imposed on the full table.',
            'nonresonant_bins': '3.0-3.1 and 3.1-3.2 GeV have J/psi contribution subtracted.'},
        'cross_region_correlations': 'Luminosity, radiative correction and efficiency uncertainties are shared between regions. No cross-region covariance is invented.',
        'cross_section_definition': 'Published Born cross sections; no additional undressing or final-state-radiation conversion performed.',
        'outputs_sha256': {name: checksum(output_directory/name) for name in names.values()},
        'invocation': '.venv/bin/python scripts/extract_babar_2021_three_pion.py'}
    manifest_path = output_directory/'babar_2021_3pi_provenance.json'
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True)+'\n')
    return [output_directory/name for name in names.values()]+[manifest_path]


class ExtractionTests(unittest.TestCase):
    """Check independent bin/error facts and covariance rounding agreement."""

    def test_source_tables_and_covariance(self):
        with tarfile.open(SOURCE/'arxiv-source.tar') as archive:
            text = archive.extractfile('babar-pub-2104.tex').read().decode()
        low, high = table_rows(text, 'table_cross', 3), table_rows(text, 'tab_cshigh', 2)
        self.assertEqual(low.shape, (114, 6))
        self.assertEqual(high.shape, (72, 5))
        np.testing.assert_allclose(low[0], [.620, .630, .192, .086, .036, .018])
        np.testing.assert_allclose(low[-1], [1.09, 1.1, 3.920, .243, .051, .170])
        np.testing.assert_allclose(high[0], [1.1, 1.125, 4.29, .23, .22])
        np.testing.assert_allclose(high[-1], [3.4, 3.5, .031, .016, .003])
        for filename, column in [('rec_cov_cr.dat', 3), ('matrix_cov_syst_sum.dat', 5)]:
            matrix = read_covariance(SOURCE/filename, 114)
            self.assertLess(validate_covariance_diagonal(matrix, low[:, column]), .000501)

    def test_noncontiguous_table_rejected(self):
        malformed = r'\label{example} 1--2 & 3 $\pm$ 1 & 3--4 & 5 $\pm$ 2 \end{tabular}'
        with self.assertRaisesRegex(ValueError, 'Noncontiguous'):
            table_rows(malformed, 'example', 1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-directory', type=Path, default=SOURCE)
    parser.add_argument('--output-directory', type=Path, default=OUTPUT)
    arguments = parser.parse_args()
    for path in extract(arguments.source_directory, arguments.output_directory):
        print(path.relative_to(ROOT) if path.is_relative_to(ROOT) else path)


if __name__ == '__main__':
    main()
