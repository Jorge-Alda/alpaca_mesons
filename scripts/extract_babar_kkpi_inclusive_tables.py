"""Regenerate BaBar:2007ceh inclusive Tables 2/3 through 1.78 GeV.

Run: .venv/bin/python scripts/extract_babar_kkpi_inclusive_tables.py
The statistical errors are published values; total corrections are already
applied. The Dalitz-component 15% systematic is not assigned to these tables.
"""
import json
from pathlib import Path

DATA = {'babar_2007_kskpi_inclusive_below_1p78': {'csv': 'energy_low_GeV,energy_high_GeV,cross_section_nb,statistical_error_nb\n'
                                                  '1.24,1.28,0.037,0.034\n'
                                                  '1.28,1.32,0.042,0.024\n'
                                                  '1.32,1.36,0.014,0.014\n'
                                                  '1.36,1.38,0.107,0.054\n'
                                                  '1.38,1.40,0.132,0.059\n'
                                                  '1.40,1.42,0.470,0.111\n'
                                                  '1.42,1.44,0.542,0.118\n'
                                                  '1.44,1.46,0.990,0.167\n'
                                                  '1.46,1.48,1.279,0.187\n'
                                                  '1.48,1.50,1.231,0.180\n'
                                                  '1.50,1.52,1.673,0.205\n'
                                                  '1.52,1.54,1.964,0.220\n'
                                                  '1.54,1.56,2.564,0.249\n'
                                                  '1.56,1.58,2.440,0.253\n'
                                                  '1.58,1.60,3.269,0.282\n'
                                                  '1.60,1.62,3.397,0.283\n'
                                                  '1.62,1.64,2.987,0.268\n'
                                                  '1.64,1.66,3.807,0.302\n'
                                                  '1.66,1.68,3.760,0.291\n'
                                                  '1.68,1.70,3.498,0.284\n'
                                                  '1.70,1.72,3.046,0.262\n'
                                                  '1.72,1.74,2.686,0.247\n'
                                                  '1.74,1.76,2.319,0.223\n'
                                                  '1.76,1.78,1.500,0.181\n',
                                           'metadata': {'bibtex_key': 'BaBar:2007ceh',
                                                        'arxiv_id': '0710.4451v3',
                                                        'source_url': 'https://arxiv.org/html/0710.4451v3',
                                                        'source_table': 2,
                                                        'source_section': 'IV.2',
                                                        'observable': 'inclusive e+e- -> K0S K+ '
                                                                      'pi- plus charge conjugate '
                                                                      'cross section',
                                                        'energy_unit': 'GeV',
                                                        'cross_section_unit': 'nb',
                                                        'number_of_bins': 24,
                                                        'selection': 'All Table 2 bins from '
                                                                     'threshold through 1.76-1.78 '
                                                                     'GeV',
                                                        'bin_interpretation': 'Published cross '
                                                                              'section bin '
                                                                              'estimates, not '
                                                                              'point values or bin '
                                                                              'integrals',
                                                        'statistical_errors': 'Table 2 quoted '
                                                                              'errors; bin '
                                                                              'covariance '
                                                                              'unavailable',
                                                        'systematic_relative': 0.051,
                                                        'systematic_source': 'Table 1 total for '
                                                                             'inclusive K0S K+/- '
                                                                             'pi-/+ cross section',
                                                        'systematic_correlation': 'No bin '
                                                                                  'correlation '
                                                                                  'matrix '
                                                                                  'supplied; fit '
                                                                                  'uses an '
                                                                                  'explicit common '
                                                                                  'normalization '
                                                                                  'nuisance and a '
                                                                                  'diagonal '
                                                                                  'sensitivity '
                                                                                  'variant',
                                                        'published_corrections': 'The +4.2% '
                                                                                 'correction is '
                                                                                 'already applied; '
                                                                                 'do not apply it '
                                                                                 'again',
                                                        'fit_qualification': 'Inclusive '
                                                                             'charge-summed cross '
                                                                             'section; suitable '
                                                                             'for coherent full K '
                                                                             'Kbar pi amplitude, '
                                                                             'distinct from the '
                                                                             'Table 7 K*K isobar '
                                                                             'component'}},
 'babar_2007_kkpi0_inclusive_below_1p78': {'csv': 'energy_low_GeV,energy_high_GeV,cross_section_nb,statistical_error_nb\n'
                                                  '1.36,1.38,0.208,0.120\n'
                                                  '1.38,1.40,0.100,0.089\n'
                                                  '1.40,1.42,0.042,0.059\n'
                                                  '1.42,1.44,0.241,0.092\n'
                                                  '1.44,1.46,0.193,0.091\n'
                                                  '1.46,1.48,0.370,0.093\n'
                                                  '1.48,1.50,0.688,0.122\n'
                                                  '1.50,1.52,0.385,0.092\n'
                                                  '1.52,1.54,0.485,0.102\n'
                                                  '1.54,1.56,0.610,0.113\n'
                                                  '1.56,1.58,0.556,0.109\n'
                                                  '1.58,1.60,0.911,0.130\n'
                                                  '1.60,1.62,1.030,0.132\n'
                                                  '1.62,1.64,0.882,0.127\n'
                                                  '1.64,1.66,0.767,0.116\n'
                                                  '1.66,1.68,0.749,0.114\n'
                                                  '1.68,1.70,0.413,0.082\n'
                                                  '1.70,1.72,0.296,0.071\n'
                                                  '1.72,1.74,0.093,0.046\n'
                                                  '1.74,1.76,0.116,0.046\n'
                                                  '1.76,1.78,0.073,0.038\n',
                                           'metadata': {'bibtex_key': 'BaBar:2007ceh',
                                                        'arxiv_id': '0710.4451v3',
                                                        'source_url': 'https://arxiv.org/html/0710.4451v3',
                                                        'source_table': 3,
                                                        'source_section': 'V',
                                                        'observable': 'inclusive e+e- -> K+ K- pi0 '
                                                                      'cross section',
                                                        'energy_unit': 'GeV',
                                                        'cross_section_unit': 'nb',
                                                        'number_of_bins': 21,
                                                        'selection': 'All Table 3 bins from '
                                                                     'threshold through 1.76-1.78 '
                                                                     'GeV',
                                                        'bin_interpretation': 'Published cross '
                                                                              'section bin '
                                                                              'estimates, not '
                                                                              'point values or bin '
                                                                              'integrals',
                                                        'statistical_errors': 'Table 3 quoted '
                                                                              'errors; bin '
                                                                              'covariance '
                                                                              'unavailable',
                                                        'systematic_relative': 0.065,
                                                        'systematic_source': 'Table 4 total for '
                                                                             'inclusive K+ K- pi0 '
                                                                             'cross section',
                                                        'systematic_correlation': 'No bin '
                                                                                  'correlation '
                                                                                  'matrix '
                                                                                  'supplied; fit '
                                                                                  'uses an '
                                                                                  'explicit common '
                                                                                  'normalization '
                                                                                  'nuisance and a '
                                                                                  'diagonal '
                                                                                  'sensitivity '
                                                                                  'variant',
                                                        'published_corrections': 'The +4.6% '
                                                                                 'correction is '
                                                                                 'already applied; '
                                                                                 'do not apply it '
                                                                                 'again',
                                                        'fit_qualification': 'Inclusive '
                                                                             'charge-conjugation-even '
                                                                             'final state; fit '
                                                                             'with the coherent '
                                                                             'full K Kbar pi '
                                                                             'amplitude, distinct '
                                                                             'from the Table 7 K*K '
                                                                             'isobar component'}}}


def write_data(directory):
    directory.mkdir(parents=True,exist_ok=True)
    paths=[]
    for basename,record in DATA.items():
        csv_path=directory/(basename+'.csv')
        csv_path.write_text(record['csv'])
        json_path=directory/(basename+'.json')
        json_path.write_text(json.dumps(record['metadata'],indent=2)+'\n')
        paths.extend((csv_path,json_path))
    return paths


if __name__=='__main__':
    root=Path(__file__).resolve().parents[1]
    for path in write_data(root/'data/vector_sector'):
        print(path.relative_to(root))
