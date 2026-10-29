"""Focused data-boundary and provenance checks for the Table 7 transcript."""

import csv
import json
from pathlib import Path
import tempfile
import unittest

from extract_babar_kkpi_isovector_data import write_data


class PublishedTableTests(unittest.TestCase):
    def test_table_column_transition_and_endpoint(self):
        """Keep all bins, including the last partly tau-accessible bin."""
        with tempfile.TemporaryDirectory() as directory:
            csv_path, json_path = write_data(Path(directory))
            with csv_path.open() as input_file:
                rows = list(csv.DictReader(input_file))
            self.assertEqual(len(rows), 19)
            self.assertEqual(rows[0]["energy_low_GeV"], "1.40")
            self.assertEqual(rows[0]["statistical_error_nb"], "0.67")
            self.assertEqual(rows[11]["cross_section_nb"], "1.52")
            self.assertEqual(rows[12]["cross_section_nb"], "1.83")
            self.assertEqual(rows[-1]["energy_high_GeV"], "1.78")
            self.assertEqual(rows[-1]["cross_section_nb"], "1.44")
            for previous_bin, next_bin in zip(rows, rows[1:]):
                self.assertEqual(previous_bin["energy_high_GeV"], next_bin["energy_low_GeV"])
            metadata = json.loads(json_path.read_text())
            self.assertEqual(metadata["dalitz_fit_systematic_relative"], 0.15)
            self.assertIn("Assumption", metadata["dalitz_fit_systematic_correlation"]["optional_model_status"])


if __name__ == "__main__":
    unittest.main()
