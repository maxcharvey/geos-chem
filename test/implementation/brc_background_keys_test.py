#!/usr/bin/env python3
"""Regression: BrC restart defaults must use the species reader's actual key."""
from pathlib import Path
import re
import unittest

ROOT = Path(__file__).resolve().parents[2]


class BackgroundKeysTest(unittest.TestCase):
    def test_reader_and_template_agree(self):
        text = (ROOT / 'run/shared/species_database.yml').read_text()
        self.assertNotIn('BackgroundVV:', text)
        blocks = dict(re.findall(r'^([A-Za-z0-9_]+):\n(.*?)(?=^[A-Za-z0-9_]+:|\Z)', text, re.M | re.S))
        expected = {n: 1e-18 for n in ['FSOAS', 'BRCSOA', 'WTC', 'DBRCPOA', 'NPBRCPOA', 'PBRCPOA']}
        expected['FSOAP'] = 1e-20
        for name, value in expected.items():
            with self.subTest(species=name):
                keys = re.findall(r'^  Background_VV:\s*(\S+)', blocks[name], re.M)
                self.assertEqual(len(keys), 1)
                self.assertEqual(float(keys[0]), value)
        reader = (ROOT / 'Headers/species_database_mod.F90').read_text()
        self.assertIn('Background_VV', reader)


if __name__ == '__main__':
    unittest.main()
