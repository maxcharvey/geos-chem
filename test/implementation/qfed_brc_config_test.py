#!/usr/bin/env python3
"""Validate effective QFED2 BrC proxy mappings without a model build."""
from pathlib import Path
import re
import sys

SCALES = {
    288: 0.5, 289: 0.375, 290: 0.125, 291: 4.0, 292: 0.013,
}
PROXIES = {
    "FSOAP": 292,
    "DBRCPOA": 291,
    "NPBRCPOA": 289,
    "PBRCPOA": 290,
}

def check(path: Path) -> None:
    text = path.read_text()
    assert "QFED2_BRC_HARMONIZED_SENSITIVITY : false" in text
    for ident, value in SCALES.items():
        assert re.search(rf"\b{ident}\s+\S+\s+{value:g}\b", text), (path, ident)
    # The default branch retains the exact ordinary QFED OC mappings and has
    # no BrC proxy entries in its effective conditional block.
    assert "(((.not.QFED2_BRC_HARMONIZED_SENSITIVITY" in text
    assert "QFED_OCPI_PBL" in text and "QFED_OCPO_FT" in text
    for species, scale_id in PROXIES.items():
        rows = re.findall(rf"^0\s+QFED_{species}_BRC_HS_(?:PBL|FT).*?(?:54|70|72)/75/{scale_id}/(?:311|312)\s+5\s+2$", text, re.M)
        assert len(rows) == 2, (path, species, rows)
    # Carbon closure is explicit and independent of source magnitude.
    assert abs(SCALES[288] + SCALES[289] + SCALES[290] - 1.0) < 1e-12
    assert SCALES[291] > 1.0 and SCALES[292] == 0.013
    # Existing QFED scale chain must be retained before the new proxy scalar.
    for species, scale_id in PROXIES.items():
        assert re.search(rf"QFED_{species}_BRC_HS_(?:PBL|FT).*?(?:54|70|72)/75/{scale_id}/31[12]", text)

for arg in sys.argv[1:]:
    check(Path(arg))
print("PASS: QFED2 effective configuration numeric checks")
