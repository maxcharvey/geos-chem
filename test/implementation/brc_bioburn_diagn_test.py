#!/usr/bin/env python3
"""Static contract checks for BrC biomass-burning diagnostic templates."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
TEMPLATES = (
    ROOT / "run/GCClassic/HEMCO_Diagn.rc.templates/HEMCO_Diagn.rc.fullchem",
    ROOT / "run/GCClassic/HEMCO_Diagn.rc.templates/HEMCO_Diagn.rc.aerosol",
)
SPECIES = ("FSOAP", "DBRCPOA", "NPBRCPOA", "PBRCPOA")
DEFAULT_SELECTOR = ("112", "-1", "-1")


def entries(path):
    """Return active diagnostic entries keyed by name."""
    result = {}
    for line in path.read_text().splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        fields = stripped.split()
        if len(fields) < 7 or not fields[0].startswith("Emis"):
            continue
        if fields[0] in result:
            raise AssertionError(f"{path}: duplicate diagnostic {fields[0]}")
        result[fields[0]] = fields
    return result


def check_template(path):
    text = path.read_text()
    for selector in ("111/-1/-1", "165/-1/-1", "0/5/3", "0/5/2"):
        assert selector in text, f"{path}: missing selector guidance {selector}"
    assert "does not\n# rewrite these selectors" in text, (
        f"{path}: must state that run-directory creation does not rewrite selectors"
    )

    records = entries(path)
    observed_selectors = set()
    for species in SPECIES:
        column = records[f"Emis{species}_BioBurn"]
        profile = records[f"Emis{species}_BioBurn3D"]
        assert column[1] == profile[1] == species, f"{path}: {species} species mismatch"
        assert tuple(column[2:5]) == DEFAULT_SELECTOR, f"{path}: {species} column selector"
        assert tuple(profile[2:5]) == DEFAULT_SELECTOR, f"{path}: {species} profile selector"
        assert column[5] == "2", f"{path}: {species} BioBurn must be two-dimensional"
        assert profile[5] == "3", f"{path}: {species} BioBurn3D must be three-dimensional"
        assert column[6] == profile[6] == "kg/m2/s", f"{path}: {species} units"
        observed_selectors.add(tuple(column[2:5]))

        total = records[f"Emis{species}_Total"]
        total_column = records[f"Emis{species}_TotalColumn"]
        assert tuple(total[2:5]) == tuple(total_column[2:5]) == ("-1", "-1", "-1")
        assert total[5] == "3" and total_column[5] == "2", (
            f"{path}: {species} total diagnostic dimensions changed"
        )

    assert observed_selectors == {DEFAULT_SELECTOR}, f"{path}: inconsistent BioBurn source"


def main():
    for template in TEMPLATES:
        check_template(template)
        print(f"PASS {template.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
