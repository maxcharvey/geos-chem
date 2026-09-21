#!/usr/bin/env python3
"""Parse effective QFED2 rows for default and BrC-sensitivity configurations."""
from dataclasses import dataclass
from pathlib import Path
import re
import sys

FIELDS = ("source", "variable", "date", "cycle", "dimension", "units",
          "species", "scales", "category", "hierarchy")
PROXY_CHAINS = {
    "FSOAP": ("54", "75", "292"), "DBRCPOA": ("75", "291"),
    "NPBRCPOA": ("75", "289"), "PBRCPOA": ("75", "290"),
    "OCPI": ("72", "75", "288"), "OCPO": ("73", "75", "288"),
}


@dataclass(frozen=True)
class Row:
    name: str
    source: str
    variable: str
    date: str
    cycle: str
    dimension: str
    units: str
    species: str
    scales: tuple[str, ...]
    category: str
    hierarchy: str


def condition_active(condition: str, flags: dict[str, bool]) -> bool:
    condition = condition.strip()
    return not flags[condition[5:]] if condition.startswith(".not.") else flags[condition]


def qfed_section(text: str) -> list[str]:
    match = re.search(r"# --- QFED2 biomass burning \(v2\.5r1\) ---\n"
                      r"(.*?)^\)\)\)QFED2\s*$", text, re.M | re.S)
    assert match, "missing QFED2 block"
    return match.group(1).splitlines()


def effective_rows(text: str, harmonized: bool) -> dict[str, Row]:
    """Evaluate QFED2 conditionals and resolve HEMCO '-' row inheritance."""
    flags = {"QFED2": True, "QFED2_BRC_HARMONIZED_SENSITIVITY": harmonized}
    active, previous, rows = [True], {}, {}
    for raw in qfed_section(text):
        line = raw.strip()
        if line.startswith("((("):
            active.append(active[-1] and condition_active(line[3:], flags))
            continue
        if line.startswith(")))"):
            assert len(active) > 1, line
            active.pop()
            continue
        if not active[-1] or not line.startswith("0 QFED_"):
            continue
        tokens = line.split()
        assert len(tokens) == 12, line
        values = dict(zip(FIELDS, tokens[2:]))
        for field in FIELDS[:6]:
            if values[field] == "-":
                assert field in previous, (line, field)
                values[field] = previous[field]
        previous = values.copy()
        row = Row(tokens[1], values["source"], values["variable"], values["date"],
                  values["cycle"], values["dimension"], values["units"],
                  values["species"], tuple(values["scales"].split("/")),
                  values["category"], values["hierarchy"])
        assert row.name not in rows, row.name
        rows[row.name] = row
    return rows


def scalar_factors(text: str) -> dict[str, float]:
    return {ident: float(value) for ident, value in re.findall(
        r"^\s*(\d+)\s+\S+\s+([0-9.]+)\s+-\s+-\s+-\s+xy", text, re.M)}


def scaled(row: Row, factors: dict[str, float], tod: float) -> float:
    value = 1.
    for factor in row.scales:
        value *= tod if factor == "75" else factors[factor]
    return value


def assert_row(row: Row, target: str, chain: tuple[str, ...], level: str) -> None:
    assert row.species == target and row.scales[:-1] == chain and row.scales[-1] == level, row
    assert row.units == "kg/m2/s" and row.category == "5" and row.hierarchy == "2", row
    assert row.dimension == ("xyL=1:PBL" if level == "311" else "xyL=PBL:5500m"), row
    source = "co" if target == "FSOAP" else "bc" if target == "DBRCPOA" else "oc"
    assert row.source.endswith(f"qfed2.emis_{source}.006.$YYYY$MM$DD.nc4"), row


def check(text: str, has_pog: bool) -> None:
    factors = scalar_factors(text)
    expected = {"54": 1.05, "70": .2, "71": .8, "72": .5, "73": .5,
                "74": 1.27, "76": .49, "77": .51, "288": .5, "289": .375,
                "290": .125, "291": 4., "292": .013, "311": .65, "312": .35}
    for ident, value in expected.items():
        assert factors[ident] == value, (ident, factors.get(ident))
    off, on = effective_rows(text, False), effective_rows(text, True)
    assert not any("BRC_HS" in name for name in off), off
    for name, species in (("QFED_OCPI_PBL", "OCPI"), ("QFED_OCPO_FT", "OCPO")):
        assert off[name].species == species
    for species, chain in PROXY_CHAINS.items():
        for suffix, level in (("PBL", "311"), ("FT", "312")):
            assert_row(on[f"QFED_{species}_BRC_HS_{suffix}"], species, chain, level)
    if has_pog:
        for suffix, level in (("PBL", "311"), ("FT", "312")):
            for number, split in (("1", "76"), ("2", "77")):
                name = f"QFED_POG{number}_{suffix}"
                assert off[name] == on[name], name
                assert_row(on[name], f"POG{number}", ("74", split, "75"), level)

    # A non-unit TOD factor proves every chain uses the operational scale 75.
    tod, pbl = 2.3, factors["311"]
    oc_rows = (on["QFED_OCPI_BRC_HS_PBL"], on["QFED_OCPO_BRC_HS_PBL"],
               on["QFED_NPBRCPOA_BRC_HS_PBL"], on["QFED_PBRCPOA_BRC_HS_PBL"])
    assert abs(sum(scaled(row, factors, tod) for row in oc_rows) - tod * pbl) < 1e-12
    assert abs(scaled(on["QFED_DBRCPOA_BRC_HS_PBL"], factors, tod) - 4. * tod * pbl) < 1e-12
    assert abs(scaled(on["QFED_FSOAP_BRC_HS_PBL"], factors, tod) - .013 * factors["54"] * tod * pbl) < 1e-12


def must_fail(mutated: str, has_pog: bool) -> None:
    try:
        check(mutated, has_pog)
    except (AssertionError, KeyError):
        return
    raise AssertionError("deliberate proxy-chain mutation passed")


def check_mutations(text: str, has_pog: bool) -> None:
    for before, after in (("DBRCPOA 75/291/311", "DBRCPOA 70/291/311"),
                          ("NPBRCPOA 75/289/311", "NPBRCPOA 72/289/311"),
                          ("FSOAP 54/75/292/311", "FSOAP 54/292/311")):
        assert before in text, before
        must_fail(text.replace(before, after, 1), has_pog)


if len(sys.argv) < 2:
    raise SystemExit("usage: qfed_brc_config_test.py HEMCO_Config.rc [...]")

for arg in sys.argv[1:]:
    source = Path(arg).read_text()
    has_pog = Path(arg).name.endswith("fullchem")
    check(source, has_pog)
    check_mutations(source, has_pog)
print("PASS: QFED2 parsed effective-row and numeric-chain checks")
