#!/usr/bin/env python
"""
Merge DownstreameR_report.xlsx (Mix_vs_HMOs sheet) with dbCAN overview.
Extract catalytic CAZy families (GH/GT/CE/PL) from HMMER and add them as columns.
"""

from __future__ import annotations

import re
from collections import OrderedDict
from pathlib import Path

import pandas as pd
import csv


OVERVIEW_FILE = r"Data/dbCAN3/overview_Mix_vs_HMOs.txt"
EXCEL_FILE = r"Data/DownstreameR_report.xlsx"
EXCEL_SHEET = "Mix_vs_HMOs"
OUTPUT_FILE = r"Data/dbCAN3/Mix_vs_HMOs.merged_dbcan.csv"
JOIN_COL = "protein_accession"

UNIPROT_RE = re.compile(r"\|([A-Z0-9]+)\|")
FAMILY_RE = re.compile(r"(GH|GT|CE|PL)(\d+)")


def extract_accession(value: str) -> str | None:
    if value is None:
        return None
    s = str(value).strip()
    if not s:
        return None
    m = UNIPROT_RE.search(s)
    if m:
        return m.group(1)
    return s


def extract_families(text: str | None) -> dict[str, str | None]:
    if not text:
        return {"GH": None, "GT": None, "CE": None, "PL": None}
    s = str(text)
    ordered: dict[str, OrderedDict[str, None]] = {
        "GH": OrderedDict(),
        "GT": OrderedDict(),
        "CE": OrderedDict(),
        "PL": OrderedDict(),
    }
    for fam, num in FAMILY_RE.findall(s):
        ordered[fam][f"{fam}{num}"] = None
    out: dict[str, str | None] = {}
    for fam in ("GH", "GT", "CE", "PL"):
        vals = list(ordered[fam].keys())
        out[fam] = "+".join(vals) if vals else None
    return out


def main() -> None:
    overview_path = Path(OVERVIEW_FILE)
    excel_path = Path(EXCEL_FILE)

    overview = pd.read_csv(overview_path, sep="\t")
    overview.columns = [c.strip() for c in overview.columns]

    if "Gene ID" not in overview.columns:
        raise ValueError("overview file missing 'Gene ID' column")

    overview["protein_accession"] = overview["Gene ID"].apply(extract_accession)

    def pick_source(row: pd.Series) -> str | None:
        for col in ("HMMER", "DIAMOND", "dbCAN_sub"):
            value = row.get(col)
            if pd.notna(value) and str(value).strip() not in {"", "N", "NA", "-"}:
                return str(value)
        return None

    families = overview.apply(lambda row: extract_families(pick_source(row)), axis=1)
    fam_df = pd.DataFrame(list(families))
    fam_df.columns = ["gh", "gt", "ce", "pl"]
    overview = pd.concat([overview, fam_df], axis=1)

    keep_cols = [
        "protein_accession",
        "gh",
        "gt",
        "ce",
        "pl",
        "HMMER",
        "dbCAN_sub",
        "DIAMOND",
        "Signalp",
        "#ofTools",
        "EC#",
    ]
    keep_cols = [c for c in keep_cols if c in overview.columns]
    overview = overview[keep_cols]

    expr = pd.read_excel(excel_path, sheet_name=EXCEL_SHEET)
    expr.columns = [str(c).strip() for c in expr.columns]
    if JOIN_COL not in expr.columns:
        raise ValueError(f"Excel sheet missing '{JOIN_COL}' column")

    expr[JOIN_COL] = expr[JOIN_COL].astype(str).str.split(";").str[0]

    merged = expr.merge(overview, how="left", left_on=JOIN_COL, right_on="protein_accession")

    merged.to_csv(
        OUTPUT_FILE,
        index=False,
        quoting=csv.QUOTE_ALL,   # quote every field
        escapechar="\\",
        encoding="utf-8",
    )

    print(f"Wrote: {OUTPUT_FILE} (rows={len(merged)})")


if __name__ == "__main__":
    main()
