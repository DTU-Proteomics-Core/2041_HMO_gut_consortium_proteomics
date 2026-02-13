#!/usr/bin/env python
"""
Create a volcano plot from a pre-merged HMOs_vs_Fiber dbCAN file.
Only GH/CE/PL families are shown (GT excluded).
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from adjustText import adjust_text


UNIPROT_RE = re.compile(r"\b(?:sp|tr)?\|?([A-NR-Z0-9]{5,10})\|?")


def clean_names(cols: list[str]) -> list[str]:
    out = []
    for c in cols:
        c = c.strip().lower()
        c = re.sub(r"[\s\.-]+", "_", c)
        c = re.sub(r"__+", "_", c)
        c = re.sub(r"^_+|_+$", "", c)
        out.append(c)
    return out


def extract_accession(value: str | float | int) -> str | None:
    if pd.isna(value):
        return None
    s = str(value).strip()
    if not s:
        return None
    m = UNIPROT_RE.search(s)
    if m:
        return m.group(1)
    m2 = re.search(r"\b([A-NR-Z0-9]{5,10})\b", s)
    return m2.group(1) if m2 else None


def find_accession_column(df: pd.DataFrame) -> str:
    candidates = {
        "protein_accession",
        "accession",
        "acc",
        "protein_id",
        "id",
        "protein",
        "geneid",
        "gene_id",
        "query",
        "query_name",
    }
    for col in df.columns:
        if col in candidates:
            return col
    # fallback: pick column with most extractable accessions
    best_col = df.columns[0]
    best_count = -1
    for col in df.columns:
        count = df[col].apply(extract_accession).notna().sum()
        if count > best_count:
            best_col = col
            best_count = count
    return best_col


def load_expression(path: Path) -> pd.DataFrame:
    if path.suffix.lower() in {".xlsx", ".xls"}:
        df = pd.read_excel(path)
    else:
        try:
            df = pd.read_csv(path)
        except pd.errors.ParserError:
            df = pd.read_csv(path, engine="python", on_bad_lines="warn")

    df.columns = clean_names([str(c) for c in df.columns])

    # canonicalize common column names
    alias_map = {
        "protein_accession": {"protein_accession", "accession", "acc", "protein_id", "id"},
        "logfc": {"logfc", "log2foldchange", "log2_fc", "log2fold_change", "fold_change"},
        "fdr": {"fdr", "padj", "adj_pval", "p_adj", "p.adj", "qvalue"},
    }
    for canonical, aliases in alias_map.items():
        if canonical in df.columns:
            continue
        for col in df.columns:
            if col in aliases:
                df = df.rename(columns={col: canonical})
                break

    acc_col = find_accession_column(df)
    df["protein_accession_raw"] = df[acc_col].astype(str)
    df["protein_accession"] = df["protein_accession_raw"].apply(extract_accession)

    if "logfc" not in df.columns:
        # fallback: any column containing 'log' or 'fold'
        for col in df.columns:
            if re.search(r"log|fold", col):
                df["logfc"] = pd.to_numeric(df[col], errors="coerce")
                break
    if "fdr" not in df.columns:
        for col in df.columns:
            if re.search(r"fdr|padj|qvalue|p_adj|p.adj", col):
                df["fdr"] = pd.to_numeric(df[col], errors="coerce")
                break

    return df


def build_cazy_from_columns(df: pd.DataFrame) -> pd.Series:
    if not any(col in df.columns for col in ("gh", "ce", "pl")):
        return df.get("cazy")

    def combine_row(row: pd.Series) -> str | None:
        parts: list[str] = []
        for col in ("gh", "ce", "pl"):
            value = row.get(col)
            if pd.isna(value):
                continue
            for token in str(value).split("+"):
                token = token.strip()
                if token and token not in parts:
                    parts.append(token)
        return "+".join(parts) if parts else None

    return df.apply(combine_row, axis=1)


def parse_species_from_fasta(value: str | float | int) -> str | None:
    if pd.isna(value):
        return None
    s = str(value).strip()
    if not s:
        return None
    match = re.search(r"([A-Z][a-z]+_[a-z0-9]+)", s)
    if match:
        return match.group(1).replace("_", " ")
    tokens = [t for t in re.split(r"[_/\\|\s]+", s) if t]
    if not tokens:
        return None
    last = tokens[-1].replace("_", " ")
    if last.lower() in {"id", "taxonomy", "taxonomy id"}:
        return None
    return last


def infer_species(df: pd.DataFrame) -> pd.Series:
    if "strain_species" in df.columns:
        series = df["strain_species"]
    elif "species" in df.columns:
        series = df["species"]
    else:
        fasta_col = None
        for col in ("fasta_file", "protein_file", "fasta"):
            if col in df.columns:
                fasta_col = col
                break
        if fasta_col:
            series = df[fasta_col].apply(parse_species_from_fasta)
        else:
            series = pd.Series([None] * len(df), index=df.index)

    series = series.astype("string")
    series = series.where(series.notna() & (series.str.strip() != ""), None)
    series = series.where(series.str.lower() != "nan", None)
    return series


def load_species_colors(path: Path) -> tuple[dict[str, str], list[str]]:
    df = pd.read_csv(path, sep="\t")
    df.columns = clean_names([str(c) for c in df.columns])
    if "species" not in df.columns or "color" not in df.columns:
        raise ValueError("species colors file must have 'species' and 'color' columns")
    mapping: dict[str, str] = {}
    order: list[str] = []
    for _, row in df.iterrows():
        species = str(row["species"]).strip()
        color = str(row["color"]).strip()
        if not species or not color:
            continue
        if species.lower() == "other":
            continue
        mapping[species] = color
        order.append(species)
    return mapping, order


def darken_hex(color: str, factor: float = 0.75) -> str:
    if not color:
        return color
    if not isinstance(color, str):
        return color
    s = color.strip()
    if not s.startswith("#") or len(s) != 7:
        return color
    try:
        r = int(s[1:3], 16)
        g = int(s[3:5], 16)
        b = int(s[5:7], 16)
    except ValueError:
        return color
    r = max(0, min(255, int(r * factor)))
    g = max(0, min(255, int(g * factor)))
    b = max(0, min(255, int(b * factor)))
    return f"#{r:02X}{g:02X}{b:02X}"


def format_species_label(species: str) -> str:
    if species == "Bifidobacterium longum":
        return r"$\it{Bifidobacterium\ longum}$ subsp. $\it{infantis}$"
    return rf"$\it{{{species.replace(' ', '\\ ')}}}$"


def volcano_plot(
    df: pd.DataFrame,
    outpath: Path,
    lfc_cut: float,
    fdr_cut: float,
    label_n: int,
    title: str,
    xlabel: str,
    species_colors_path: Path | None = None,
) -> None:
    df = df.copy()
    if "cazy" not in df.columns or any(col in df.columns for col in ("gh", "ce", "pl")):
        df["cazy"] = build_cazy_from_columns(df)
    cazy_str = df["cazy"].astype(str).str.strip()
    cazy_upper = cazy_str.str.upper()
    bad_cazy = {"", "NA", "N/A", "NONE", "NULL", "NAN", "-"}
    df["cazy"] = cazy_str.where(~cazy_upper.isin(bad_cazy), np.nan)
    df.loc[df["cazy"].astype(str).str.fullmatch(r"GH0", na=False), "cazy"] = np.nan
    df.loc[cazy_upper.str.startswith("CBM", na=False), "cazy"] = np.nan
    df.loc[cazy_upper.str.startswith("GT", na=False), "cazy"] = np.nan
    df["neglog10_fdr"] = -np.log10(df["fdr"].astype(float))
    df["significant"] = (df["fdr"] < fdr_cut) & (df["logfc"].abs() > lfc_cut)
    df["has_cazy"] = df["cazy"].notna() & (df["cazy"].astype(str) != "")
    df["strain_species"] = infer_species(df)

    plt.figure(figsize=(10, 8))
    plt.scatter(
        df["logfc"],
        df["neglog10_fdr"],
        s=15,
        facecolors="none",
        edgecolors="#9e9e9e",
        alpha=0.5,
        linewidths=0.6,
    )

    if species_colors_path is not None:
        species_colors, species_order = load_species_colors(species_colors_path)
    else:
        species_counts = (
            df.loc[df["has_cazy"] & df["strain_species"].notna(), "strain_species"]
            .value_counts()
            .head(5)
        )
        species_to_color = list(species_counts.index)
        palette = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e"]
        species_colors = dict(zip(species_to_color, palette[: len(species_to_color)]))
        species_order = list(species_colors.keys())

    legend_handles: list[Line2D] = []
    legend_labels: list[str] = []
    colored_artists = []
    species_markers = {
        "Bifidobacterium longum": ">",
        "Bifidobacterium adolescentis": "d",
        "Roseburia inulinivorans": "p",
        "Roseburia intestinalis": "p",
        "Roseburia hominis": "p",
    }
    if species_colors:
        preferred_order = [
            "Bifidobacterium longum",
            "Bifidobacterium adolescentis",
            "Roseburia inulinivorans",
            "Roseburia intestinalis",
            "Roseburia hominis",
        ]
        ordered_species = [
            species
            for species in preferred_order
            if species in species_colors
        ]
        ordered_species.extend(
            [species for species in species_order if species not in ordered_species]
        )
        for species in ordered_species:
            color = species_colors.get(species)
            if not color:
                continue
            marker = species_markers.get(species, "o")
            subset = df[df["has_cazy"] & (df["strain_species"] == species)]
            if subset.empty:
                continue
            colored_artists.append(
                plt.scatter(
                subset["logfc"],
                subset["neglog10_fdr"],
                s=30,
                c=color,
                alpha=0.9,
                linewidths=0,
                marker=marker,
                )
            )
            legend_handles.append(
                Line2D(
                    [0],
                    [0],
                    marker=marker,
                    color="none",
                    markerfacecolor=color,
                    markeredgecolor=color,
                    markersize=6,
                )
            )
            legend_labels.append(format_species_label(species))

        if species_colors_path is not None:
            legend_handles.append(
                Line2D(
                    [0],
                    [0],
                    marker="o",
                    color="none",
                    markerfacecolor="none",
                    markeredgecolor="#9e9e9e",
                    markeredgewidth=1.5,
                    markersize=6,
                )
            )
            legend_labels.append("Non-CAZy")
   
    plt.axvline(-lfc_cut, color="#666666", linestyle="--", linewidth=1)
    plt.axvline(lfc_cut, color="#666666", linestyle="--", linewidth=1)
    plt.axhline(-np.log10(fdr_cut), color="#666666", linestyle="--", linewidth=1)

    ax = plt.gca()
    ax.set_xlim(-10, 10)
    x_min, x_max = -10, 10
    y_max = df["neglog10_fdr"].max()
    lfc_label = 0.585
    x_offset = (x_max - x_min) * 0.012
    ax.text(
        -lfc_cut - x_offset,
        y_max * 0.95,
        f"log$_2$(FC) = {-lfc_label:.3f}",
        ha="right",
        va="bottom",
        fontsize=9,
        fontweight="bold",
    )
    ax.text(
        lfc_cut + x_offset,
        y_max * 0.95,
        f"log$_2$(FC) = {lfc_label:.3f}",
        ha="left",
        va="bottom",
        fontsize=9,
        fontweight="bold",
    )
    ax.text(
        x_min + (x_max - x_min) * 0.03,
        -np.log10(fdr_cut),
        f"FDR = {fdr_cut:.2f}",
        ha="left",
        va="bottom",
        fontsize=9,
        fontweight="bold",
    )

    # label top hits by score among non-GT CAZy
    df["score"] = df["logfc"].abs() * df["neglog10_fdr"]
    cazy_str = df["cazy"].astype(str)
    non_gt = ~cazy_str.str.startswith("GT", na=False)
    top = df[df["has_cazy"] & non_gt].sort_values("score", ascending=False).head(label_n)

    gh_mask = df["has_cazy"] & cazy_str.str.contains(
        r"\bGH10\b|\bGH26\b|\bGH29\b|\bGH33\b|\bGH95\b|\bGH112\b|\bGH3\b|\bGH109\b|\bGH49\b",
        regex=True,
        na=False,
    )
    gh_targets = df[gh_mask]
    rose_inu_mask = (df["has_cazy"] & (df["strain_species"] == "Roseburia inulinivorans"))
    rose_inu_targets = df[rose_inu_mask]
    label_df = pd.concat([top, gh_targets, rose_inu_targets], ignore_index=True).drop_duplicates()
    texts = []
    label_points = []
    for _, row in label_df.iterrows():
        tokens = str(row.get("cazy") or "").split("+")
        targets = [t for t in tokens if t in {"GH10", "GH26", "GH29", "GH33", "GH95", "GH112", "GH3", "GH109", "GH49"}]
        annotation_tokens = [t for t in tokens if t.startswith(("GH", "PL", "CE"))]
        label = "+".join(targets) if targets else (row.get("cazy") or row.get("protein_accession"))
        if pd.isna(label):
            continue
        x = row["logfc"]
        y = row["neglog10_fdr"]
        label_color = "#333333"
        label_weight = "normal"
        if annotation_tokens and species_colors_path is not None:
            species = row.get("strain_species")
            if isinstance(species, str) and species in species_colors:
                label_color = darken_hex(species_colors[species], factor=0.65)
        if annotation_tokens:
            label_weight = "bold"
        texts.append(
            plt.text(
                x,
                y,
                str(label),
                fontsize=7,
                zorder=5,
                color=label_color,
                fontweight=label_weight,
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="none", alpha=0.5),
            )
        )
        label_points.append((x, y))
    if texts:
        adjust_text(texts)
        for text, (x, y) in zip(texts, label_points):
            ax.annotate(
                "",
                xy=(x, y),
                xytext=text.get_position(),
                textcoords="data",
                arrowprops=dict(arrowstyle="-", color="#888888", lw=0.6, shrinkA=6),
                zorder=4,
            )

    plt.title(title, fontsize=18, fontweight="bold")
    plt.xlabel(xlabel, fontsize=12, fontweight="bold")
    plt.ylabel(
        r"log$_{10}$(FDR)",
        fontsize=12,
        fontweight="bold",
    )

    if legend_handles:
        plt.legend(
            legend_handles,
            legend_labels,
            frameon=False,
            loc="upper center",
            bbox_to_anchor=(0.5, -0.08),
            ncol=2,
        )
    plt.tight_layout()
    plt.savefig(outpath, dpi=1000)
    plt.close()


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Plot volcano from merged dbCAN data.")
    p.add_argument("--lfc", type=float, default=0.585, help="log2FC cutoff")
    p.add_argument("--fdr", type=float, default=0.05, help="FDR cutoff")
    p.add_argument("--label-n", type=int, default=50, help="Number of CAZy hits to label")
    p.add_argument("--title", default="HMOs vs Fibre", help="Plot title")
    p.add_argument(
        "--species-colors",
        default=None,
        help="TSV file with columns: species, color (hex). Colors only those species.",
    )
    return p.parse_args()


def main() -> None:
    PLOTS = [
        {
            "expr": r"Data/dbCAN3/HMOs_vs_Fiber.merged_dbcan.csv",
            "out": r"Data/dbCAN3/HMOs_vs_Fibre_supplementary_figure_xx.png",
            "title": "HMOs vs Fibre",
            "xlabel": r"Growth on Fibre  $\leftarrow$  log$_2$(FC)  $\rightarrow$  Growth on HMOs",
        },
        {
            "expr": r"Data/dbCAN3/Mix_vs_Fiber.merged_dbcan.csv",
            "out": r"Data/dbCAN3/Mix_vs_Fibre_supplementary_figure_xx.png",
            "title": "HMO/Fibre vs Fibre",
            "xlabel": r"Growth on Fibre  $\leftarrow$  log$_2$(FC)  $\rightarrow$  Growth on HMO/Fibre mix",
        },
        {
            "expr": r"Data/dbCAN3/Mix_vs_HMOs.merged_dbcan.csv",
            "out": r"Data/dbCAN3/Mix_vs_HMOs_supplementary_figure_xx.png",
            "title": "HMO/Fibre vs HMOs",
            "xlabel": r"Growth on HMOs  $\leftarrow$  log$_2$(FC)  $\rightarrow$  Growth on HMO/Fibre mix",
        },
    ]

    args = parse_args()

    # pick species colors file once
    if args.species_colors:
        species_colors_path = Path(args.species_colors)
    else:
        default_colors = Path(r"Data/dbCAN3/species.color_legend.txt")
        species_colors_path = default_colors if default_colors.exists() else None

    # loop over all comparisons
    for cfg in PLOTS:
        expr = load_expression(Path(cfg["expr"]))

        volcano_plot(expr, Path(cfg["out"]), args.lfc, args.fdr, args.label_n, cfg["title"], cfg["xlabel"], species_colors_path)

        print(f"Rows: {len(expr)}")
        cazy_series = expr["cazy"] if "cazy" in expr.columns else build_cazy_from_columns(expr)
        print(f"CAZy-annotated rows: {pd.Series(cazy_series).notna().sum()}")
        print(f"Wrote: {cfg['out']}\n")



if __name__ == "__main__":
    main()
