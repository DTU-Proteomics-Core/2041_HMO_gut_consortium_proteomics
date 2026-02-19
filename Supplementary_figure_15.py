#!/usr/bin/env python
"""
Create a 3-panel volcano plot (A4) from pre-merged dbCAN files.
Only GH/CE/PL families are shown (GT excluded).
Title is shown inside each panel enwrapped in a box,
and the x-axis label "log2(FC)" is only shown on the bottom panel.
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
    if not color or not isinstance(color, str):
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
    return rf"$\it{{{species.replace(' ', r'\ ')}}}$"


def prepare_volcano_df(df: pd.DataFrame, lfc_cut: float, fdr_cut: float) -> pd.DataFrame:
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
    return df


def build_legend_handles(
    species_colors: dict[str, str],
    species_order: list[str],
    include_non_cazy: bool,
) -> tuple[list[Line2D], list[str]]:
    legend_handles: list[Line2D] = []
    legend_labels: list[str] = []
    species_markers = {
        "Bifidobacterium longum": ">",
        "Bifidobacterium adolescentis": "d",
        "Roseburia inulinivorans": "p",
        "Roseburia intestinalis": "p",
        "Roseburia hominis": "p",
    }
    preferred_order = [
        "Bifidobacterium longum",
        "Bifidobacterium adolescentis",
        "Roseburia inulinivorans",
        "Roseburia intestinalis",
        "Roseburia hominis",
    ]
    ordered_species = [s for s in preferred_order if s in species_colors]
    ordered_species.extend([s for s in species_order if s not in ordered_species])

    for species in ordered_species:
        color = species_colors.get(species)
        if not color:
            continue
        marker = species_markers.get(species, "o")
        legend_handles.append(
            Line2D([0], [0], marker=marker, color="none",
                   markerfacecolor=color, markeredgecolor=color, markersize=6)
        )
        legend_labels.append(format_species_label(species))

    if include_non_cazy:
        legend_handles.append(
            Line2D([0], [0], marker="o", color="none",
                   markerfacecolor="none", markeredgecolor="#9e9e9e",
                   markeredgewidth=1.5, markersize=6)
        )
        legend_labels.append("Non-CAZy")

    return legend_handles, legend_labels



def volcano_plot_ax(
    ax: plt.Axes,
    df: pd.DataFrame,
    lfc_cut: float,
    fdr_cut: float,
    label_n: int,
    xlabel: str,        
    direction_text: str,  
    species_colors: dict[str, str],
    species_order: list[str],
    use_species_colors: bool,
    panel_label: str,
    y_max: float,
    box_x: float = 0.5, 
    box_y: float = 0.035,
) -> None:
    ax.scatter(
        df["logfc"], df["neglog10_fdr"],
        s=15, facecolors="none", edgecolors="#9e9e9e",
        alpha=0.5, linewidths=0.6,
    )

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
        ordered_species = [s for s in preferred_order if s in species_colors]
        ordered_species.extend([s for s in species_order if s not in ordered_species])
        for species in ordered_species:
            color = species_colors.get(species)
            if not color:
                continue
            marker = species_markers.get(species, "o")
            subset = df[df["has_cazy"] & (df["strain_species"] == species)]
            if subset.empty:
                continue
            ax.scatter(
                subset["logfc"], subset["neglog10_fdr"],
                s=30, c=color, alpha=0.9, linewidths=0, marker=marker,
            )

    ax.axvline(-lfc_cut, color="#666666", linestyle="--", linewidth=1)
    ax.axvline(lfc_cut, color="#666666", linestyle="--", linewidth=1)
    ax.axhline(-np.log10(fdr_cut), color="#666666", linestyle="--", linewidth=1)

    ax.set_xlim(-10, 10)
    ax.set_ylim(0, y_max*1.01)

    x_min, x_max = -10, 10
    lfc_label = 0.585
    x_offset = (x_max - x_min) * 0.006


    ax.text(
        -lfc_cut - x_offset, y_max * 0.035,
        f"log$_2$(FC) = {-lfc_label:.3f}",
        ha="right", va="bottom", fontsize=8, fontweight="bold",
    )
    ax.text(
        lfc_cut + x_offset, y_max * 0.035,
        f"log$_2$(FC) = {lfc_label:.3f}",
        ha="left", va="bottom", fontsize=8, fontweight="bold",
    )
    ax.text(
        x_min + (x_max - x_min) * 0.02, -np.log10(fdr_cut),
        f"FDR = {fdr_cut:.2f}",
        ha="left", va="bottom", fontsize=8, fontweight="bold",
    )

    df["score"] = df["logfc"].abs() * df["neglog10_fdr"]
    cazy_str = df["cazy"].astype(str)
    non_gt = ~cazy_str.str.startswith("GT", na=False)
    top = df[df["has_cazy"] & non_gt].sort_values("score", ascending=False).head(label_n)

    gh_mask = df["has_cazy"] & cazy_str.str.contains(
        r"\bGH10\b|\bGH26\b|\bGH29\b|\bGH33\b|\bGH95\b|\bGH112\b|\bGH109\b|\bGH49\b|\bGH2\b",
        regex=True, na=False,
    )
    gh_targets = df[gh_mask]
    rose_inu_targets = df[df["has_cazy"] & (df["strain_species"] == "Roseburia inulinivorans")]

    label_df = pd.concat([top, gh_targets, rose_inu_targets], ignore_index=True).drop_duplicates()

    texts = []
    label_points = []
    target_set = {"GH10", "GH26", "GH29", "GH33", "GH95", "GH112", "GH109", "GH49", "GH2"}

    for _, row in label_df.iterrows():
        tokens = str(row.get("cazy") or "").split("+")
        targets = [t for t in tokens if t in target_set]
        annotation_tokens = [t for t in tokens if t.startswith(("GH", "PL", "CE"))]
        label = "+".join(targets) if targets else (row.get("cazy") or row.get("protein_accession"))
        if pd.isna(label):
            continue

        x = row["logfc"]
        y = row["neglog10_fdr"]

        label_color = "#333333"
        label_weight = "bold" if annotation_tokens else "normal"
        if annotation_tokens and use_species_colors:
            species = row.get("strain_species")
            if isinstance(species, str) and species in species_colors:
                label_color = darken_hex(species_colors[species], factor=0.65)

        texts.append(
            ax.text(
                x, y, str(label),
                fontsize=7, zorder=5,
                color=label_color, fontweight=label_weight,
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="none", alpha=0.5),
            )
        )
        label_points.append((x, y))

    if texts:
        adjust_text(texts, ax=ax)
        for text, (x, y) in zip(texts, label_points):
            ax.annotate(
                "",
                xy=(x, y),
                xytext=text.get_position(),
                textcoords="data",
                arrowprops=dict(arrowstyle="-", color="#888888", lw=0.6, shrinkA=6),
                zorder=4,
            )

    ax.set_xlabel(xlabel, fontsize=10, fontweight="bold")
    ax.set_ylabel(r"log$_{10}$(FDR)", fontsize=10, fontweight="bold")
    ax.tick_params(axis="both", labelsize=9)

 
    ax.text(
        0.01, 0.98, panel_label,
        transform=ax.transAxes,
        ha="left", va="top",
        fontsize=13, fontweight="bold",
    )
    ax.text( 
        box_x, box_y,
        direction_text,
        transform=ax.transAxes,
        ha="center", va="bottom",
        fontsize=10, fontweight="bold",
        bbox=dict(boxstyle="square,pad=0.25", facecolor="white", edgecolor="black",linewidth=0.8),
        zorder=10,
    )



def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Plot volcano from merged dbCAN data.")
    p.add_argument("--lfc", type=float, default=0.585, help="log2FC cutoff")
    p.add_argument("--fdr", type=float, default=0.05, help="FDR cutoff")
    p.add_argument("--label-n", type=int, default=50, help="Number of CAZy hits to label")
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
            "direction": r"    Fibre  $\leftarrow$  Substrate  $\rightarrow$  HMOs     ",
            "panel": "a",
            "box_x":0.50, "box_y":0.93
        },
        {
            "expr": r"Data/dbCAN3/Mix_vs_Fiber.merged_dbcan.csv",
            "direction": r"       Fibre  $\leftarrow$  Substrate  $\rightarrow$  HMO/Fibre ",
            "panel": "b",
            "box_x":0.50, "box_y":0.925
        },
        {
            "expr": r"Data/dbCAN3/Mix_vs_HMOs.merged_dbcan.csv",
            "direction": r"      HMOs  $\leftarrow$  Substrate  $\rightarrow$  HMO/Fibre ",
            "panel": "c",
            "box_x":0.50, "box_y":0.905
        },
    ]

    args = parse_args()

    if args.species_colors:
        species_colors_path = Path(args.species_colors)
    else:
        default_colors = Path(r"Data/dbCAN3/species.color_legend.txt")
        species_colors_path = default_colors if default_colors.exists() else None

    prepared: list[tuple[dict, pd.DataFrame]] = []
    y_max_raw: list[float] = []
    for cfg in PLOTS:
        expr = load_expression(Path(cfg["expr"]))
        df_prepped = prepare_volcano_df(expr, args.lfc, args.fdr)
        prepared.append((cfg, df_prepped))

        m = df_prepped["neglog10_fdr"].max()
        m = max(m, -np.log10(args.fdr))
        y_max_raw.append(m if np.isfinite(m) else 1.0)

    if species_colors_path is None:
        raise FileNotFoundError(
            "No species color file found. Provide --species-colors or place species.color_legend.txt in the folder."
        )
    species_colors, species_order = load_species_colors(species_colors_path)
    use_species_colors = True
    include_non_cazy = True

    legend_handles, legend_labels = build_legend_handles(species_colors, species_order, include_non_cazy)

    height_ratios = [max(v, 0.1) for v in y_max_raw]
    y_max_scaled = [v * 1.06 for v in y_max_raw]  

    fig, axes = plt.subplots(
        3, 1,
        figsize=(8.27, 11.69),
        sharex=True,
        gridspec_kw={"height_ratios": [1.7, 1.6, 1.2]},
    )

    for i, (ax, (cfg, df_prepped)) in enumerate(zip(axes, prepared)):
        xlabel = r"log$_2$(FC)" if i == len(prepared) - 1 else ""
        volcano_plot_ax(
            ax=ax,
            df=df_prepped,
            lfc_cut=args.lfc,
            fdr_cut=args.fdr,
            label_n=args.label_n,
            xlabel=xlabel,
            direction_text=cfg["direction"],
            box_x=cfg.get("box_x", 0.5),
            box_y=cfg.get("box_y", 0.95),
            species_colors=species_colors,
            species_order=species_order,
            use_species_colors=use_species_colors,
            panel_label=cfg["panel"],
            y_max=y_max_scaled[i],
        )
        if i < 2:  
            ax.tick_params(axis="x", labelbottom=False)

        print(f"Rows: {len(df_prepped)}")
        print(f"CAZy-annotated rows: {df_prepped['has_cazy'].sum()}")

    if legend_handles:
        fig.legend(
            legend_handles,
            legend_labels,
            frameon=False,
            loc="lower center",
            bbox_to_anchor=(0.5, 0.01),
            ncol=3,
            columnspacing=1.4,
            handletextpad=0.6,
        )

    fig.tight_layout(rect=[0, 0.045, 1, 1])
    fig.subplots_adjust(hspace=0.04)

    outpath = Path("Data/dbCAN3/Supplementary_figure_15.png")
    fig.savefig(outpath, dpi=1000)
    plt.close(fig)
    print(f"Wrote: {outpath}")


if __name__ == "__main__":
    main()