#!/usr/bin/env python

from __future__ import annotations
import argparse
import re
from pathlib import Path
from typing import cast
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.lines import Line2D
from adjustText import adjust_text

plt.style.use("seaborn-v0_8-whitegrid")
plt.rcParams["mathtext.fontset"] = "dejavusans"
plt.rcParams["font.family"] = "DejaVu Sans"

UNIPROT_RE = re.compile(r"\b(?:sp|tr)?\|?([A-NR-Z0-9]{5,10})\|?")

HIGHLIGHT_GH_MARKERS = {
    #HMO-degrading GHs
    "GH112": "d",
    "GH136": "d",
    "GH20": "d",
    "GH29": "d",
    "GH95": "d",
    "GH33": "d",
    "GH181": "d",
    #Arabinoxylan-degrading GHs
    "GH43": "*",
    "GH10": "*",
    "GH74": "*",
    "CE17": "*",
    #Mannan-degrading GHs
    "GH130": "s",
    "GH26": "s",
    "CE2": "s",

}

GH_REGEX = r"\b(?:" + "|".join(HIGHLIGHT_GH_MARKERS.keys()) + r")\b"

SPECIES_COLORS = {
    "Bifidobacterium longum": "#ffd92f",
    "Bifidobacterium adolescentis": "#fc7465",
    "Roseburia inulinivorans": "#60a4d4",
    "Roseburia intestinalis": "#70d396",
    "Roseburia hominis": "#a299d9",
}

LINE_GREY = "#bdbdbd"

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


def load_expression(path: Path, sheet: str | int | None = 0) -> pd.DataFrame:
    if path.suffix.lower() in {".xlsx", ".xls"}:
        sheet_name = 0 if sheet is None else sheet
        df = pd.read_excel(path, sheet_name=sheet_name)
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
    acc_series = cast(pd.Series, df[acc_col])
    df["protein_accession_raw"] = acc_series.astype(str)
    df["protein_accession"] = cast(pd.Series, df["protein_accession_raw"]).map(extract_accession)

    if "logfc" not in df.columns:
        for col in df.columns:
            if re.search(r"log|fold", col):
                df["logfc"] = pd.to_numeric(cast(pd.Series, df[col]), errors="coerce")
                break
    if "fdr" not in df.columns:
        for col in df.columns:
            if re.search(r"fdr|padj|qvalue|p_adj|p.adj", col):
                df["fdr"] = pd.to_numeric(cast(pd.Series, df[col]), errors="coerce")
                break

    return df


def build_cazy_from_columns(df: pd.DataFrame) -> pd.Series:
    if not any(col in df.columns for col in ("gh", "ce", "pl")):
        cazy = df.get("cazy")
        if isinstance(cazy, pd.Series):
            return cazy
        return pd.Series([None] * len(df), index=df.index, dtype="string")

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
        return r"$\it{B.\ infantis}$"
    if species == "Bifidobacterium adolescentis":
        return r"$\it{B.\ adolescentis}$"
    if species == "Roseburia inulinivorans":
        return r"$\it{R.\ inulinivorans}$"
    if species == "Roseburia intestinalis":
        return r"$\it{R.\ intestinalis}$"
    if species == "Roseburia hominis":
        return r"$\it{R.\ hominis}$"
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

def build_species_legend_handles(species_colors):
    handles, labels = [], []

    def blank():
        return Line2D([], [], linestyle="None", marker=None, color="none")

    def dot(sp):
        return Line2D([0], [0], marker="o", linestyle="None",
                      color="none",
                      markerfacecolor=species_colors[sp],
                      markeredgecolor=species_colors[sp],
                      markersize=7)


    handles += [blank(), dot("Bifidobacterium longum"), blank(), blank()]
    labels  += ["HMO utiliser",
                format_species_label("Bifidobacterium longum"),
                "", ""]

    handles += [blank(), dot("Bifidobacterium adolescentis"), blank(), blank()]
    labels  += ["Fibre utiliser",
                format_species_label("Bifidobacterium adolescentis"),
                "", ""]

    handles += [
        blank(),
        dot("Roseburia intestinalis"),
        dot("Roseburia hominis"),
        dot("Roseburia inulinivorans"),
    ]
    labels += [
        "Dual HMO-fibre utiliser",
        format_species_label("Roseburia intestinalis"),
        format_species_label("Roseburia hominis"),
        format_species_label("Roseburia inulinivorans"),
    ]

    return handles, labels


def build_marker_legend_handles() -> tuple[list[Line2D], list[str]]:
    blank = Line2D([], [], linestyle="None", marker=None, color="none")

    handles = [
        # Column 1
        Line2D([0], [0], marker="d",  color="none", 
               markerfacecolor="white", markeredgecolor="black",
               markersize=6, linestyle="None"),
        Line2D([0], [0], marker="*", color="none",
               markerfacecolor="white", markeredgecolor="black",
               markersize=9, linestyle="None"),
        Line2D([0], [0], marker="s", color="none",
               markerfacecolor="white", markeredgecolor="black",
               markersize=6, linestyle="None"),

        # Column 2
        Line2D([0], [0], marker="x", color="none",
               markerfacecolor="white", markeredgecolor="black",
               markersize=6, linestyle="None"),
        Line2D([0], [0], marker="o", color="none",
               markerfacecolor="white", markeredgecolor=LINE_GREY,
               markersize=6, linestyle="None"),
        blank,
    ]

    labels = [
        "HMOs",
        "Arabinoxylan",
        "Galactomannan",
        "Others",
        "Non-CAZymes",
        "",
    ]

    return handles, labels


def volcano_plot_ax(
    ax: Axes,
    df: pd.DataFrame,
    lfc_cut: float,
    fdr_cut: float,
    label_n: int,
    xlabel: str,  
    ylabel: str,      
    direction_text: str,  
    species_colors: dict[str, str],
    species_order: list[str],
    use_species_colors: bool,
    panel_label: str,
    y_max: float,
    box_x: float = 0.5, 
    box_y: float = 0.035,
) -> None:
    

    non_cazy_df = df[~df["has_cazy"]]

    ax.scatter(
        non_cazy_df["logfc"], non_cazy_df["neglog10_fdr"],
        s=15,
        facecolors="none",
        edgecolors="#9e9e9e",
        alpha=0.4,
        linewidths=0.6,
        marker="o",
        zorder=1,
    )

    if species_colors:
        preferred_order = [
            "Bifidobacterium longum",
            "Bifidobacterium adolescentis",
            "Roseburia intestinalis",
            "Roseburia hominis",
            "Roseburia inulinivorans",
        ]
        ordered_species = [s for s in preferred_order if s in species_colors]
        ordered_species.extend([s for s in species_order if s not in ordered_species])
        for species in ordered_species:
            color = species_colors.get(species)
            if not color:
                continue

            species_df = df[df["has_cazy"] & (df["strain_species"] == species)]
            if species_df.empty:
                continue

            normal = species_df[~species_df["cazy"].astype(str).str.contains(GH_REGEX, na=False)]

            if not normal.empty:
                ax.scatter(
                    normal["logfc"], normal["neglog10_fdr"],
                    s=35,
                    c=color,
                    alpha=0.7,
                    marker="x",
                    zorder=2,
                )

            highlight = species_df[species_df["cazy"].astype(str).str.contains(GH_REGEX, na=False)]

            for gh, marker in HIGHLIGHT_GH_MARKERS.items():
                sub = highlight[highlight["cazy"].astype(str).str.contains(rf"\b{gh}\b", na=False)]
                if sub.empty:
                    continue

                size_map = {
                    "*": 95,  
                    "s": 40,  
                    "d": 55,  
                }

                ax.scatter(
                    sub["logfc"], sub["neglog10_fdr"],
                    s=size_map.get(marker, 45),
                    facecolors=color,
                    alpha=0.9,
                    edgecolors="black",
                    linewidths=0.7,
                    marker=marker,
                    zorder=3,
                )

    ax.axvline(-lfc_cut, color=LINE_GREY, linestyle="--", alpha=0.8, linewidth=1)
    ax.axvline(lfc_cut, color=LINE_GREY, linestyle="--", alpha=0.8, linewidth=1)
    ax.axhline(-np.log10(fdr_cut), color=LINE_GREY, linestyle="--", alpha=0.8, linewidth=1)

    ax.grid(False)
    ax.set_axisbelow(True)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.set_xlim(-10, 10)
    ax.set_ylim(0, y_max*1.01)

    plt.rcParams.update({
        "font.size": 9,
        "axes.labelsize": 9,    
        "axes.titlesize": 9,})


    x_min, x_max = -10, 10
    lfc_label = 0.585
    x_offset = (x_max - x_min) * 0.007

    ax.text(
        -lfc_cut - x_offset, 0.15,  
        f"–{lfc_label:.3f}",
        ha="right", va="bottom",
        fontsize=7, color="#6f6f6f",
        bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="none", alpha=0.7),
    )
    ax.text(
        lfc_cut + x_offset, 0.15,
        f"{lfc_label:.3f}",
        ha="left", va="bottom",
        fontsize=7, color="#6f6f6f",
        bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="none", alpha=0.7),
    )

    ax.text(
        x_min + (x_max - x_min) * 0.02,
        -np.log10(fdr_cut),
        rf"$P$" + f".adj = {fdr_cut:.2f}",
        ha="left", va="bottom",
        fontsize=7,
        color="#6f6f6f",
    )

    df = df.copy()
    df["score"] = df["logfc"].abs() * df["neglog10_fdr"]
    cazy_str = df["cazy"].astype(str)
    non_gt = ~cazy_str.str.startswith("GT", na=False)

    sig_cazy = df[
        df["has_cazy"] & non_gt &
        df["significant"] &
        cazy_str.str.contains(r"\b(?:GH|CE|PL)\d+\b", regex=True, na=False)
    ]
    if label_n and len(sig_cazy) > label_n:
        sig_cazy = sig_cazy.nlargest(label_n, "score")

    highlighted = df[
        df["has_cazy"] & non_gt &
        cazy_str.str.contains(GH_REGEX, regex=True, na=False)
    ]

    label_df = pd.concat([sig_cazy, highlighted]).drop_duplicates()

        
    texts = []
    label_points = []

    for _, row in label_df.iterrows():

        tokens = str(row.get("cazy") or "").split("+")

        annotation_tokens = [
            t for t in tokens
            if t.startswith(("GH", "CE", "PL"))
        ]

        if not annotation_tokens:
            continue

        label = "+".join(annotation_tokens)
        if pd.isna(label):
            continue

        x = row["logfc"]
        y = row["neglog10_fdr"]

        species = row.get("strain_species")

        label_color = "#333333"
        if annotation_tokens and use_species_colors:
            species = row.get("strain_species")
            label_color = "#333333"
            if isinstance(species, str) and species in species_colors:
                label_color = darken_hex(species_colors[species], factor=0.75)

        texts.append(
            ax.text(
                x, y, str(label),
                fontsize=7, zorder=5,
                color=label_color, fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.15", facecolor="white", edgecolor = "none", alpha=0.7),
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
                arrowprops=dict(arrowstyle="-", color="black", lw=0.6, shrinkA=6),
                zorder=4,
            )

    ax.set_xlabel(xlabel, fontsize=10)
    ax.set_ylabel(ylabel, fontsize=10)
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
        fontsize=9, fontweight="bold",
        bbox=dict(boxstyle="square,pad=0.25", facecolor="white", edgecolor=LINE_GREY, linewidth=0.8),
        zorder=10,
    )


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Plot volcano from merged dbCAN data.")
    p.add_argument("--lfc", type=float, default=0.585, help="log2FC cutoff")
    p.add_argument("--fdr", type=float, default=0.05, help="FDR cutoff")
    p.add_argument("--label-n", type=int, default=15, help="Number of CAZy hits to label")
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

    species_colors = SPECIES_COLORS
    species_order = list(SPECIES_COLORS.keys())

    prepared: list[tuple[dict, pd.DataFrame]] = []
    y_max_raw: list[float] = []
    for cfg in PLOTS:
        expr = load_expression(Path(cfg["expr"]))
        df_prepped = prepare_volcano_df(expr, args.lfc, args.fdr)
        prepared.append((cfg, df_prepped))

        m = df_prepped["neglog10_fdr"].max()
        m = max(m, -np.log10(args.fdr))
        y_max_raw.append(m if np.isfinite(m) else 1.0)

    use_species_colors = True

    species_handles, species_labels = build_species_legend_handles(
        species_colors,
    )

    marker_handles, marker_labels = build_marker_legend_handles()

    y_max_scaled = [v * 1.06 for v in y_max_raw]  

    fig, axes = plt.subplots(
        3, 1,
        #figsize=(8.27, 10.69),
        figsize=(8.27, 9.5),
        sharex=True,
        gridspec_kw={"height_ratios": [1.7, 1.6, 1.2]},
    )

    for i, (ax, (cfg, df_prepped)) in enumerate(zip(axes, prepared)):
        xlabel = r"log$_2$(Fold Change)" if i == len(prepared) - 1 else ""
        ylabel = r"–log$_{10}$($\it{P}$ adjusted)" if i == 1 else ""

        volcano_plot_ax(
            ax=ax,
            df=df_prepped,
            lfc_cut=args.lfc,
            fdr_cut=args.fdr,
            label_n=args.label_n,
            xlabel=xlabel,
            ylabel=ylabel,
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

    species_legend = fig.legend(
        species_handles,
        species_labels,
        frameon=False,
        loc="lower left",
        bbox_to_anchor=(0.06, 0.001),
        ncol=3,
        fontsize=9,
        handletextpad=0.5,
        columnspacing=1.8,
    )

    fig.legend(
        marker_handles,
        marker_labels,
        title="CAZyme specificity",
        frameon=False,
        loc="lower left",
        bbox_to_anchor=(0.62, 0.001),
        ncol=2,
        fontsize=9,
        title_fontproperties={"weight": "bold", "size": 9},
        handletextpad=0.6,
        columnspacing=1.2,
    )

    for i in [0, 4, 8]:
        species_legend.get_texts()[i].set_fontweight("bold")

    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.subplots_adjust(hspace=0.04)

    outpath = Path("Data/dbCAN3/Figure_6.png")
    fig.savefig(outpath, dpi=1000)
    plt.close(fig)
    print(f"Wrote: {outpath}")


if __name__ == "__main__":
    main()