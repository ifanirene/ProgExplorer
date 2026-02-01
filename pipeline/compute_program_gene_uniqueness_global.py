"""
@description
This component computes per-program gene uniqueness scores across all programs.
It is responsible for reading a top-genes loading CSV and generating a table
with both loading score and global uniqueness score (no family grouping).

Key features:
- Global TF-IDF-style uniqueness across all programs
- Coverage summary for programs in the input file

@dependencies
- pandas: DataFrame I/O and grouping
- numpy: numeric calculations

@examples
- python scripts/src/analysis/compute_program_gene_uniqueness_global.py \
    --top-genes data/Topic_analysis_k100_new/FB_moi15_seq2_loading_gene_k100_top300.csv \
    --output data/Topic_analysis_k100_new/FB_moi15_seq2_loading_gene_k100_top300_with_global_uniqueness.csv
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Tuple

import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compute global TF-IDF-style uniqueness scores for top-K loading "
            "genes across all programs."
        )
    )
    parser.add_argument(
        "--top-genes",
        type=Path,
        required=True,
        help="Top-K loading genes CSV with Name, Score, RowID or program_id.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Output CSV path with global uniqueness scores.",
    )
    return parser.parse_args()


def load_top_genes(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    if "RowID" in df.columns:
        id_col = "RowID"
    elif "program_id" in df.columns:
        id_col = "program_id"
    else:
        raise ValueError("Top-genes file missing RowID or program_id column.")

    required = {"Name", "Score", id_col}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Top-genes file missing columns: {sorted(missing)}")

    df = df.copy()
    df["Score"] = pd.to_numeric(df["Score"], errors="coerce")
    df[id_col] = pd.to_numeric(df[id_col], errors="coerce").astype("Int64")
    df = df.rename(columns={id_col: "program_id"})
    return df.dropna(subset=["Name", "Score", "program_id"])


def compute_global_uniqueness(top_genes: pd.DataFrame) -> Tuple[pd.DataFrame, int]:
    total_programs = top_genes["program_id"].nunique()
    gene_counts = top_genes.groupby("Name")["program_id"].nunique().astype(float)
    idf = np.log((total_programs + 1.0) / (gene_counts + 1.0))

    scored = top_genes.copy()
    scored["UniquenessScore"] = scored["Score"] * scored["Name"].map(idf)
    return scored, total_programs


def summarize_coverage(top_genes: pd.DataFrame, total_programs: int) -> None:
    print(f"Total programs in top-genes file: {total_programs}")
    print(f"Total genes in top-genes file: {len(top_genes)}")


def main() -> None:
    args = parse_args()
    top_genes = load_top_genes(args.top_genes)
    scored, total_programs = compute_global_uniqueness(top_genes)

    summarize_coverage(top_genes, total_programs)

    output_df = scored[
        ["Name", "Score", "program_id", "UniquenessScore"]
    ].copy()
    output_df.to_csv(args.output, index=False)
    print(f"Wrote global uniqueness scores to: {args.output}")


if __name__ == "__main__":
    main()
