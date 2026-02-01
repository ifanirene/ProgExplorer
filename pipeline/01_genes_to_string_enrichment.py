"""
/**
 * @description 
 * Unified CLI to create program-wise top-N gene lists from a cNMF loading CSV and
 * run STRING functional enrichment per program with optional figure generation.
 * 
 * It merges functionality from `annotate_cnmf_programs_string.py` and
 * `run_string_enrichment.py` into subcommands:
 * - extract:  read loading CSV → save JSON {program_id: [genes]} and overview CSV
 * - enrich:   read JSON → call STRING API → write full and filtered CSVs, optionally figures
 * - all:      extract → enrich (convenience)
 * 
 * Key features:
 * - Robust HTTP with retries and pacing
 * - Full unfiltered CSV and Process/KEGG filtered CSV (<500 background genes)
 * - Direct retrieval of enrichment figures from STRING API
 * - Writes a global UniquenessScore gene table for downstream steps
 * 
 * @dependencies
 * - pandas, requests
 * 
 * @examples
 * - Extract top 100 genes per program (RowID or program_id column):
 *   python pipeline/01_genes_to_string_enrichment.py extract \
 *     --input input/genes/FB_moi15_seq2_loading_gene_k100_top300.csv \
 *     --n-top 100 \
 *     --json-out input/enrichment/genes_top100.json \
 *     --csv-out input/enrichment/genes_overview_top100.csv
 * 
 * - Run enrichment and figures:
 *   python pipeline/01_genes_to_string_enrichment.py enrich \
 *     --genes-json input/enrichment/genes_top100.json \
 *     --species 10090 \
 *     --out-csv-full input/enrichment/string_enrichment_full.csv \
 *     --out-csv-filtered input/enrichment/string_enrichment_filtered_process_kegg.csv \
 *     --figures-dir input/enrichment/enrichment_figures
 * 
 * - End-to-end:
 *   python pipeline/01_genes_to_string_enrichment.py all \
 *     --input input/genes/FB_moi15_seq2_loading_gene_k100_top300.csv \
 *     --n-top 100 \
 *     --json-out input/enrichment/genes_top100.json \
 *     --csv-out input/enrichment/genes_overview_top100.csv \
 *     --species 10090 \
 *     --out-csv-full input/enrichment/string_enrichment_full.csv \
 *     --out-csv-filtered input/enrichment/string_enrichment_filtered_process_kegg.csv
 */
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Set

import numpy as np
import pandas as pd
import requests


logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

"""
@description
Configuration loader for the topic annotation workflow.
It is responsible for reading JSON/YAML configs and returning a dict for
per-step defaults with CLI override support.

Key features:
- Supports JSON and YAML (if PyYAML is installed).
- Returns an empty dict when no config is provided.

@dependencies
- json: Built-in JSON parser
- yaml (optional): YAML parser when available

@examples
- cfg = load_config("configs/example_config.yaml")
"""


def load_config(config_path: Optional[str]) -> Dict[str, Any]:
    if not config_path:
        return {}
    path = Path(config_path)
    if not path.exists():
        raise SystemExit(f"Config file not found: {path}")

    suffix = path.suffix.lower()
    if suffix in {".yaml", ".yml"}:
        try:
            import yaml  # type: ignore
        except ImportError as exc:
            raise SystemExit("PyYAML is required for YAML configs.") from exc
        data = yaml.safe_load(path.read_text(encoding="utf-8"))
    else:
        data = json.loads(path.read_text(encoding="utf-8"))

    if not isinstance(data, dict):
        raise SystemExit("Config must be a mapping at the top level.")
    return data


"""
@description
Utility helpers for merging CLI arguments with config defaults.
It is responsible for identifying CLI overrides and applying test-mode filters.

Key features:
- Detects explicitly provided CLI flags for override precedence.
- Applies optional test-mode program filters when configured.

@dependencies
- sys: Access to raw CLI arguments
"""


def get_cli_overrides(argv: List[str]) -> Set[str]:
    overrides: Set[str] = set()
    for token in argv:
        if token.startswith("--"):
            name = token[2:]
            if "=" in name:
                name = name.split("=", 1)[0]
            overrides.add(name.replace("-", "_"))
    return overrides


def parse_topics(value: Optional[object]) -> Optional[Set[int]]:
    if value is None:
        return None
    if isinstance(value, list):
        return {int(v) for v in value}
    if isinstance(value, str):
        items = [item.strip() for item in value.split(",") if item.strip()]
        return {int(v) for v in items}
    return None


def apply_test_mode(
    args: argparse.Namespace, config: Dict[str, Any], cli_overrides: Set[str]
) -> argparse.Namespace:
    test_cfg = config.get("test", {}) if isinstance(config.get("test", {}), dict) else {}
    enabled = bool(test_cfg.get("enabled") or config.get("test_mode"))
    if not enabled:
        return args

    topics = test_cfg.get("topics") or test_cfg.get("programs") or config.get("test_programs")
    if hasattr(args, "topics") and "topics" not in cli_overrides and not getattr(args, "topics", None):
        if topics is not None:
            args.topics = topics
    return args


def apply_config_overrides(
    args: argparse.Namespace,
    config: Dict[str, Any],
    cli_overrides: Set[str],
) -> argparse.Namespace:
    steps_cfg = config.get("steps", {}) if isinstance(config.get("steps", {}), dict) else {}
    step_cfg = steps_cfg.get("string_enrichment", {})
    if isinstance(step_cfg, dict) and args.command in step_cfg:
        step_cfg = step_cfg.get(args.command, {})
    if not isinstance(step_cfg, dict):
        return args

    for key, value in step_cfg.items():
        dest = key.replace("-", "_")
        if dest in cli_overrides:
            continue
        if hasattr(args, dest):
            setattr(args, dest, value)
    return args


DEFAULT_ENRICH_DIR = Path("input/enrichment")
DEFAULT_GENES_JSON_TEMPLATE = "genes_top{n_top}.json"
DEFAULT_GENES_OVERVIEW_TEMPLATE = "genes_overview_top{n_top}.csv"
DEFAULT_STRING_FULL = "string_enrichment_full.csv"
DEFAULT_STRING_FILTERED = "string_enrichment_filtered_process_kegg.csv"

"""
@description
Default output path resolver for the STRING enrichment CLI.
It is responsible for filling in output paths when CLI args are omitted
so the pipeline can run from a single input CSV.

Key features:
- Defaults enrichment outputs to input/enrichment/
- Uses n_top to name gene list outputs
- Honors explicit CLI/config values

@dependencies
- pathlib: Path composition for defaults
"""


def apply_default_paths(args: argparse.Namespace) -> argparse.Namespace:
    n_top = getattr(args, "n_top", None) or 100

    if hasattr(args, "json_out") and not getattr(args, "json_out", None):
        args.json_out = str(DEFAULT_ENRICH_DIR / DEFAULT_GENES_JSON_TEMPLATE.format(n_top=n_top))
    if hasattr(args, "csv_out") and not getattr(args, "csv_out", None):
        args.csv_out = str(DEFAULT_ENRICH_DIR / DEFAULT_GENES_OVERVIEW_TEMPLATE.format(n_top=n_top))
    if hasattr(args, "out_csv_full") and not getattr(args, "out_csv_full", None):
        args.out_csv_full = str(DEFAULT_ENRICH_DIR / DEFAULT_STRING_FULL)
    if hasattr(args, "out_csv_filtered") and not getattr(args, "out_csv_filtered", None):
        args.out_csv_filtered = str(DEFAULT_ENRICH_DIR / DEFAULT_STRING_FILTERED)

    return args


"""
@description
Helpers for caching STRING enrichment results per program.
They are responsible for reading and writing cached JSON payloads to avoid
re-querying STRING for previously processed programs.

Key features:
- Simple per-program JSON cache files
- Graceful fallback when cache is missing or invalid

@dependencies
- json: read/write cached enrichment payloads
- pathlib: cache path handling
"""


def cache_path(cache_dir: Path, program_id: int) -> Path:
    return cache_dir / f"program_{program_id}_enrichment.json"


def load_cached_results(cache_dir: Path, program_id: int) -> Optional[List[Dict[str, Any]]]:
    cache_file = cache_path(cache_dir, program_id)
    if not cache_file.exists():
        return None
    try:
        data = json.loads(cache_file.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        logger.warning("Cache file is invalid JSON: %s", cache_file)
        return None
    if not isinstance(data, list):
        logger.warning("Cache file has unexpected format: %s", cache_file)
        return None
    return data


def write_cached_results(cache_dir: Path, program_id: int, results: List[Dict[str, Any]]) -> None:
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_file = cache_path(cache_dir, program_id)
    cache_file.write_text(json.dumps(results, indent=2), encoding="utf-8")


# --------------------------- Extract top genes (CSV) --------------------------

def ensure_parent_dir(path_str: str) -> None:
    path = Path(path_str)
    if path.parent and not path.parent.exists():
        path.parent.mkdir(parents=True, exist_ok=True)


def resolve_program_id_column(df: pd.DataFrame) -> str:
    if "RowID" in df.columns:
        return "RowID"
    if "program_id" in df.columns:
        return "program_id"
    raise ValueError("Input is missing required column: RowID or program_id")


def normalize_program_id(value: object) -> object:
    try:
        return int(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return value


def extract_top_genes_by_program(
    df: pd.DataFrame, n_top: int, id_col: str
) -> Dict[str, List[str]]:
    required_cols = {"Name", "Score", id_col}
    missing = required_cols.difference(df.columns)
    if missing:
        raise ValueError(f"Input is missing required columns: {sorted(missing)}")

    top_map: Dict[str, List[str]] = {}
    for program_id, sub in df.groupby(id_col, sort=True):
        program_id_norm = normalize_program_id(program_id)
        program_key = str(program_id_norm)
        sub_sorted = sub.sort_values("Score", ascending=False).head(n_top)
        genes = [str(g) for g in sub_sorted["Name"].dropna().tolist()]
        seen = set()
        unique_genes: List[str] = []
        for g in genes:
            if g not in seen:
                seen.add(g)
                unique_genes.append(g)
        top_map[program_key] = unique_genes
    return top_map


"""
@description
This component builds a gene loading table with global UniquenessScore values.
It is responsible for normalizing program identifiers, computing TF-IDF-style
uniqueness across all programs, and exporting a table for downstream steps.

Key features:
- Accepts RowID or program_id inputs.
- Preserves family_id if available.
- Computes UniquenessScore when missing.

@dependencies
- numpy: IDF calculation
- pandas: DataFrame manipulation

@examples
- uniqueness_df = build_uniqueness_table(df, id_col="RowID")
"""


def default_uniqueness_output(input_path: Path) -> Path:
    suffix = input_path.suffix or ".csv"
    stem = input_path.stem
    if stem.endswith("_with_uniqueness"):
        return input_path
    return input_path.with_name(f"{stem}_with_uniqueness{suffix}")


def build_uniqueness_table(df: pd.DataFrame, id_col: str) -> pd.DataFrame:
    required_cols = {"Name", "Score", id_col}
    missing = required_cols.difference(df.columns)
    if missing:
        raise ValueError(f"Input is missing required columns: {sorted(missing)}")

    work = df.copy()
    if "program_id" not in work.columns:
        work["program_id"] = work[id_col]

    if "UniquenessScore" not in work.columns or work["UniquenessScore"].isna().all():
        work["Score"] = pd.to_numeric(work["Score"], errors="coerce")
        work["program_id"] = pd.to_numeric(work["program_id"], errors="coerce")
        valid = work.dropna(subset=["Name", "Score", "program_id"]).copy()
        if valid.empty:
            raise ValueError("No valid rows to compute UniquenessScore.")

        valid["program_id"] = valid["program_id"].astype(int)
        total_programs = valid["program_id"].nunique()
        gene_counts = valid.groupby("Name")["program_id"].nunique().astype(float)
        idf = np.log((total_programs + 1.0) / (gene_counts + 1.0))
        valid["UniquenessScore"] = valid["Score"] * valid["Name"].map(idf)

        work["UniquenessScore"] = np.nan
        work.loc[valid.index, "UniquenessScore"] = valid["UniquenessScore"]

    columns = ["Name", "Score", "program_id"]
    if "family_id" in work.columns:
        columns.append("family_id")
    columns.append("UniquenessScore")
    out_df = work[columns].copy()
    out_df.dropna(subset=["Name", "Score", "program_id", "UniquenessScore"], inplace=True)
    return out_df


def build_overview_long_table(
    df: pd.DataFrame, top_map: Dict[str, List[str]], id_col: str
) -> pd.DataFrame:
    records = []
    sub_indexed_cache: Dict[object, pd.DataFrame] = {}
    for program_id_str, genes in top_map.items():
        program_id = normalize_program_id(program_id_str)
        if program_id not in sub_indexed_cache:
            sub_indexed_cache[program_id] = (
                df[df[id_col] == program_id].set_index("Name")
            )
        sub_idx = sub_indexed_cache[program_id]
        for rank, gene in enumerate(genes, start=1):
            score = float(sub_idx.loc[gene, "Score"]) if gene in sub_idx.index else float("nan")
            records.append(
                {"program_id": program_id, "rank": rank, "gene": gene, "score": score}
            )
    out_df = pd.DataFrame.from_records(records)
    if not out_df.empty:
        out_df.sort_values(["program_id", "rank"], inplace=True)
    return out_df


def cmd_extract(args: argparse.Namespace) -> int:
    if not args.input or not args.json_out or not args.csv_out:
        logger.error("--input, --json-out, and --csv-out are required.")
        return 2
    input_path = Path(args.input)
    if not input_path.exists():
        logger.error(f"Input not found: {input_path}")
        return 2
    df = pd.read_csv(input_path)
    logger.info(f"Loaded input: {input_path} with shape {df.shape} and columns {list(df.columns)}")

    id_col = resolve_program_id_column(df)
    uniqueness_out = args.gene_loading_out or str(default_uniqueness_output(input_path))
    uniqueness_df = build_uniqueness_table(df, id_col=id_col)
    ensure_parent_dir(uniqueness_out)
    uniqueness_df.to_csv(uniqueness_out, index=False)
    logger.info(
        "Wrote uniqueness CSV: %s (rows=%s)",
        uniqueness_out,
        uniqueness_df.shape[0],
    )

    top_map = extract_top_genes_by_program(df=df, n_top=args.n_top, id_col=id_col)
    allowed_topics = parse_topics(args.topics)
    if allowed_topics:
        top_map = {
            pid: genes for pid, genes in top_map.items() if int(pid) in allowed_topics
        }
    logger.info(f"Extracted gene lists for {len(top_map)} programs")

    ensure_parent_dir(args.json_out)
    with open(args.json_out, "w", encoding="utf-8") as f:
        json.dump(top_map, f, indent=2)
    logger.info(f"Wrote JSON: {args.json_out}")

    overview_df = build_overview_long_table(df, top_map, id_col=id_col)
    ensure_parent_dir(args.csv_out)
    overview_df.to_csv(args.csv_out, index=False)
    logger.info(f"Wrote overview CSV: {args.csv_out} (rows={overview_df.shape[0]})")
    return 0


# ----------------------------- STRING enrichment -----------------------------

STRING_ENRICH_ENDPOINT = "https://string-db.org/api/json/enrichment"


def call_string_enrichment(genes: List[str], species: int, retries: int = 3, sleep_between: float = 0.6) -> List[Dict[str, Any]]:
    identifiers_value = "\r".join(genes)
    params = {"identifiers": identifiers_value, "species": species, "caller_identity": "topic_analysis_string_enrichment"}

    attempt = 0
    while attempt <= retries:
        try:
            response = requests.get(STRING_ENRICH_ENDPOINT, params=params, timeout=60)
            if response.status_code == 200:
                try:
                    data = response.json()
                except Exception as json_err:
                    logger.error(f"Failed to parse JSON (n={len(genes)}): {json_err}")
                    data = []
                return data if isinstance(data, list) else []
            else:
                logger.warning(f"STRING returned status {response.status_code}: {response.text[:200]}")
        except requests.RequestException as e:
            logger.warning(f"HTTP error on STRING request (attempt {attempt+1}/{retries+1}): {e}")

        attempt += 1
        time.sleep(min(2.0, sleep_between * (attempt + 1)))

    return []


def build_full_csv(program_to_results: Dict[str, List[Dict[str, Any]]]) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    for pid, terms in program_to_results.items():
        for t in terms:
            rows.append(
                {
                    "program_id": int(pid),
                    "category": str(t.get("category", "NA")),
                    "term": str(t.get("term", t.get("description", "NA"))),
                    "term_id": str(t.get("term_id", "NA")),
                    "description": str(t.get("description", t.get("term", "NA"))),
                    "fdr": float(t.get("fdr", float("nan"))),
                    "p_value": float(t.get("p_value", float("nan"))),
                    "number_of_genes": int(t.get("number_of_genes", 0)),
                    "number_of_genes_in_background": int(t.get("number_of_genes_in_background", 0)),
                    "ncbiTaxonId": int(t.get("ncbiTaxonId", 0)),
                    "inputGenes": "|".join(t.get("inputGenes", [])) if t.get("inputGenes") else "",
                }
            )
    df = pd.DataFrame(rows)
    if not df.empty:
        df.sort_values(["program_id", "fdr", "p_value"], inplace=True)
    return df


def filter_process_kegg(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df
    category_mask = df["category"].str.contains("Process|KEGG", case=False, na=False)
    background_mask = df["number_of_genes_in_background"] < 500
    filtered_df = df[category_mask & background_mask].copy()
    if not filtered_df.empty:
        filtered_df.sort_values(["program_id", "fdr", "p_value"], inplace=True)
    return filtered_df


def cmd_enrich(args: argparse.Namespace) -> int:
    if not args.genes_json:
        logger.error("--genes-json is required.")
        return 2
    if args.figures_only and not args.figures_dir:
        logger.error("--figures-dir is required when --figures-only is set.")
        return 2
    if not args.figures_only and (not args.out_csv_full or not args.out_csv_filtered):
        logger.error("--out-csv-full and --out-csv-filtered are required unless --figures-only is set.")
        return 2

    genes_path = Path(args.genes_json)
    if not genes_path.exists():
        logger.error(f"Genes JSON not found: {genes_path}")
        return 2

    program_to_genes: Dict[str, List[str]] = json.loads(genes_path.read_text(encoding="utf-8"))
    allowed_topics = parse_topics(args.topics)
    if allowed_topics:
        program_to_genes = {
            pid: genes
            for pid, genes in program_to_genes.items()
            if int(pid) in allowed_topics
        }
    logger.info(f"Loaded gene lists for {len(program_to_genes)} programs from {genes_path}")

    existing_full_df: Optional[pd.DataFrame] = None
    existing_programs: Set[int] = set()
    if args.resume and not args.figures_only and args.out_csv_full:
        existing_full_path = Path(args.out_csv_full)
        if existing_full_path.exists():
            try:
                existing_full_df = pd.read_csv(existing_full_path)
                if "program_id" in existing_full_df.columns:
                    existing_programs = set(
                        existing_full_df["program_id"].dropna().astype(int).tolist()
                    )
                logger.info(
                    "Resume enabled: found %d programs in existing full CSV.",
                    len(existing_programs),
                )
            except Exception as exc:
                logger.warning("Failed to read existing full CSV: %s", exc)

    cache_dir = Path(args.cache_dir) if args.cache_dir else None
    if cache_dir:
        cache_dir.mkdir(parents=True, exist_ok=True)

    program_to_results: Dict[str, List[Dict[str, Any]]] = {}
    total = len(program_to_genes)
    for idx, program_id in enumerate(sorted(program_to_genes.keys(), key=lambda x: int(x)), start=1):
        program_id_int = int(program_id)
        genes = [g for g in program_to_genes[program_id] if isinstance(g, str) and g.strip()]
        skip_enrichment = (
            args.resume and not args.force_refresh and program_id_int in existing_programs
        )

        results: Optional[List[Dict[str, Any]]] = None
        if not args.figures_only and not skip_enrichment:
            if cache_dir and not args.force_refresh:
                results = load_cached_results(cache_dir, program_id_int)
                if results is not None:
                    logger.info(
                        "Program %s: using cached enrichment results (%d terms).",
                        program_id,
                        len(results),
                    )
            if results is None:
                logger.info(
                    "[%d/%d] STRING enrichment for program %s with %d genes ...",
                    idx,
                    total,
                    program_id,
                    len(genes),
                )
                results = call_string_enrichment(
                    genes=genes,
                    species=args.species,
                    retries=args.retries,
                    sleep_between=args.sleep,
                )
                logger.info("Program %s: retrieved %d enriched terms", program_id, len(results))
                if cache_dir:
                    write_cached_results(cache_dir, program_id_int, results)
                time.sleep(args.sleep)

        if not args.figures_only and results is not None and not skip_enrichment:
            program_to_results[program_id] = results

    if not args.figures_only:
        df_full = build_full_csv(program_to_results)
        if existing_full_df is not None:
            df_full = pd.concat([existing_full_df, df_full], ignore_index=True)
        if not df_full.empty:
            df_full.sort_values(["program_id", "fdr", "p_value"], inplace=True)

        if args.out_csv_full:
            out_csv_full_path = Path(args.out_csv_full)
            out_csv_full_path.parent.mkdir(parents=True, exist_ok=True)
            df_full.to_csv(out_csv_full_path, index=False)
            logger.info(
                "Wrote full unfiltered CSV with %d rows to %s",
                len(df_full),
                out_csv_full_path,
            )

        if args.out_csv_filtered:
            out_csv_filtered_path = Path(args.out_csv_filtered)
            df_filtered = filter_process_kegg(df_full)
            out_csv_filtered_path.parent.mkdir(parents=True, exist_ok=True)
            df_filtered.to_csv(out_csv_filtered_path, index=False)
            logger.info(
                "Wrote filtered CSV (Process/KEGG) with %d rows to %s",
                len(df_filtered),
                out_csv_filtered_path,
            )

    # Optional: retrieve enrichment figures from STRING API
    if args.figures_dir:
        figures_dir = Path(args.figures_dir)
        figures_dir.mkdir(parents=True, exist_ok=True)

        def download_string_enrichment_figure(
            genes: List[str],
            species: int,
            category: str,
            output_path: Path,
            retries: int = 3,
        ) -> bool:
            """
            Download enrichment figure directly from STRING API.

            Args:
                genes: List of gene identifiers
                species: NCBI taxonomy ID
                category: Enrichment category (e.g., "Process", "KEGG")
                output_path: Path to save the figure
                retries: Number of retry attempts

            Returns:
                True if successful, False otherwise
            """
            if not genes:
                return False

            if args.resume and output_path.exists():
                return True

            # STRING enrichment figure endpoint
            base_url = "https://string-db.org/api/image/enrichmentfigure"

            # Prepare parameters
            identifiers_value = "\r".join(genes)
            params = {
                "identifiers": identifiers_value,
                "species": species,
                "category": category,
                "caller_identity": "topic_analysis_string_enrichment",
            }

            attempt = 0
            while attempt <= retries:
                try:
                    response = requests.get(base_url, params=params, timeout=120)
                    if response.status_code == 200:
                        # Check if response is actually an image
                        content_type = response.headers.get("content-type", "")
                        if "image" in content_type:
                            output_path.parent.mkdir(parents=True, exist_ok=True)
                            with open(output_path, "wb") as f:
                                f.write(response.content)
                            return True
                        else:
                            logger.warning(
                                "STRING returned non-image content for category %s: %s",
                                category,
                                content_type,
                            )
                            return False
                    else:
                        logger.warning(
                            "STRING figure API returned status %s", response.status_code
                        )
                except requests.RequestException as e:
                    logger.warning(
                        "HTTP error downloading figure (attempt %d/%d): %s",
                        attempt + 1,
                        retries + 1,
                        e,
                    )

                attempt += 1
                time.sleep(min(3.0, 1.0 * (attempt + 1)))

            return False

        for program_id, genes in program_to_genes.items():
            genes_list = [g for g in genes if isinstance(g, str) and g.strip()]

            # Download Process enrichment figure
            ok_p = download_string_enrichment_figure(
                genes_list,
                args.species,
                "Process",
                figures_dir / f"program_{program_id}_process_enrichment.png",
            )

            # Download KEGG enrichment figure
            ok_k = download_string_enrichment_figure(
                genes_list,
                args.species,
                "KEGG",
                figures_dir / f"program_{program_id}_kegg_enrichment.png",
            )

            logger.info(
                "Program %s: figures - Process=%s, KEGG=%s",
                program_id,
                "✓" if ok_p else "✗",
                "✓" if ok_k else "✗",
            )

            # Add delay between programs to avoid overwhelming the API
            time.sleep(0.5)

    return 0


# -------------------------------- Entry points -------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Extract program gene lists and run STRING enrichment.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # extract
    p_extract = subparsers.add_parser("extract", help="Extract top-N genes per program from loading CSV")
    p_extract.add_argument("--config", help="Path to config file (YAML or JSON)")
    p_extract.add_argument(
        "--input",
        help="Path to loading CSV (columns: Name, Score, RowID or program_id)",
    )
    p_extract.add_argument("--n-top", type=int, default=100, help="Number of top genes per program to extract")
    p_extract.add_argument("--json-out", help="Output JSON {program_id: [genes...]}")
    p_extract.add_argument("--csv-out", help="Output overview CSV")
    p_extract.add_argument(
        "--gene-loading-out",
        help="Output gene loading CSV with UniquenessScore (default: <input>_with_uniqueness.csv)",
    )
    p_extract.add_argument(
        "--topics",
        type=str,
        help="Comma-separated list of program IDs to include (e.g. '1,2,3')",
    )
    p_extract.set_defaults(func=cmd_extract)

    # enrich
    p_enrich = subparsers.add_parser("enrich", help="Run STRING enrichment for program gene lists from JSON")
    p_enrich.add_argument("--config", help="Path to config file (YAML or JSON)")
    p_enrich.add_argument("--genes-json", help="Path to JSON mapping {program_id: [genes...]}")
    p_enrich.add_argument("--species", type=int, default=10090, help="NCBI/STRING species id (default: 10090 mouse)")
    p_enrich.add_argument("--out-csv-full", help="Full unfiltered CSV output path")
    p_enrich.add_argument("--out-csv-filtered", help="Filtered CSV (Process/KEGG only, background<500)")
    p_enrich.add_argument("--figures-dir", help="Directory to save enrichment figures")
    p_enrich.add_argument(
        "--figures-only",
        action="store_true",
        help="Only download enrichment figures; skip enrichment CSV generation",
    )
    p_enrich.add_argument(
        "--cache-dir",
        help="Directory to cache per-program STRING enrichment JSON",
    )
    p_enrich.add_argument(
        "--resume",
        action="store_true",
        help="Skip programs already present in output CSVs and existing figures",
    )
    p_enrich.add_argument(
        "--force-refresh",
        action="store_true",
        help="Ignore cached results and re-query STRING",
    )
    p_enrich.add_argument("--sleep", type=float, default=0.6, help="Sleep seconds between API calls")
    p_enrich.add_argument("--retries", type=int, default=3, help="Retries per program on HTTP failures")
    p_enrich.add_argument(
        "--topics",
        type=str,
        help="Comma-separated list of program IDs to include (e.g. '1,2,3')",
    )
    p_enrich.set_defaults(func=cmd_enrich)

    # all
    p_all = subparsers.add_parser("all", help="Run extract then enrich")
    p_all.add_argument("--config", help="Path to config file (YAML or JSON)")
    # extract args
    p_all.add_argument("--input")
    p_all.add_argument("--n-top", type=int, default=100)
    p_all.add_argument("--json-out")
    p_all.add_argument("--csv-out")
    p_all.add_argument("--gene-loading-out")
    # enrich args
    p_all.add_argument("--species", type=int, default=10090)
    p_all.add_argument("--out-csv-full")
    p_all.add_argument("--out-csv-filtered")
    p_all.add_argument("--figures-dir")
    p_all.add_argument("--figures-only", action="store_true")
    p_all.add_argument("--cache-dir")
    p_all.add_argument("--resume", action="store_true")
    p_all.add_argument("--force-refresh", action="store_true")
    p_all.add_argument("--sleep", type=float, default=0.6)
    p_all.add_argument("--retries", type=int, default=3)
    p_all.add_argument(
        "--topics",
        type=str,
        help="Comma-separated list of program IDs to include (e.g. '1,2,3')",
    )

    def run_all(args: argparse.Namespace) -> int:
        rc = cmd_extract(args)
        if rc != 0:
            return rc
        enrich_args = argparse.Namespace(
            genes_json=args.json_out,
            species=args.species,
            out_csv_full=args.out_csv_full,
            out_csv_filtered=args.out_csv_filtered,
            figures_dir=args.figures_dir,
            figures_only=args.figures_only,
            cache_dir=args.cache_dir,
            resume=args.resume,
            force_refresh=args.force_refresh,
            sleep=args.sleep,
            retries=args.retries,
            topics=args.topics,
            func=cmd_enrich,
        )
        return cmd_enrich(enrich_args)

    p_all.set_defaults(func=run_all)

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    config = load_config(getattr(args, "config", None))
    cli_overrides = get_cli_overrides(sys.argv)
    args = apply_config_overrides(args, config, cli_overrides)
    args = apply_test_mode(args, config, cli_overrides)
    args = apply_default_paths(args)
    rc = args.func(args)
    raise SystemExit(rc)


if __name__ == "__main__":
    main()
