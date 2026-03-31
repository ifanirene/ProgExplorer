"""
Unit tests for flexible regulator input parsing and ranking.
"""

from pathlib import Path
import importlib.util
import sys

import pandas as pd


repo_root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(repo_root / "pipeline"))

from column_mapper import standardize_regulator_results


def _load_module(module_name: str, relative_path: str):
    spec = importlib.util.spec_from_file_location(
        module_name, repo_root / relative_path
    )
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


pipeline_step02 = _load_module("pipeline_step02", "pipeline/02_fetch_ncbi_data.py")
report_module = _load_module("report_module", "pipeline/05_generate_html_report.py")


def test_standardize_regulator_results_supports_aliases_and_adj_threshold():
    df = pd.DataFrame(
        {
            "program_name": ["Program_5", "Program_5"],
            "target_gene_names": ["Kdr", "Flt1"],
            "log2fc": [-1.2, 0.4],
            "adj_pval": [0.001, 0.2],
        }
    )

    standardized = standardize_regulator_results(df, significance_threshold=0.05)

    assert standardized["program_id"].tolist() == [5, 5]
    assert standardized["grna_target"].tolist() == ["Kdr", "Flt1"]
    assert standardized["p_value"].tolist() == [0.001, 0.2]
    assert standardized["significant"].tolist() == [True, False]


def test_standardize_regulator_results_prefers_explicit_significant_column():
    df = pd.DataFrame(
        {
            "response_id": ["X7", "X7"],
            "grna_target": ["Foxo1", "Eng"],
            "log_2_fold_change": [-0.7, 0.8],
            "adj_pval": [0.2, 0.001],
            "significant": ["FALSE", "TRUE"],
        }
    )

    standardized = standardize_regulator_results(df, significance_threshold=0.01)

    assert standardized["program_id"].tolist() == [7, 7]
    assert standardized["significant"].tolist() == [False, True]


def test_get_top_regulators_uses_separate_positive_negative_limits():
    regulator_data = {
        5: pd.DataFrame(
            {
                "grna_target": ["A", "B", "C", "D", "E"],
                "log_2_fold_change": [-2.0, -1.0, 0.5, 1.5, 2.5],
                "p_value": [0.001] * 5,
                "significant": [True] * 5,
            }
        )
    }

    result = pipeline_step02.get_top_regulators(
        regulator_data,
        5,
        top_n_positive=2,
        top_n_negative=1,
    )

    assert [row["gene"] for row in result["positive"]] == ["A", "B"]
    assert [row["gene"] for row in result["negative"]] == ["E"]


def test_generate_report_prefers_adjusted_pvalues_for_volcano(tmp_path):
    summary_csv = tmp_path / "summary.csv"
    summary_csv.write_text("Topic,Name\n5,Program 5\n", encoding="utf-8")

    annotations_dir = tmp_path / "annotations"
    annotations_dir.mkdir()
    (annotations_dir / "topic_5_annotation.md").write_text(
        "**Program label:** Test\n**Brief Summary:** Example\n",
        encoding="utf-8",
    )

    enrichment_dir = tmp_path / "enrichment"
    enrichment_dir.mkdir()

    volcano_csv = tmp_path / "regulators.csv"
    pd.DataFrame(
        {
            "program_name": ["Program_5"],
            "target_gene_names": ["Kdr"],
            "log2fc": [-1.0],
            "p_value": [1e-10],
            "adj_pval": [1e-4],
        }
    ).to_csv(volcano_csv, index=False)

    gene_loading_csv = tmp_path / "genes.csv"
    gene_loading_csv.write_text("Name,Score,program_id\nKdr,1.0,5\n", encoding="utf-8")

    output_html = tmp_path / "report.html"
    report_module.generate_report(
        str(summary_csv),
        str(annotations_dir),
        str(enrichment_dir),
        str(volcano_csv),
        str(gene_loading_csv),
        str(output_html),
        regulator_significance_threshold=0.05,
    )

    html = output_html.read_text(encoding="utf-8")
    assert '"g": "Kdr"' in html
    assert '"p": 4.0' in html
