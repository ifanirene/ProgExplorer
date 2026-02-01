import argparse
import importlib.util
import json
from pathlib import Path

import pandas as pd


def load_module():
    module_path = Path(__file__).resolve().parents[1] / "pipeline" / "01_genes_to_string_enrichment.py"
    spec = importlib.util.spec_from_file_location("genes_to_string_enrichment", module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def write_genes_json(path: Path, mapping: dict[str, list[str]]) -> None:
    path.write_text(json.dumps(mapping, indent=2), encoding="utf-8")


def test_figures_only_requires_figures_dir(tmp_path):
    module = load_module()
    genes_json = tmp_path / "genes.json"
    write_genes_json(genes_json, {"1": ["GeneA", "GeneB"]})

    args = argparse.Namespace(
        genes_json=str(genes_json),
        figures_only=True,
        figures_dir=None,
        out_csv_full=None,
        out_csv_filtered=None,
        species=10090,
        sleep=0.0,
        retries=0,
        topics=None,
        cache_dir=None,
        resume=False,
        force_refresh=False,
    )

    rc = module.cmd_enrich(args)
    assert rc == 2


def test_cache_used_without_network(tmp_path, monkeypatch):
    module = load_module()
    genes_json = tmp_path / "genes.json"
    write_genes_json(genes_json, {"1": ["GeneA", "GeneB"]})

    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    cached_payload = [
        {
            "category": "Process",
            "term": "test",
            "term_id": "T1",
            "description": "test",
            "fdr": 0.01,
            "p_value": 0.02,
            "number_of_genes": 5,
            "number_of_genes_in_background": 100,
            "ncbiTaxonId": 10090,
            "inputGenes": ["GeneA", "GeneB"],
        }
    ]
    (cache_dir / "program_1_enrichment.json").write_text(
        json.dumps(cached_payload), encoding="utf-8"
    )

    def fail_call(*_args, **_kwargs):
        raise AssertionError("call_string_enrichment should not be invoked when cache exists")

    monkeypatch.setattr(module, "call_string_enrichment", fail_call)

    out_full = tmp_path / "full.csv"
    out_filtered = tmp_path / "filtered.csv"
    args = argparse.Namespace(
        genes_json=str(genes_json),
        figures_only=False,
        figures_dir=None,
        out_csv_full=str(out_full),
        out_csv_filtered=str(out_filtered),
        species=10090,
        sleep=0.0,
        retries=0,
        topics=None,
        cache_dir=str(cache_dir),
        resume=False,
        force_refresh=False,
    )

    rc = module.cmd_enrich(args)
    assert rc == 0
    df_full = pd.read_csv(out_full)
    assert set(df_full["program_id"]) == {1}


def test_resume_merges_existing_full(tmp_path, monkeypatch):
    module = load_module()
    genes_json = tmp_path / "genes.json"
    write_genes_json(genes_json, {"1": ["GeneA"], "2": ["GeneB"]})

    existing_full = tmp_path / "full.csv"
    pd.DataFrame(
        [
            {
                "program_id": 1,
                "category": "Process",
                "term": "existing",
                "term_id": "E1",
                "description": "existing",
                "fdr": 0.05,
                "p_value": 0.1,
                "number_of_genes": 3,
                "number_of_genes_in_background": 100,
                "ncbiTaxonId": 10090,
                "inputGenes": "GeneA",
            }
        ]
    ).to_csv(existing_full, index=False)

    called = []

    def fake_call(genes, **_kwargs):
        called.append(tuple(genes))
        return [
            {
                "category": "Process",
                "term": "new",
                "term_id": "N1",
                "description": "new",
                "fdr": 0.01,
                "p_value": 0.02,
                "number_of_genes": 4,
                "number_of_genes_in_background": 100,
                "ncbiTaxonId": 10090,
                "inputGenes": genes,
            }
        ]

    monkeypatch.setattr(module, "call_string_enrichment", fake_call)

    out_filtered = tmp_path / "filtered.csv"
    args = argparse.Namespace(
        genes_json=str(genes_json),
        figures_only=False,
        figures_dir=None,
        out_csv_full=str(existing_full),
        out_csv_filtered=str(out_filtered),
        species=10090,
        sleep=0.0,
        retries=0,
        topics=None,
        cache_dir=None,
        resume=True,
        force_refresh=False,
    )

    rc = module.cmd_enrich(args)
    assert rc == 0
    df_full = pd.read_csv(existing_full)
    assert set(df_full["program_id"]) == {1, 2}
    assert called == [("GeneB",)]
