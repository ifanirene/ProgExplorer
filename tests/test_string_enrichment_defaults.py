import argparse
import importlib.util
from pathlib import Path


def load_module():
    module_path = Path(__file__).resolve().parents[1] / "pipeline" / "01_genes_to_string_enrichment.py"
    spec = importlib.util.spec_from_file_location("genes_to_string_enrichment", module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_defaults_for_all_command():
    module = load_module()
    args = argparse.Namespace(
        command="all",
        input="input/genes/FB_moi15_seq2_loading_gene_k100_top300.csv",
        n_top=300,
        json_out=None,
        csv_out=None,
        species=10090,
        out_csv_full=None,
        out_csv_filtered=None,
        figures_dir=None,
        sleep=0.6,
        retries=3,
        topics=None,
    )

    args = module.apply_default_paths(args)

    assert args.json_out == "input/enrichment/genes_top300.json"
    assert args.csv_out == "input/enrichment/genes_overview_top300.csv"
    assert args.out_csv_full == "input/enrichment/string_enrichment_full.csv"
    assert args.out_csv_filtered == "input/enrichment/string_enrichment_filtered_process_kegg.csv"
    assert args.figures_dir is None
