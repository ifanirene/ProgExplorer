import importlib.util
from pathlib import Path

import pandas as pd


def load_module():
    module_path = Path(__file__).resolve().parents[1] / "pipeline" / "01_genes_to_string_enrichment.py"
    spec = importlib.util.spec_from_file_location("genes_to_string_enrichment", module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_build_uniqueness_table_from_rowid():
    module = load_module()
    df = pd.DataFrame(
        {
            "Name": ["GeneA", "GeneB", "GeneA", "GeneC"],
            "Score": [1.0, 0.5, 0.8, 0.2],
            "RowID": [1, 1, 2, 2],
        }
    )
    out = module.build_uniqueness_table(df, id_col="RowID")
    assert "program_id" in out.columns
    assert "UniquenessScore" in out.columns
    assert out["UniquenessScore"].notna().all()
    assert out["program_id"].nunique() == 2


def test_default_uniqueness_output_suffix():
    module = load_module()
    path = Path("input/genes/example.csv")
    out = module.default_uniqueness_output(path)
    assert out.name == "example_with_uniqueness.csv"
