import importlib.util
from pathlib import Path

import pandas as pd


def load_module():
    module_path = Path(__file__).resolve().parents[1] / "pipeline" / "03_submit_and_monitor_batch.py"
    spec = importlib.util.spec_from_file_location("submit_and_monitor_batch", module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_load_gene_table_computes_uniqueness(tmp_path):
    df = pd.DataFrame(
        {
            "Name": ["GeneA", "GeneB", "GeneA", "GeneC"],
            "Score": [1.0, 0.5, 0.8, 0.2],
            "RowID": [1, 1, 2, 2],
        }
    )
    csv_path = tmp_path / "genes.csv"
    df.to_csv(csv_path, index=False)

    module = load_module()
    loaded = module.load_gene_table(csv_path)

    assert "program_id" in loaded.columns
    assert "UniquenessScore" in loaded.columns
    assert loaded["UniquenessScore"].notna().any()
