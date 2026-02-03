"""
Test script for flexible column name matching.
"""
from pathlib import Path
import sys

# Add pipeline directory to path
sys.path.insert(0, str(Path(__file__).parent.parent / "pipeline"))

import pandas as pd
from column_mapper import ColumnMapper

# Test gene loading formats
print("=" * 60)
print("Testing Gene Loading Column Matching")
print("=" * 60)

test_cases = [
    # Standard format
    {'Name': ['Gene1', 'Gene2'], 'Score': [0.8, 0.6], 'program_id': [1, 1]},
    # Alternative format 1
    {'Gene': ['Gene1', 'Gene2'], 'Loading': [0.8, 0.6], 'RowID': [1, 1]},
    # Alternative format 2  
    {'gene_name': ['Gene1', 'Gene2'], 'score': [0.8, 0.6], 'Topic': [1, 1]},
    # Mixed case
    {'GENE': ['Gene1', 'Gene2'], 'SCORE': [0.8, 0.6], 'topic_id': [1, 1]},
]

for i, case in enumerate(test_cases, 1):
    df = pd.DataFrame(case)
    print(f"\nTest Case {i}: {list(df.columns)}")
    try:
        mapper = ColumnMapper(df)
        gene_col = mapper.get_column('gene')
        score_col = mapper.get_column('score')
        prog_col = mapper.get_column('program_id')
        print(f"  ✓ Mapped: gene={gene_col}, score={score_col}, program_id={prog_col}")
    except Exception as e:
        print(f"  ✗ Error: {e}")

# Test cell-type enrichment formats
print("\n" + "=" * 60)
print("Testing Cell-Type Enrichment Column Matching")
print("=" * 60)

ct_test_cases = [
    # Standard format
    {'cell_type': ['Artery'], 'program': ['Program_5'], 'log2_fc_in_vs_out': [1.5], 'fdr': [0.01]},
    # Alternative format 1
    {'cluster': ['Artery'], 'topic': ['5'], 'log2FC': [1.5], 'p_adj': [0.01]},
    # Alternative format 2
    {'CellType': ['Artery'], 'Program': ['Program_5'], 'Log2FoldChange': [1.5], 'FDR': [0.01]},
    # Alternative format 3
    {'Cluster': ['Artery'], 'RowID': ['5'], 'lfc': [1.5], 'qvalue': [0.01]},
]

for i, case in enumerate(ct_test_cases, 1):
    df = pd.DataFrame(case)
    print(f"\nTest Case {i}: {list(df.columns)}")
    try:
        mapper = ColumnMapper(df)
        ct_col = mapper.get_column('cell_type')
        prog_col = mapper.get_column('program_id')
        fc_col = mapper.get_column('log2_fc')
        fdr_col = mapper.get_column('fdr')
        print(f"  ✓ Mapped: cell_type={ct_col}, program_id={prog_col}, log2_fc={fc_col}, fdr={fdr_col}")
    except Exception as e:
        print(f"  ✗ Error: {e}")

print("\n" + "=" * 60)
print("All Tests Complete!")
print("=" * 60)
