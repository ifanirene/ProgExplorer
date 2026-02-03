# Flexible Input Format Support

The ProgExplorer pipeline supports diverse input file formats through intelligent column name matching. This makes it compatible with outputs from different analysis tools (cNMF, Seurat, Scanpy, etc.) without requiring file reformatting.

## How It Works

The `column_mapper.py` module provides:
- **Case-insensitive matching**: `Gene`, `gene`, `GENE` all work
- **Alias support**: Multiple common names for each column type
- **Clear error messages**: Shows expected column names when validation fails

## Supported Column Names

### Gene Loading Files

| Standard Name | Supported Aliases |
|---------------|-------------------|
| **Gene name** | Name, Gene, gene_name, GeneName, gene_symbol, GeneSymbol, Symbol, gene_id, GeneID |
| **Score/Loading** | Score, score, Loading, loading, Weight, weight, Value, value, Loading_Score, gene_score |
| **Program ID** | program_id, programid, RowID, row_id, Topic, topic, topic_id, Program, program_number, topic_number, K, k, Component, component, Factor, factor |

### Cell-Type Enrichment Files

| Standard Name | Supported Aliases |
|---------------|-------------------|
| **Cell type** | cell_type, celltype, Cluster, cluster, cell_cluster, CellCluster, annotation, cell_annotation, Type, type, celltype_id, cluster_id |
| **Program ID** | program, program_id, programid, RowID, row_id, Topic, topic, topic_id, K, k |
| **Log2 fold change** | log2_fc, log2FC, log2_fold_change, log2FoldChange, log2_fc_in_vs_out, log2fc_in_vs_out, log2_fold, log2fold, LFC, lfc, L2FC, l2fc, fold_change, FoldChange, fc, FC, log2ratio, log2_ratio |
| **FDR/p-value** | fdr, FDR, fdr_corrected, q_value, qvalue, qval, Q, q, padj, p_adj, p_adjusted, adjusted_pvalue, p_value, pvalue, pval, P, p |

## Examples

### Gene Loading Formats

All of these work:

```csv
# Format 1: Standard cNMF output
Name,Score,program_id
Gja1,0.85,5
Notch1,0.72,5

# Format 2: Alternative names
Gene,Loading,RowID
Gja1,0.85,5
Notch1,0.72,5

# Format 3: Seurat-style
gene_name,score,Topic
Gja1,0.85,5
Notch1,0.72,5

# Format 4: All caps
GENE,SCORE,TOPIC_ID
Gja1,0.85,5
Notch1,0.72,5
```

### Cell-Type Enrichment Formats

All of these work:

```csv
# Format 1: Standard
cell_type,program,log2_fc_in_vs_out,fdr
Artery,Program_5,1.5,0.001

# Format 2: Seurat cluster style
cluster,topic,log2FC,p_adj
Artery,5,1.5,0.001

# Format 3: Mixed case
CellType,Program,Log2FoldChange,FDR
Artery,Program_5,1.5,0.001

# Format 4: Alternative names
Cluster,RowID,lfc,qvalue
Artery,5,1.5,0.001
```

## Program ID Formats

The pipeline automatically parses these program identifier formats:
- `Program_5`, `program_5`
- `Topic_5`, `topic_5`
- `P_5`, `p_5`
- `5_5` (regulator file format)
- `Program5`, `ProgramNumber5`
- `Topic5`, `TopicNumber5`
- `55` (plain integer)
- `5` (plain integer)

## Testing Compatibility

Test if your files are compatible:

```bash
python tests/test_column_mapper.py
```

This will show you all supported column name variants.

## Implementation Details

The `ColumnMapper` class:
1. Creates a lowercase mapping of all column names
2. Tries to match against known aliases
3. Returns the actual column name from the DataFrame
4. Provides helpful error messages with expected aliases when columns are missing

## Adding New Aliases

To add support for new column name variants, edit `pipeline/column_mapper.py` and add aliases to the `ALIASES` dictionary:

```python
ALIASES = {
    'gene': [
        'gene', 'name', 'gene_name', ...,
        'your_new_alias',  # Add here
    ],
    ...
}
```
