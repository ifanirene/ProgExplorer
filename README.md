# Topic Annotation Pipeline

Automated pipeline for annotating gene programs using LLM analysis, STRING enrichment, and literature evidence.

## Quick Start

Run the full pipeline with a configuration file:

```bash
python pipeline/run_pipeline.py --config configs/pipeline_config.yaml
```

Or specify inputs directly:

```bash
python pipeline/run_pipeline.py \
  --gene-loading input/genes/FB_moi15_seq2_loading_gene_k100_top300.csv \
  --celltype-enrichment input/celltype/fp_seq2_cnmf_celltype_l2_program_enrichment.csv \
  --regulator-file input/regulators/sceptre_discovery_analysis_results.csv \
  --output-dir results/output/my_run
```

Test with specific topics (faster):

```bash
python pipeline/run_pipeline.py --config configs/pipeline_config.yaml --topics 5,6,8,11,18
```

## Prerequisites

### Environment Setup

Install dependencies using conda:

```bash
conda env create -f configs/environment.yaml
conda activate progexplorer
```

Or manually install: `pandas`, `requests`, `numpy`, `markdown`, `pyyaml`, `anthropic`, `matplotlib`, `seaborn`, `tqdm`, `pillow`

### API Keys

Set environment variables in `.env`:

- `ANTHROPIC_API_KEY` (required for LLM annotations)
- `NCBI_API_KEY` (optional, increases PubMed rate limits)

### Optional: Vertex AI

For Vertex AI backend instead of Anthropic:
- Install Google Cloud SDK
- Configure appropriate GCP permissions

## Required Inputs

The pipeline needs 3 input files:

| File | Description | Required Columns |
|------|-------------|------------------|
| **Gene loading** | Gene loading matrix from cNMF or similar | `Name`, `Score`, `RowID` or `program_id` |
| **Cell-type enrichment** | Raw cell-type enrichment results | `cell_type`, `program`, `log2_fc_in_vs_out`, `fdr` |
| **Regulator file** | SCEPTRE or Perturb-seq regulator results | Standard SCEPTRE output format |

**Note**: `program` column must use format `Program_<id>` (e.g., `Program_5`).

## Configuration

Edit `configs/pipeline_config.yaml` to set inputs and parameters:

```yaml
input:
  gene_loading: input/genes/FB_moi15_seq2_loading_gene_k100_top300.csv
  celltype_enrichment: input/celltype/fp_seq2_cnmf_celltype_l2_program_enrichment.csv
  regulator_file: input/regulators/sceptre_discovery_analysis_results.csv

output_dir: results/output/pipeline_run

topics: null  # null = all topics, or list specific: [5, 6, 8, 11, 18]
species: 10090  # 10090 = mouse, 9606 = human
context: '(endothelial OR endothelium OR "vascular endothelial")'

llm_backend: anthropic  # "anthropic" or "vertex"
llm_model: claude-4-sonnet-20250514
llm_max_tokens: 8192
llm_wait: true  # Wait for batch to complete
```

## Pipeline Steps

The pipeline runs 6 steps automatically:

1. **String enrichment** - Extract top genes, compute UniquenessScore, run STRING enrichment, generate cell-type summary
2. **Literature fetch** - Fetch Harmonizome gene summaries and PubMed literature (25 papers/program)
3. **Batch prepare** - Generate LLM prompts with all evidence
4. **Batch submit** - Submit to Anthropic or Vertex AI batch API
5. **Parse results** - Extract LLM annotations and generate summary
6. **HTML report** - Create interactive HTML report with figures

## Partial Pipeline Runs

Stop at a specific step:

```bash
# Prepare batch without submitting
python pipeline/run_pipeline.py --config configs/pipeline_config.yaml --stop-after batch_prepare
```

Resume from a specific step:

```bash
# Resume from parsing (Anthropic)
python pipeline/run_pipeline.py --config configs/pipeline_config.yaml --start-from parse_results

# Resume from parsing (Vertex AI - requires GCS prefix)
python pipeline/run_pipeline.py --config configs/pipeline_config.yaml \
  --start-from parse_results \
  --gcs-prefix gs://bucket/batch/prediction-model-<timestamp>/
```

Available step names: `string_enrichment`, `literature_fetch`, `batch_prepare`, `batch_submit`, `parse_results`, `html_report`

## Output Structure

```
results/output/my_run/
├── genes_top.json                      # Program → gene list mapping
├── gene_loading_with_uniqueness.csv    # Gene table with UniquenessScore
├── celltype_summary.csv                # Auto-generated cell-type summary
├── literature_context.json             # Gene summaries & literature
├── string_enrichment/
│   ├── enrichment_filtered.csv         # Process/KEGG enrichment terms
│   └── figures/                        # Enrichment bar charts (PNG)
├── llm_batches/
│   ├── batch_request.json              # LLM prompts
│   └── batch_request_results.jsonl     # LLM responses
└── annotations/
    ├── topic_*_annotation.md           # Per-topic annotations
    ├── summary.csv                     # Topic names and summaries
    └── report.html                     # Interactive HTML report
```

## Advanced Options

### Use Full Gene Summaries

By default, short Harmonizome descriptions are used. For longer summaries with PMID references:

```yaml
full_summaries: true  # In config
```

Or:

```bash
python pipeline/02_fetch_ncbi_data.py ... --full-summaries
```

### Switch LLM Backend

**Anthropic (default):**

```yaml
llm_backend: anthropic
```

```bash
# Manual batch operations
python pipeline/03_submit_and_monitor_batch.py submit batch_request.json --wait
python pipeline/03_submit_and_monitor_batch.py check --batch-id <BATCH_ID>
python pipeline/03_submit_and_monitor_batch.py results --batch-id <BATCH_ID>
```

**Vertex AI:**

```yaml
llm_backend: vertex
```

```bash
python pipeline/03_submit_and_monitor_batch.py submit-vertex batch_request.json --wait
python pipeline/03_submit_and_monitor_batch.py check-vertex --job-name <JOB_NAME>
```

### Run Individual Scripts

For debugging or custom workflows, run pipeline steps independently:

```bash
# Step 1: STRING enrichment + cell-type summary
python pipeline/01_genes_to_string_enrichment.py all \
  --input input/genes/FB_moi15_seq2_loading_gene_k100_top300.csv \
  --celltype-enrichment input/celltype/fp_seq2_cnmf_celltype_l2_program_enrichment.csv \
  --topics 5,6,8

# Step 2: Literature fetch
python pipeline/02_fetch_ncbi_data.py \
  --input results/output/gene_loading_with_uniqueness.csv \
  --json-out results/output/literature_context.json \
  --topics 5,6,8

# Step 3: Prepare batch
python pipeline/03_submit_and_monitor_batch.py prepare \
  --gene-file results/output/gene_loading_with_uniqueness.csv \
  --celltype-file results/output/celltype_summary.csv \
  --output-file results/output/batch_request.json
```

Use `--help` on any script for full options.

## Troubleshooting

**Column name errors**: Ensure your gene loading file has `Name`, `Score`, and either `RowID` or `program_id`.

**Cell-type format errors**: The `program` column must use format `Program_<id>` (e.g., `Program_5`, not just `5`).

**API rate limits**: Set `NCBI_API_KEY` for higher PubMed throughput. Anthropic batch API has no rate limits.

**Missing dependencies**: Run `conda env create -f configs/environment.yaml` or manually install packages listed in Prerequisites.
