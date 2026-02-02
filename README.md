# Topic Annotation Pipeline

Automated pipeline for gene program annotation using LLM-based analysis with STRING enrichment and literature evidence.

---

## Quick Start

The pipeline requires **3 input files** and runs all steps automatically:

```bash
# Full pipeline run
python pipeline/run_pipeline.py \
  --gene-loading input/genes/FB_moi15_seq2_loading_gene_k100_top300.csv \
  --celltype-enrichment input/celltype/fp_seq2_cnmf_celltype_l2_program_enrichment.csv \
  --regulator-file input/regulators/sceptre_discovery_analysis_results.csv \
  --output-dir results/output/my_run

# Or use the config file
python pipeline/run_pipeline.py --config configs/pipeline_config.yaml

# Test with specific topics (faster)
python pipeline/run_pipeline.py --config configs/pipeline_config.yaml --topics 5,6,8,11,18
```

---

## Pipeline Steps

| Step | Name | Description |
|------|------|-------------|
| 1 | `string_enrichment` | Extract top genes, run STRING enrichment |
| 2 | `literature_fetch` | Fetch gene summaries & PubMed literature (25 papers/program) |
| 3a | `batch_prepare` | Prepare LLM batch request |
| 3b | `batch_submit` | Submit to Anthropic Batch API (default) or Vertex AI |
| 4 | `parse_results` | Parse LLM responses, generate summary |
| 5 | `html_report` | Generate interactive HTML report |

---

## Required Inputs

3 files are required:

| File | Description |
|------|-------------|
| `input/genes/FB_moi15_seq2_loading_gene_k100_top300.csv` | Gene loading matrix (columns: Name, Score, RowID or program_id) |
| `input/celltype/fp_seq2_cnmf_celltype_l2_program_enrichment.csv` | Raw cell-type enrichment (columns: cell_type, program, log2_fc_in_vs_out, fdr) |
| `input/regulators/sceptre_discovery_analysis_results.csv` | SCEPTRE regulator results (Perturb-seq) |

### Column compatibility (key files)

| File | Required columns |
|------|------------------|
| Gene loading (`gene_loading`) | `Name`, `Score`, `RowID` **or** `program_id` |
| Cell-type enrichment (`celltype_enrichment`) | `cell_type`, `program` (format `Program_<id>`), `log2_fc_in_vs_out`, `fdr` |
| STRING enrichment outputs | `program_id`, `category`, `term`, `fdr`, `p_value`, `inputGenes` |
| Cell-type summary (auto-generated) | `program`, `highly_cell_type_specific`, `moderately_enriched`, `weakly_enriched`, `significantly_lower_expression` |

---

## Configuration

Edit `configs/pipeline_config.yaml`:

```yaml
input:
  gene_loading: input/genes/FB_moi15_seq2_loading_gene_k100_top300.csv
  celltype_enrichment: input/celltype/fp_seq2_cnmf_celltype_l2_program_enrichment.csv
  regulator_file: input/regulators/sceptre_discovery_analysis_results.csv

output_dir: results/output/pipeline_run

# Limit to specific topics (null = all)
topics: null

# Species for STRING-DB (10090 = mouse, 9606 = human)
species: 10090

# Literature search context
context: '(endothelial OR endothelium OR "vascular endothelial")'

# LLM settings
llm_backend: anthropic           # "anthropic" (default) or "vertex"
llm_model: claude-4-sonnet-20250514
llm_max_tokens: 8192
llm_wait: true                   # Wait for batch completion

# Gene summary options
# full_summaries: true           # Use full HTML summaries (~3x longer)
```

---

## LLM Backend Options

### Anthropic Batch API (Default, Recommended)

Uses Anthropic's native Batch API directly. Faster and simpler for most use cases.

```bash
# Automatic via pipeline
python pipeline/run_pipeline.py --config configs/pipeline_config.yaml

# Manual submission
python pipeline/03_submit_and_monitor_batch.py submit batch_request.json --wait

# Check status
python pipeline/03_submit_and_monitor_batch.py check --batch-id <BATCH_ID>

# Download results
python pipeline/03_submit_and_monitor_batch.py results --batch-id <BATCH_ID>
```

**Required**: `ANTHROPIC_API_KEY` environment variable.

### Vertex AI (Alternative)

Uses Google Cloud Vertex AI for batch processing.

```yaml
# In config:
llm_backend: vertex
```

```bash
# Manual submission
python pipeline/03_submit_and_monitor_batch.py submit-vertex batch_request.json --wait

# Check status
python pipeline/03_submit_and_monitor_batch.py check-vertex --job-name <JOB_NAME>
```

**Required**: Google Cloud SDK configured with appropriate permissions.

---

## Gene Summary Options

Gene summaries are fetched from **Harmonizome** (curated NCBI-like descriptions).

| Option | Description | Length |
|--------|-------------|--------|
| Default | Short API descriptions | ~400-600 chars |
| `--full-summaries` | Full HTML summaries with PMID references | ~2000+ chars |
| `--gene-summary-source ncbi` | Direct NCBI Entrez summaries | Variable |

```bash
# Use full summaries (richer context, slower)
python pipeline/02_fetch_ncbi_data.py ... --full-summaries

# Or in config:
full_summaries: true
```

---

## Partial Runs

Run specific steps:

```bash
# Stop after preparing batch (steps 1-3a, no submission)
python pipeline/run_pipeline.py --config configs/pipeline_config.yaml --stop-after batch_prepare

# Resume from parsing - Anthropic (local results file)
python pipeline/run_pipeline.py --config configs/pipeline_config.yaml \
  --start-from parse_results \
  --results-file results/output/my_run/llm_batches/batch_request_results.jsonl

# Resume from parsing - Vertex AI (GCS prefix)
python pipeline/run_pipeline.py --config configs/pipeline_config.yaml \
  --start-from parse_results \
  --gcs-prefix gs://perturbseq/batch/prediction-model-2024...
```

Available steps for `--start-from` and `--stop-after`:
- `string_enrichment`
- `literature_fetch`
- `batch_prepare`
- `batch_submit`
- `parse_results`
- `html_report`

---

## Output Structure

All outputs are organized under the specified `output_dir`:

```
results/output/my_run/
├── genes_top.json                     # Program → gene list mapping
├── genes_overview.csv                 # Gene loading overview
├── gene_loading_with_uniqueness.csv   # Gene table with UniquenessScore
├── celltype_summary.csv               # Auto-generated cell-type summary
├── literature_context.json            # Gene summaries & literature evidence
├── string_enrichment/
│   ├── enrichment_full.csv           # All STRING enrichment terms
│   ├── enrichment_filtered.csv       # Filtered (Process/KEGG only)
│   ├── figures/                      # Enrichment figures (PNG)
│   └── cache/                        # API response cache
├── llm_batches/
│   └── batch_request.json            # LLM batch request payload
└── annotations/
    ├── topic_*_annotation.md         # Per-topic markdown files
    ├── summary.csv                   # Topic names and summaries
    └── report.html                   # Interactive HTML report
```

---

## Prerequisites

- Python environment with: `pandas`, `requests`, `numpy`, `markdown`, `pyyaml`, `anthropic`, `matplotlib`, `seaborn`, `tqdm`, `pillow`
- Google Cloud SDK (only for Vertex AI batch submission)
- API keys in `.env`:
  - `ANTHROPIC_API_KEY` (required for Anthropic backend)
  - `NCBI_API_KEY` (optional, increases PubMed throughput)

Recommended env spec: `configs/environment.yaml` (conda) covers core pipeline deps and optional Vertex AI support.

---

## Individual Step Scripts

For advanced usage, each step can be run independently:

| Script | Purpose |
|--------|---------|
| `pipeline/01_genes_to_string_enrichment.py` | STRING enrichment |
| `pipeline/02_fetch_ncbi_data.py` | Literature fetching |
| `pipeline/03_submit_and_monitor_batch.py` | Batch preparation/submission |
| `pipeline/04_parse_and_summarize.py` | Result parsing |
| `pipeline/05_generate_html_report.py` | HTML report generation |

Run `python pipeline/<script>.py --help` for detailed options.

### Batch Submission Commands

```bash
# Prepare batch request only
python pipeline/03_submit_and_monitor_batch.py prepare \
  --gene-file input/genes/FB_moi15_seq2_loading_gene_k100_top300.csv \
  --output-file results/output/batch_request.json

# Submit to Anthropic (default)
python pipeline/03_submit_and_monitor_batch.py submit batch_request.json --wait

# Submit to Vertex AI (alternative)
python pipeline/03_submit_and_monitor_batch.py submit-vertex batch_request.json --wait
```

---

## HTML Report Features

- **Program stats**: Top-loading genes, unique genes, cell-type enrichment
- **LLM annotations**: Full markdown with functional modules
- **Enrichment figures**: KEGG + Process enrichment panels
- **Interactive volcano plot**: Priority gene labels
- **Full-text search**: Search across all annotation content
- **Keyboard navigation**: Arrow keys (←/→) to switch programs
