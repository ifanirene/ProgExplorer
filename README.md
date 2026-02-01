# Topic Annotation Workflow

Automated pipeline for gene program annotation using LLM-based analysis with STRING enrichment and NCBI literature evidence.

---

## Scripts Overview

| Step | Script | Purpose |
|------|--------|---------|
| 1 | `01_genes_to_string_enrichment.py` | Extract top genes, run STRING enrichment |
| 2 | `02_fetch_ncbi_data.py` | Fetch gene summaries (NCBI or Harmonizome) & literature |
| 3 | `03_submit_and_monitor_batch.py` | Prepare/submit/monitor batch jobs (Anthropic + Vertex AI) |
| 4 | `04_parse_and_summarize.py` | Parse results, generate summary CSV |
| 5 | `05_generate_html_report.py` | Generate interactive HTML report |

Utility: `extract_prompts_to_md.py` - Extract prompts to markdown files

---

## Prerequisites

- Conda env with `pandas`, `requests`, `python-dotenv`, `markdown`, `numpy`
- `ANTHROPIC_API_KEY` in `.env`
- `NCBI_API_KEY` in `.env` (optional, increases throughput)
- Google Cloud SDK for Vertex AI batch jobs

## Example Inputs & Outputs

This repo ships full inputs under `input/` with a simple layout:

- `input/genes/FB_moi15_seq2_loading_gene_k100_top300.csv`
- `input/regulators/sceptre_discovery_analysis_results.csv`
- `input/celltype/program_celltype_annotations_summary.csv`
- `input/enrichment/string_enrichment_filtered_process_kegg.csv`
- `input/enrichment/enrichment_figures/`
- `input/volcano/Discovery_FP_moi15_seq2_thresh10_k100_default.csv`

Uniqueness scores are computed automatically if the gene CSV lacks
`UniquenessScore`, so the pipeline can start from the raw loading file.

---

## Workflow Dependency Graph

```
Step 1: STRING Enrichment (01_genes_to_string_enrichment.py)
    ↓
Step 2: NCBI Data Fetching (02_fetch_ncbi_data.py) ←── Required for prompts
    ↓
Step 3: Batch Submission (03_submit_and_monitor_batch.py)
    ↓
Step 4: Parse Results (04_parse_and_summarize.py)
    ↓
Step 5: Generate HTML Report (05_generate_html_report.py)
```

---

## Quick Start (100 Programs)

```bash
# Step 1: STRING enrichment (accepts RowID or program_id)
python pipeline/01_genes_to_string_enrichment.py all \
  --input input/genes/FB_moi15_seq2_loading_gene_k100_top300.csv \
  --n-top 300 --species 10090 \
  --json-out results/output/genes_top300.json \
  --csv-out results/output/genes_overview.csv \
  --out-csv-filtered results/output/enrichment_filtered.csv \
  --figures-dir results/output/enrichment_figures

# Step 2: Fetch literature context (REQUIRED before batch submission)
python pipeline/02_fetch_ncbi_data.py \
  --input input/genes/FB_moi15_seq2_loading_gene_k100_top300.csv \
  --csv-out results/output/ncbi_context.csv \
  --json-out results/output/ncbi_context.json \
  --api-key "$NCBI_API_KEY"

# Optional: use Harmonizome gene descriptions instead of NCBI summaries
python pipeline/02_fetch_ncbi_data.py \
  --input input/genes/FB_moi15_seq2_loading_gene_k100_top300.csv \
  --csv-out results/output/ncbi_context.csv \
  --json-out results/output/ncbi_context.json \
  --gene-summary-source harmonizome

# Step 3: Prepare batch & submit to Vertex AI
python pipeline/03_submit_and_monitor_batch.py prepare \
  --gene-file input/genes/FB_moi15_seq2_loading_gene_k100_top300.csv \
  --enrichment-file input/enrichment/string_enrichment_filtered_process_kegg.csv \
  --ncbi-file results/output/ncbi_context.json \
  --regulator-file input/regulators/sceptre_discovery_analysis_results.csv \
  --celltype-dir input/celltype \
  --output-file results/output/batch_request.json

python pipeline/03_submit_and_monitor_batch.py submit-vertex \
  results/output/batch_request.json

# Step 4: Parse results
python pipeline/04_parse_and_summarize.py \
  --gcs-prefix gs://your-bucket/batch-results/ \
  --markdown-dir results/output/annotations/ \
  --summary-csv results/output/annotations/summary.csv

# Step 5: Generate HTML report
python pipeline/05_generate_html_report.py \
  --summary-csv results/output/annotations/summary.csv \
  --annotations-dir results/output/annotations \
  --enrichment-dir input/enrichment/enrichment_figures \
  --volcano-csv input/volcano/Discovery_FP_moi15_seq2_thresh10_k100_default.csv \
  --gene-loading-csv input/genes/FB_moi15_seq2_loading_gene_k100_top300.csv \
  --output-html results/output/annotations/report.html
```

---

## Minimal Enrichment (Single Input CSV)

If you only provide the gene loading CSV, Step 1 now defaults its outputs to
`input/enrichment/` so you can regenerate enrichment inputs in one command:

```bash
python pipeline/01_genes_to_string_enrichment.py all \
  --input input/genes/FB_moi15_seq2_loading_gene_k100_top300.csv \
  --n-top 300
```

Defaults produced:
- `input/enrichment/genes_top300.json`
- `input/enrichment/genes_overview_top300.csv`
- `input/enrichment/string_enrichment_full.csv`
- `input/enrichment/string_enrichment_filtered_process_kegg.csv`

Optional figures:
```bash
python pipeline/01_genes_to_string_enrichment.py all \
  --input input/genes/FB_moi15_seq2_loading_gene_k100_top300.csv \
  --n-top 300 \
  --figures-dir input/enrichment/enrichment_figures
```

---

## Faster Re-runs (Cache/Resume/Figures-Only)

You can avoid re-calling STRING by caching per-program JSON responses and/or
resuming from existing output CSVs:

```bash
python pipeline/01_genes_to_string_enrichment.py enrich \
  --genes-json input/enrichment/genes_top300.json \
  --out-csv-full input/enrichment/string_enrichment_full.csv \
  --out-csv-filtered input/enrichment/string_enrichment_filtered_process_kegg.csv \
  --cache-dir input/enrichment/string_cache \
  --resume
```

Only download figures (skip enrichment CSVs):
```bash
python pipeline/01_genes_to_string_enrichment.py enrich \
  --genes-json input/enrichment/genes_top300.json \
  --figures-dir input/enrichment/enrichment_figures \
  --figures-only
```

Force a fresh STRING query (ignore cache):
```bash
python pipeline/01_genes_to_string_enrichment.py enrich \
  --genes-json input/enrichment/genes_top300.json \
  --out-csv-full input/enrichment/string_enrichment_full.csv \
  --out-csv-filtered input/enrichment/string_enrichment_filtered_process_kegg.csv \
  --cache-dir input/enrichment/string_cache \
  --force-refresh
```

Species can be overridden via CLI (`--species 9606`) or config (`species: 9606`).

---

## Step 3: Batch Submission Options

### Anthropic Direct API
```bash
python pipeline/03_submit_and_monitor_batch.py submit --batch-json batch.json
python pipeline/03_submit_and_monitor_batch.py monitor --batch-id <ID> --output-file results.jsonl
```

### Vertex AI (Claude on GCP)
```bash
python pipeline/03_submit_and_monitor_batch.py submit-vertex batch.json
python pipeline/03_submit_and_monitor_batch.py check-vertex --job-name <JOB_NAME>
```

---

## Config-Driven Workflow

All steps can read defaults from a config file (YAML or JSON). CLI flags always
override config values. See the template in:

- `configs/example_config.yaml`

Example:

```bash
# Step 1 (all) using config defaults
python pipeline/01_genes_to_string_enrichment.py all \
  --config configs/example_config.yaml

# Step 2 using config defaults
python pipeline/02_fetch_ncbi_data.py \
  --config configs/example_config.yaml

# Step 3 prepare using config defaults
python pipeline/03_submit_and_monitor_batch.py prepare \
  --config configs/example_config.yaml
```

### Test Mode in Config

You can enable test mode and restrict to specific programs:

```yaml
test:
  enabled: true
  topics: [1, 2, 3, 4, 5]
```

Steps that support topic limits (`01`, `02`, `03`) will use these topics unless
the CLI provides `--topics` (CLI wins).

## HTML Report Features

- **Stats at top**: Top-loading genes, unique genes, cell-type enrichment
- **LLM annotations**: Full markdown with functional modules
- **Enrichment figures**: KEGG + Process enrichment panels
- **Interactive volcano plot**: 1:1 ratio, priority gene labels
- **Full-text search**: Searches across all annotation content
- **Keyboard navigation**: Arrow keys (←/→) to switch programs
