#!/usr/bin/env python3
"""
@description
Unified pipeline orchestrator for topic annotation workflow.
Runs all 5 steps sequentially with a simplified configuration.

Users only need to specify:
- gene_loading: Path to the gene loading matrix CSV
- celltype_annotations: Path to cell-type annotations summary CSV  
- output_dir: Directory for all outputs (intermediates auto-managed)

Key features:
- Single command to run the entire pipeline
- Auto-generates all intermediate file paths
- Supports --start-from and --stop-after for partial runs
- Test mode with --topics to limit to specific programs
- Resume capability with cached API results

@dependencies
- All step scripts in pipeline/

@examples
- Full run:
  python pipeline/run_pipeline.py --config configs/pipeline_config.yaml

- Partial run (steps 1-3 only):
  python pipeline/run_pipeline.py --config configs/pipeline_config.yaml --stop-after batch_submit

- Resume from step 4:
  python pipeline/run_pipeline.py --config configs/pipeline_config.yaml --start-from parse_results

- Test mode (5 topics only):
  python pipeline/run_pipeline.py \
    --gene-loading input/genes/FB_moi15_seq2_loading_gene_k100_top300.csv \
    --celltype-annotations input/celltype/program_celltype_annotations_summary.csv \
    --output-dir results/output/test_run \
    --topics 5,6,8,11,18
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

from dotenv import load_dotenv

# Load environment variables from .env file (supports NCBI_API_KEY, etc.)
load_dotenv()

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# Pipeline steps in order
STEPS = [
    "string_enrichment",   # Step 1: Extract genes + STRING enrichment
    "literature_fetch",    # Step 2: Fetch NCBI/Harmonizome data
    "batch_prepare",       # Step 3a: Prepare LLM batch request
    "batch_submit",        # Step 3b: Submit to Vertex AI
    "parse_results",       # Step 4: Parse LLM responses
    "html_report",         # Step 5: Generate HTML report
]


@dataclass
class PipelineConfig:
    """Simplified pipeline configuration.
    
    User provides:
    - gene_loading: input gene loading CSV
    - celltype_annotations: cell-type summary CSV
    - regulator_file: SCEPTRE regulator results CSV
    - output_dir: base output directory
    
    All intermediate paths are auto-generated.
    """
    # Required inputs (user must provide these)
    gene_loading: Path
    celltype_annotations: Path
    regulator_file: Path
    output_dir: Path
    
    # Optional settings with defaults
    topics: Optional[List[int]] = None  # None = all topics
    species: int = 10090  # Mouse by default
    context: str = '(endothelial OR endothelium OR "vascular endothelial")'
    n_top_genes: int = 300
    top_loading: int = 20
    top_unique: int = 10
    top_enrichment: int = 3
    genes_per_term: int = 10
    
    # LLM settings
    llm_backend: str = "anthropic"  # "anthropic" (default) or "vertex"
    llm_model: str = "claude-4-sonnet-20250514"  # Anthropic model name
    llm_max_tokens: int = 8192
    llm_wait: bool = True  # Wait for batch completion
    
    # Vertex AI settings (only used if llm_backend="vertex")
    vertex_bucket: str = "gs://perturbseq/batch"
    
    # Gene summary options
    full_summaries: bool = False  # Use full HTML summaries (3x longer)
    
    # Resume/caching
    resume: bool = True
    
    # Auto-generated paths (computed from output_dir)
    _paths: Dict[str, Path] = field(default_factory=dict, repr=False)
    
    def __post_init__(self):
        """Generate all intermediate paths from output_dir."""
        self.gene_loading = Path(self.gene_loading)
        self.celltype_annotations = Path(self.celltype_annotations)
        self.output_dir = Path(self.output_dir)
        if self.regulator_file:
            self.regulator_file = Path(self.regulator_file)
        
        # Create output subdirectories
        out = self.output_dir
        self._paths = {
            # Step 1 outputs
            "genes_json": out / "genes_top.json",
            "genes_overview": out / "genes_overview.csv",
            "enrichment_full": out / "string_enrichment" / "enrichment_full.csv",
            "enrichment_filtered": out / "string_enrichment" / "enrichment_filtered.csv",
            "enrichment_figures": out / "string_enrichment" / "figures",
            "enrichment_cache": out / "string_enrichment" / "cache",
            "gene_loading_with_uniqueness": out / "gene_loading_with_uniqueness.csv",
            
            # Step 2 outputs
            "ncbi_json": out / "literature_context.json",
            "ncbi_csv": out / "literature_context.csv",
            
            # Step 3 outputs
            "batch_request_json": out / "llm_batches" / "batch_request.json",
            "batch_request_jsonl": out / "llm_batches" / "batch_request.jsonl",
            
            # Step 4 outputs
            "annotations_dir": out / "annotations",
            "summary_csv": out / "annotations" / "summary.csv",
            
            # Step 5 outputs
            "report_html": out / "annotations" / "report.html",
        }
    
    def get_path(self, key: str) -> Path:
        """Get an auto-generated path by key."""
        return self._paths[key]
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "PipelineConfig":
        """Create config from a dictionary (e.g., loaded from YAML)."""
        # Handle nested 'input' section
        if "input" in data:
            input_section = data["input"]
            data["gene_loading"] = input_section.get("gene_loading")
            data["celltype_annotations"] = input_section.get("celltype_annotations")
            data["regulator_file"] = input_section.get("regulator_file")
        
        # Filter to only known fields
        known_fields = {
            "gene_loading", "celltype_annotations", "output_dir",
            "regulator_file", "topics", "species", "context",
            "n_top_genes", "top_loading", "top_unique",
            "top_enrichment", "genes_per_term",
            "llm_backend", "llm_model", "llm_max_tokens", "llm_wait",
            "vertex_bucket", "full_summaries", "resume"
        }
        filtered = {k: v for k, v in data.items() if k in known_fields}
        return cls(**filtered)
    
    @classmethod
    def from_yaml(cls, path: Path) -> "PipelineConfig":
        """Load config from a YAML file."""
        try:
            import yaml
        except ImportError:
            raise SystemExit("PyYAML is required. Install with: pip install pyyaml")
        
        with open(path, "r", encoding="utf-8") as f:
            data = yaml.safe_load(f)
        return cls.from_dict(data)


def run_step_1_string_enrichment(config: PipelineConfig) -> bool:
    """Step 1: Extract top genes and run STRING enrichment."""
    logger.info("=" * 60)
    logger.info("STEP 1: STRING Enrichment")
    logger.info("=" * 60)
    
    script = Path(__file__).parent / "01_genes_to_string_enrichment.py"
    
    cmd = [
        sys.executable, str(script), "all",
        "--input", str(config.gene_loading),
        "--n-top", str(config.n_top_genes),
        "--json-out", str(config.get_path("genes_json")),
        "--csv-out", str(config.get_path("genes_overview")),
        "--species", str(config.species),
        "--out-csv-full", str(config.get_path("enrichment_full")),
        "--out-csv-filtered", str(config.get_path("enrichment_filtered")),
        "--figures-dir", str(config.get_path("enrichment_figures")),
        "--cache-dir", str(config.get_path("enrichment_cache")),
        "--gene-loading-out", str(config.get_path("gene_loading_with_uniqueness")),
    ]
    
    if config.topics:
        cmd.extend(["--topics", ",".join(map(str, config.topics))])
    if config.resume:
        cmd.append("--resume")
    
    logger.info(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=False)
    
    if result.returncode != 0:
        logger.error("Step 1 failed with return code %d", result.returncode)
        return False
    
    logger.info("Step 1 completed successfully")
    return True


def run_step_2_literature_fetch(config: PipelineConfig) -> bool:
    """Step 2: Fetch literature evidence from NCBI/Harmonizome."""
    logger.info("=" * 60)
    logger.info("STEP 2: Literature Fetch")
    logger.info("=" * 60)
    
    script = Path(__file__).parent / "02_fetch_ncbi_data.py"
    
    # Use the gene loading with uniqueness if it exists, otherwise original
    gene_input = config.get_path("gene_loading_with_uniqueness")
    if not gene_input.exists():
        gene_input = config.gene_loading
    
    cmd = [
        sys.executable, str(script),
        "--input", str(gene_input),
        "--context", config.context,
        "--json-out", str(config.get_path("ncbi_json")),
        "--csv-out", str(config.get_path("ncbi_csv")),
        "--gene-summary-source", "harmonizome",
        # Fetch summaries for exactly top_loading + top_unique genes
        "--top-loading", str(config.top_loading),
        "--top-unique", str(config.top_unique),
    ]
    
    # Optional: Use full HTML summaries (3x longer)
    if config.full_summaries:
        cmd.append("--full-summaries")
    
    if config.regulator_file:
        cmd.extend(["--regulator-file", str(config.regulator_file)])
    if config.topics:
        cmd.extend(["--topics", ",".join(map(str, config.topics))])
    
    # Pass NCBI API key if available in environment
    ncbi_api_key = os.environ.get("NCBI_API_KEY")
    if ncbi_api_key:
        cmd.extend(["--api-key", ncbi_api_key])
    
    logger.info(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=False)
    
    if result.returncode != 0:
        logger.error("Step 2 failed with return code %d", result.returncode)
        return False
    
    logger.info("Step 2 completed successfully")
    return True


def run_step_3a_batch_prepare(config: PipelineConfig) -> bool:
    """Step 3a: Prepare LLM batch request."""
    logger.info("=" * 60)
    logger.info("STEP 3a: Prepare Batch Request")
    logger.info("=" * 60)
    
    script = Path(__file__).parent / "03_submit_and_monitor_batch.py"
    
    # Use gene loading with uniqueness if available
    gene_input = config.get_path("gene_loading_with_uniqueness")
    if not gene_input.exists():
        gene_input = config.gene_loading
    
    cmd = [
        sys.executable, str(script), "prepare",
        "--gene-file", str(gene_input),
        "--celltype-dir", str(config.celltype_annotations.parent),
        "--enrichment-file", str(config.get_path("enrichment_filtered")),
        "--ncbi-file", str(config.get_path("ncbi_json")),
        "--output-file", str(config.get_path("batch_request_json")),
        "--top-loading", str(config.top_loading),
        "--top-unique", str(config.top_unique),
        "--top-enrichment", str(config.top_enrichment),
        "--genes-per-term", str(config.genes_per_term),
    ]
    
    if config.regulator_file:
        cmd.extend(["--regulator-file", str(config.regulator_file)])
    if config.topics:
        cmd.extend(["--topics", ",".join(map(str, config.topics))])
    
    logger.info(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=False)
    
    if result.returncode != 0:
        logger.error("Step 3a failed with return code %d", result.returncode)
        return False
    
    logger.info("Step 3a completed successfully")
    return True


def run_step_3b_batch_submit(config: PipelineConfig) -> Optional[str]:
    """Step 3b: Submit batch to LLM backend.
    
    Uses Anthropic Batch API by default, or Vertex AI if configured.
    Returns the output path/prefix if successful, None otherwise.
    """
    backend = config.llm_backend
    
    logger.info("=" * 60)
    logger.info(f"STEP 3b: Submit to {backend.upper()}")
    logger.info("=" * 60)
    
    script = Path(__file__).parent / "03_submit_and_monitor_batch.py"
    
    if backend == "vertex":
        # Vertex AI submission
        cmd = [
            sys.executable, str(script), "submit-vertex",
            str(config.get_path("batch_request_json")),
            "--model", config.llm_model.replace("-20250514", ""),  # Vertex uses short names
            "--bucket", config.vertex_bucket,
        ]
    else:
        # Anthropic Batch API (default)
        cmd = [
            sys.executable, str(script), "submit",
            str(config.get_path("batch_request_json")),
            "--model", config.llm_model,
            "--max-tokens", str(config.llm_max_tokens),
        ]
    
    if config.llm_wait:
        cmd.append("--wait")
    
    logger.info(f"Running: {' '.join(cmd)}")
    
    # Capture output to get the GCS output path
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        logger.error("Step 3b failed with return code %d", result.returncode)
        logger.error("stderr: %s", result.stderr)
        return None
    
    # Print stdout for visibility (since we captured it)
    if result.stdout:
        for line in result.stdout.splitlines():
            print(line)
    
    # Parse output to find results location
    import re
    output_path = None
    
    if backend == "vertex":
        # Vertex AI: look for GCS prefix
        # The submit-vertex command prints: "SUCCESS! Results available at: gs://..."
        # Vertex AI creates a subfolder like "prediction-model-<timestamp>/" containing predictions.jsonl
        output_base = None
        
        for line in result.stdout.splitlines():
            if "Results available at:" in line or "Output:" in line:
                match = re.search(r'(gs://[^\s]+)', line)
                if match:
                    output_base = match.group(1).rstrip('/')
                    break
        
        # If we found the base output path, look for prediction-model subfolder
        if output_base:
            try:
                list_cmd = ["gcloud", "storage", "ls", output_base + "/"]
                list_result = subprocess.run(list_cmd, capture_output=True, text=True)
                if list_result.returncode == 0:
                    for gcs_line in list_result.stdout.splitlines():
                        if "prediction-model-" in gcs_line:
                            output_path = gcs_line.rstrip('/')
                            break
                if not output_path:
                    output_path = output_base
            except Exception as e:
                logger.warning(f"Could not list GCS directory: {e}")
                output_path = output_base
    else:
        # Anthropic: look for local results file
        # The submit command prints: "SUCCESS! Results saved to: <path>"
        for line in result.stdout.splitlines():
            if "Results saved to:" in line:
                match = re.search(r'Results saved to:\s*(.+)$', line)
                if match:
                    output_path = match.group(1).strip()
                    break
        
        # Fallback: look for the _results.jsonl file
        if not output_path:
            results_file = config.get_path("batch_request_json").with_name(
                config.get_path("batch_request_json").stem + "_results.jsonl"
            )
            if results_file.exists():
                output_path = str(results_file)
    
    logger.info("Step 3b completed successfully")
    if output_path:
        logger.info(f"Results location: {output_path}")
    
    return output_path


def run_step_4_parse_results(config: PipelineConfig, results_path: Optional[str] = None) -> bool:
    """Step 4: Parse LLM responses and generate summary.
    
    Args:
        config: Pipeline configuration
        results_path: Path to results - either:
            - GCS prefix (gs://...) for Vertex AI
            - Local JSONL file path for Anthropic
    """
    logger.info("=" * 60)
    logger.info("STEP 4: Parse Results")
    logger.info("=" * 60)
    
    script = Path(__file__).parent / "04_parse_and_summarize.py"
    
    # Use gene loading with uniqueness if available
    gene_input = config.get_path("gene_loading_with_uniqueness")
    if not gene_input.exists():
        gene_input = config.gene_loading
    
    cmd = [
        sys.executable, str(script),
        "--markdown-dir", str(config.get_path("annotations_dir")),
        "--summary-csv", str(config.get_path("summary_csv")),
        "--gene-loading-file", str(gene_input),
    ]
    
    if results_path:
        if results_path.startswith("gs://"):
            # Vertex AI: use GCS prefix
            cmd.extend(["--gcs-prefix", results_path])
        else:
            # Anthropic: use local file
            cmd.extend(["--results-jsonl", results_path])
    
    logger.info(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=False)
    
    if result.returncode != 0:
        logger.error("Step 4 failed with return code %d", result.returncode)
        return False
    
    logger.info("Step 4 completed successfully")
    return True


def run_step_5_html_report(config: PipelineConfig) -> bool:
    """Step 5: Generate interactive HTML report."""
    logger.info("=" * 60)
    logger.info("STEP 5: Generate HTML Report")
    logger.info("=" * 60)
    
    script = Path(__file__).parent / "05_generate_html_report.py"
    
    # Use gene loading with uniqueness if available
    gene_input = config.get_path("gene_loading_with_uniqueness")
    if not gene_input.exists():
        gene_input = config.gene_loading
    
    cmd = [
        sys.executable, str(script),
        "--summary-csv", str(config.get_path("summary_csv")),
        "--annotations-dir", str(config.get_path("annotations_dir")),
        "--enrichment-dir", str(config.get_path("enrichment_figures")),
        "--gene-loading-csv", str(gene_input),
        "--output-html", str(config.get_path("report_html")),
    ]
    
    if config.regulator_file:
        cmd.extend(["--volcano-csv", str(config.regulator_file)])
    
    logger.info(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=False)
    
    if result.returncode != 0:
        logger.error("Step 5 failed with return code %d", result.returncode)
        return False
    
    logger.info("Step 5 completed successfully")
    return True


def run_pipeline(
    config: PipelineConfig,
    start_from: Optional[str] = None,
    stop_after: Optional[str] = None,
    gcs_prefix: Optional[str] = None,
) -> bool:
    """Run the full pipeline or a subset of steps.
    
    Args:
        config: Pipeline configuration
        start_from: Step to start from (skips earlier steps)
        stop_after: Step to stop after (skips later steps)
        gcs_prefix: GCS prefix for step 4 (if resuming from batch results)
    
    Returns:
        True if all requested steps completed successfully
    """
    # Validate step names
    if start_from and start_from not in STEPS:
        logger.error(f"Unknown step: {start_from}. Valid steps: {STEPS}")
        return False
    if stop_after and stop_after not in STEPS:
        logger.error(f"Unknown step: {stop_after}. Valid steps: {STEPS}")
        return False
    
    # Determine which steps to run
    start_idx = STEPS.index(start_from) if start_from else 0
    stop_idx = STEPS.index(stop_after) if stop_after else len(STEPS) - 1
    
    steps_to_run = STEPS[start_idx:stop_idx + 1]
    
    logger.info("=" * 60)
    logger.info("PIPELINE ORCHESTRATOR")
    logger.info("=" * 60)
    logger.info(f"Output directory: {config.output_dir}")
    logger.info(f"Steps to run: {steps_to_run}")
    if config.topics:
        logger.info(f"Topics: {config.topics}")
    logger.info("=" * 60)
    
    # Create output directory
    config.output_dir.mkdir(parents=True, exist_ok=True)
    
    # Track GCS prefix for step 4
    batch_gcs_prefix = gcs_prefix
    
    # Run each step
    for step in steps_to_run:
        if step == "string_enrichment":
            if not run_step_1_string_enrichment(config):
                return False
        
        elif step == "literature_fetch":
            if not run_step_2_literature_fetch(config):
                return False
        
        elif step == "batch_prepare":
            if not run_step_3a_batch_prepare(config):
                return False
        
        elif step == "batch_submit":
            result = run_step_3b_batch_submit(config)
            if result is None and not batch_gcs_prefix:
                logger.warning("Batch submission did not return GCS prefix")
                # Continue anyway - user may provide it manually
            else:
                batch_gcs_prefix = result or batch_gcs_prefix
        
        elif step == "parse_results":
            if not run_step_4_parse_results(config, batch_gcs_prefix):
                return False
        
        elif step == "html_report":
            if not run_step_5_html_report(config):
                return False
    
    logger.info("=" * 60)
    logger.info("PIPELINE COMPLETED SUCCESSFULLY")
    logger.info("=" * 60)
    logger.info(f"Report: {config.get_path('report_html')}")
    
    return True


def parse_topics_arg(value: str) -> List[int]:
    """Parse comma-separated topic IDs."""
    return [int(x.strip()) for x in value.split(",") if x.strip()]


def main():
    parser = argparse.ArgumentParser(
        description="Run the topic annotation pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Full run with config file
  python run_pipeline.py --config configs/pipeline_config.yaml
  
  # Full run with CLI arguments
  python run_pipeline.py \\
    --gene-loading input/genes/FB_moi15_seq2_loading_gene_k100_top300.csv \\
    --celltype-annotations input/celltype/program_celltype_annotations_summary.csv \\
    --output-dir results/output/my_run
  
  # Test run with specific topics
  python run_pipeline.py --config configs/pipeline_config.yaml --topics 5,6,8,11,18
  
  # Partial run: prepare batch only (steps 1-3a)
  python run_pipeline.py --config configs/pipeline_config.yaml --stop-after batch_prepare
  
  # Resume from parsing (provide GCS prefix from batch job)
  python run_pipeline.py --config configs/pipeline_config.yaml \\
    --start-from parse_results \\
    --gcs-prefix gs://perturbseq/batch/prediction-model-2024...

Steps:
  string_enrichment  Step 1: Extract genes + STRING enrichment
  literature_fetch   Step 2: Fetch NCBI/Harmonizome data
  batch_prepare      Step 3a: Prepare LLM batch request
  batch_submit       Step 3b: Submit to Vertex AI
  parse_results      Step 4: Parse LLM responses
  html_report        Step 5: Generate HTML report
"""
    )
    
    # Config file
    parser.add_argument(
        "--config", "-c",
        help="Path to YAML config file"
    )
    
    # Required inputs (can be overridden by CLI)
    parser.add_argument(
        "--gene-loading",
        help="Path to gene loading matrix CSV"
    )
    parser.add_argument(
        "--celltype-annotations",
        help="Path to cell-type annotations summary CSV"
    )
    parser.add_argument(
        "--output-dir",
        help="Output directory for all results"
    )
    
    # Optional inputs
    parser.add_argument(
        "--regulator-file",
        help="Path to SCEPTRE regulator results CSV"
    )
    
    # Step control
    parser.add_argument(
        "--start-from",
        choices=STEPS,
        help="Step to start from (skips earlier steps)"
    )
    parser.add_argument(
        "--stop-after",
        choices=STEPS,
        help="Step to stop after (skips later steps)"
    )
    parser.add_argument(
        "--gcs-prefix",
        help="GCS prefix for batch results (for resuming at parse_results)"
    )
    
    # Optional settings
    parser.add_argument(
        "--topics",
        help="Comma-separated list of topic IDs to process (default: all)"
    )
    parser.add_argument(
        "--species",
        type=int,
        default=10090,
        help="NCBI taxonomy ID (default: 10090 for mouse)"
    )
    parser.add_argument(
        "--context",
        help="Literature search context query"
    )
    parser.add_argument(
        "--no-resume",
        action="store_true",
        help="Disable resume/caching (re-query all APIs)"
    )
    parser.add_argument(
        "--no-wait",
        action="store_true",
        help="Don't wait for batch job completion"
    )
    
    args = parser.parse_args()
    
    # Load config from file or build from CLI args
    if args.config:
        config = PipelineConfig.from_yaml(Path(args.config))
    else:
        # Require CLI arguments if no config file
        if not args.gene_loading or not args.celltype_annotations or not args.regulator_file or not args.output_dir:
            parser.error(
                "Either --config or all of (--gene-loading, --celltype-annotations, --regulator-file, --output-dir) are required"
            )
        config = PipelineConfig(
            gene_loading=Path(args.gene_loading),
            celltype_annotations=Path(args.celltype_annotations),
            regulator_file=Path(args.regulator_file),
            output_dir=Path(args.output_dir),
        )
    
    # Apply CLI overrides
    if args.gene_loading:
        config.gene_loading = Path(args.gene_loading)
    if args.celltype_annotations:
        config.celltype_annotations = Path(args.celltype_annotations)
    if args.output_dir:
        config.output_dir = Path(args.output_dir)
        config.__post_init__()  # Regenerate paths
    if args.regulator_file:
        config.regulator_file = Path(args.regulator_file)
    if args.topics:
        config.topics = parse_topics_arg(args.topics)
    if args.species:
        config.species = args.species
    if args.context:
        config.context = args.context
    if args.no_resume:
        config.resume = False
    if args.no_wait:
        config.llm_wait = False
    
    # Validate inputs exist
    if not config.gene_loading.exists():
        logger.error(f"Gene loading file not found: {config.gene_loading}")
        return 1
    if not config.celltype_annotations.exists():
        logger.error(f"Cell-type annotations file not found: {config.celltype_annotations}")
        return 1
    if not config.regulator_file.exists():
        logger.error(f"Regulator file not found: {config.regulator_file}")
        return 1
    
    # Run pipeline
    success = run_pipeline(
        config=config,
        start_from=args.start_from,
        stop_after=args.stop_after,
        gcs_prefix=args.gcs_prefix,
    )
    
    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())
