"""
Extract prompts from batch request JSON and save as markdown files.

Usage:
    python tools/topic_annotation_workflow/extract_prompts_to_md.py \
        results/output/topic_annotations/fb_k100_prompt_batch_programs_1_6.json \
        results/output/topic_annotations/prompts_preview
"""

import argparse
import json
from pathlib import Path


def extract_prompts(json_file: Path, output_dir: Path) -> None:
    """Extract prompts from batch request JSON and save as markdown files."""
    output_dir.mkdir(parents=True, exist_ok=True)
    
    with open(json_file, "r", encoding="utf-8") as f:
        data = json.load(f)
    
    requests = data.get("requests", [])
    print(f"Found {len(requests)} prompts in {json_file}")
    
    for req in requests:
        custom_id = req.get("custom_id", "unknown")
        content = req.get("params", {}).get("messages", [{}])[0].get("content", "")
        
        # Extract program ID from custom_id (e.g., "topic_1_annotation" -> "1")
        program_id = custom_id.replace("topic_", "").replace("_annotation", "")
        
        output_file = output_dir / f"program_{program_id}_prompt.md"
        with open(output_file, "w", encoding="utf-8") as f:
            f.write(f"# Program {program_id} Annotation Prompt\n\n")
            f.write(content)
        
        print(f"  Saved: {output_file}")
    
    print(f"\nAll prompts saved to {output_dir}")


def main():
    parser = argparse.ArgumentParser(description="Extract prompts from batch JSON to markdown files")
    parser.add_argument("json_file", type=Path, help="Path to batch request JSON file")
    parser.add_argument("output_dir", type=Path, help="Directory to save markdown files")
    args = parser.parse_args()
    
    extract_prompts(args.json_file, args.output_dir)


if __name__ == "__main__":
    main()
