"""
Unit tests for cell-type enrichment validation and program ID extraction.

Tests the flexible program name parsing and non-fatal validation added
for cell-type enrichment inputs.
"""

import pandas as pd
import pytest
from pathlib import Path
import sys
import importlib.util

# Add pipeline directory to path
repo_root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(repo_root / "pipeline"))

# Import from the 01_genes_to_string_enrichment module
import importlib.util
spec = importlib.util.spec_from_file_location(
    "pipeline_step01",
    repo_root / "pipeline" / "01_genes_to_string_enrichment.py"
)
pipeline_step01 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pipeline_step01)

extract_program_id = pipeline_step01.extract_program_id
validate_celltype_enrichment = pipeline_step01.validate_celltype_enrichment


class TestExtractProgramId:
    """Test flexible program ID extraction from various naming formats."""
    
    def test_program_uppercase_underscore(self):
        """Test 'Program_X' format."""
        assert extract_program_id("Program_1") == 1
        assert extract_program_id("Program_100") == 100
        assert extract_program_id("Program_42") == 42
    
    def test_program_lowercase_underscore(self):
        """Test 'program_X' format."""
        assert extract_program_id("program_1") == 1
        assert extract_program_id("program_100") == 100
        assert extract_program_id("program_42") == 42
    
    def test_topic_uppercase_underscore(self):
        """Test 'Topic_X' format."""
        assert extract_program_id("Topic_1") == 1
        assert extract_program_id("Topic_100") == 100
        assert extract_program_id("Topic_42") == 42
    
    def test_topic_lowercase_underscore(self):
        """Test 'topic_X' format."""
        assert extract_program_id("topic_1") == 1
        assert extract_program_id("topic_100") == 100
        assert extract_program_id("topic_42") == 42
    
    def test_program_no_underscore(self):
        """Test 'ProgramX' format (no underscore)."""
        assert extract_program_id("Program1") == 1
        assert extract_program_id("program100") == 100
        assert extract_program_id("PROGRAM42") == 42
    
    def test_topic_no_underscore(self):
        """Test 'TopicX' format (no underscore)."""
        assert extract_program_id("Topic1") == 1
        assert extract_program_id("topic100") == 100
        assert extract_program_id("TOPIC42") == 42
    
    def test_p_prefix(self):
        """Test 'P_X' and 'PX' formats."""
        assert extract_program_id("P_1") == 1
        assert extract_program_id("p_42") == 42
        assert extract_program_id("P1") == 1
        assert extract_program_id("p42") == 42
    
    def test_x_prefix(self):
        """Test 'X_X' and 'XX' formats (regulator file format)."""
        assert extract_program_id("X1") == 1
        assert extract_program_id("X_1") == 1
        assert extract_program_id("x42") == 42
        assert extract_program_id("x_42") == 42
        assert extract_program_id("X100") == 100
    
    def test_plain_integer(self):
        """Test plain integer (as string or int)."""
        assert extract_program_id("1") == 1
        assert extract_program_id("100") == 100
        assert extract_program_id(42) == 42
        assert extract_program_id("  7  ") == 7  # with whitespace
    
    def test_invalid_formats(self):
        """Test that invalid formats return None."""
        assert extract_program_id("invalid") is None
        assert extract_program_id("Program_X") is None
        assert extract_program_id("Program-1") is None  # hyphen not supported
        assert extract_program_id("") is None
        assert extract_program_id(None) is None
        assert extract_program_id("ABC123XYZ") is None


class TestValidateCelltypeEnrichment:
    """Test cell-type enrichment validation function."""
    
    def test_valid_enrichment_data(self, tmp_path):
        """Test validation passes with valid data (enough programs)."""
        # Create valid DataFrame with enough programs to avoid warnings
        programs = [f'Program_{i}' for i in range(1, 12)]  # 11 programs
        cell_types = ['Neuron', 'Astrocyte', 'Microglia']
        
        rows = []
        for i, prog in enumerate(programs):
            for ct in cell_types:
                rows.append({
                    'cell_type': ct,
                    'program': prog,
                    'log2_fc_in_vs_out': 2.5 + (i % 3) * 0.5,
                    'fdr': 0.001 + (i % 3) * 0.01
                })
        
        df = pd.DataFrame(rows)
        
        file_path = tmp_path / "test_enrichment.csv"
        df.to_csv(file_path, index=False)
        
        # Validation should pass (returns True) - no warnings about program count
        result = validate_celltype_enrichment(df, file_path)
        assert result is True
    
    def test_missing_columns(self, tmp_path):
        """Test validation handles missing required columns."""
        # Missing 'fdr' column
        df = pd.DataFrame({
            'cell_type': ['Neuron', 'Astrocyte'],
            'program': ['Program_1', 'Program_1'],
            'log2_fc_in_vs_out': [2.5, 1.2],
        })
        
        file_path = tmp_path / "test_missing_col.csv"
        df.to_csv(file_path, index=False)
        
        # Validation should return False (warnings logged, no exception)
        result = validate_celltype_enrichment(df, file_path)
        assert result is False
    
    def test_empty_dataframe(self, tmp_path):
        """Test validation handles empty DataFrame."""
        df = pd.DataFrame(columns=['cell_type', 'program', 'log2_fc_in_vs_out', 'fdr'])
        
        file_path = tmp_path / "test_empty.csv"
        df.to_csv(file_path, index=False)
        
        # Validation should return False
        result = validate_celltype_enrichment(df, file_path)
        assert result is False
    
    def test_flexible_program_names(self, tmp_path):
        """Test validation handles various program name formats."""
        # Create enough programs to avoid "very few programs" warning
        programs = ['Program_1', 'program_2', 'Topic_3', 'P_4', '5', 
                    'Program_6', 'program_7', 'Topic_8', 'P_9', '10', '11']
        cell_types = ['Neuron', 'Astrocyte', 'Microglia']
        
        rows = []
        for i, prog in enumerate(programs):
            for ct in cell_types:  # Use all 3 cell types
                rows.append({
                    'cell_type': ct,
                    'program': prog,
                    'log2_fc_in_vs_out': 2.5 + (i % 3) * 0.5,
                    'fdr': 0.001 + (i % 3) * 0.01
                })
        
        df = pd.DataFrame(rows)
        
        file_path = tmp_path / "test_flex_names.csv"
        df.to_csv(file_path, index=False)
        
        # Validation should pass - all formats are supported
        result = validate_celltype_enrichment(df, file_path)
        assert result is True
    
    def test_invalid_program_names(self, tmp_path):
        """Test validation warns about unparseable program names."""
        df = pd.DataFrame({
            'cell_type': ['Neuron', 'Astrocyte', 'Microglia'],
            'program': ['Program_1', 'INVALID_NAME', 'Topic_3'],
            'log2_fc_in_vs_out': [2.5, 1.2, 3.0],
            'fdr': [0.001, 0.01, 0.0001]
        })
        
        file_path = tmp_path / "test_invalid_prog.csv"
        df.to_csv(file_path, index=False)
        
        # Validation should return False (warnings about unparseable names)
        result = validate_celltype_enrichment(df, file_path)
        assert result is False
    
    def test_non_numeric_values(self, tmp_path):
        """Test validation handles non-numeric log2FC and FDR values."""
        df = pd.DataFrame({
            'cell_type': ['Neuron', 'Astrocyte', 'Microglia'],
            'program': ['Program_1', 'Program_1', 'Program_2'],
            'log2_fc_in_vs_out': [2.5, 'invalid', 3.0],
            'fdr': [0.001, 0.01, 'bad_fdr']
        })
        
        file_path = tmp_path / "test_non_numeric.csv"
        df.to_csv(file_path, index=False)
        
        # Validation should return False (warnings about non-numeric values)
        result = validate_celltype_enrichment(df, file_path)
        assert result is False
    
    def test_fdr_out_of_range(self, tmp_path):
        """Test validation warns about FDR values outside [0, 1]."""
        df = pd.DataFrame({
            'cell_type': ['Neuron', 'Astrocyte', 'Microglia'],
            'program': ['Program_1', 'Program_1', 'Program_2'],
            'log2_fc_in_vs_out': [2.5, 1.2, 3.0],
            'fdr': [0.001, 1.5, -0.1]  # out of range
        })
        
        file_path = tmp_path / "test_fdr_range.csv"
        df.to_csv(file_path, index=False)
        
        # Validation should return False (warnings about out-of-range FDR)
        result = validate_celltype_enrichment(df, file_path)
        assert result is False
    
    def test_few_programs(self, tmp_path):
        """Test validation warns about very few programs."""
        df = pd.DataFrame({
            'cell_type': ['Neuron', 'Astrocyte'],
            'program': ['Program_1', 'Program_2'],
            'log2_fc_in_vs_out': [2.5, 1.2],
            'fdr': [0.001, 0.01]
        })
        
        file_path = tmp_path / "test_few_prog.csv"
        df.to_csv(file_path, index=False)
        
        # Validation should return False (warning about few programs)
        result = validate_celltype_enrichment(df, file_path)
        assert result is False


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
