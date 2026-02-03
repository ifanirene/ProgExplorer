"""
@description
Flexible column name matching utility for handling diverse input file formats.
Maps common column name variants to standardized names used throughout the pipeline.

Key features:
- Case-insensitive matching
- Support for common naming conventions (snake_case, camelCase, etc.)
- Clear error messages when required columns are missing

@dependencies
- pandas

@examples
- mapper = ColumnMapper(df)
- gene_col = mapper.get_column('gene')  # Finds 'Gene', 'gene', 'gene_name', etc.
"""

from typing import Dict, List, Optional, Set
import pandas as pd


class ColumnMapper:
    """
    Maps flexible column names to standardized names.
    
    Supports case-insensitive matching and common variants for:
    - Gene names: Name, Gene, gene_name, GeneName, etc.
    - Scores: Score, score, Loading, loading, etc.
    - Program IDs: program_id, RowID, topic, Topic, etc.
    - Cell types: cell_type, celltype, cluster, Cluster, etc.
    - Log2 fold change: log2_fc, log2FC, log2_fold_change, etc.
    - FDR/p-values: fdr, FDR, p_value, pval, etc.
    """
    
    # Define column aliases (all lowercase for matching)
    ALIASES: Dict[str, List[str]] = {
        'gene': [
            'gene', 'name', 'gene_name', 'genename', 'gene_symbol', 
            'genesymbol', 'symbol', 'gene_id', 'geneid'
        ],
        'score': [
            'score', 'loading', 'weight', 'value', 'loading_score',
            'loadingscore', 'gene_score', 'genescore'
        ],
        'program_id': [
            'program_id', 'programid', 'rowid', 'row_id', 
            'topic', 'topic_id', 'topicid', 'program', 'program_number',
            'topic_number', 'k', 'component', 'factor'
        ],
        'cell_type': [
            'cell_type', 'celltype', 'cluster', 'cell_cluster',
            'cellcluster', 'annotation', 'cell_annotation',
            'cellannotation', 'type', 'celltype_id', 'cluster_id'
        ],
        'log2_fc': [
            'log2_fc', 'log2fc', 'log2_fold_change', 'log2foldchange',
            'log2_fc_in_vs_out', 'log2fc_in_vs_out', 'log2_fold',
            'log2fold', 'lfc', 'l2fc', 'fold_change', 'foldchange',
            'fc', 'log2ratio', 'log2_ratio'
        ],
        'fdr': [
            'fdr', 'fdr_corrected', 'q_value', 'qvalue', 'qval', 'q',
            'padj', 'p_adj', 'p_adjusted', 'adjusted_pvalue', 
            'p_value', 'pvalue', 'pval', 'p'
        ]
    }
    
    def __init__(self, df: pd.DataFrame):
        """
        Initialize mapper with a DataFrame.
        
        Args:
            df: DataFrame to map columns for
        """
        self.df = df
        self.columns_lower = {col.lower(): col for col in df.columns}
        self._mapped: Dict[str, str] = {}
    
    def get_column(self, standard_name: str, required: bool = True) -> Optional[str]:
        """
        Find the actual column name that matches the standard name.
        
        Args:
            standard_name: Standard column name (e.g., 'gene', 'score', 'program_id')
            required: If True, raise ValueError when column not found
            
        Returns:
            Actual column name in the DataFrame, or None if not found and not required
            
        Raises:
            ValueError: If required=True and column not found
        """
        # Return cached mapping if available
        if standard_name in self._mapped:
            return self._mapped[standard_name]
        
        # Get aliases for this standard name
        if standard_name not in self.ALIASES:
            if required:
                raise ValueError(f"Unknown standard column name: {standard_name}")
            return None
        
        aliases = self.ALIASES[standard_name]
        
        # Try to find a matching column (case-insensitive)
        for alias in aliases:
            if alias in self.columns_lower:
                actual_col = self.columns_lower[alias]
                self._mapped[standard_name] = actual_col
                return actual_col
        
        if required:
            raise ValueError(
                f"Could not find column for '{standard_name}'. "
                f"Expected one of: {aliases[:5]}... "
                f"Found columns: {sorted(self.df.columns)}"
            )
        
        return None
    
    def get_columns(self, standard_names: List[str], required: bool = True) -> Dict[str, Optional[str]]:
        """
        Find multiple columns at once.
        
        Args:
            standard_names: List of standard column names
            required: If True, raise ValueError for any missing required columns
            
        Returns:
            Dict mapping standard names to actual column names
        """
        result = {}
        missing = []
        
        for name in standard_names:
            try:
                result[name] = self.get_column(name, required=required)
            except ValueError:
                missing.append(name)
                result[name] = None
        
        if missing and required:
            aliases_str = "\n".join([
                f"  - {name}: {', '.join(self.ALIASES.get(name, [])[:5])}"
                for name in missing
            ])
            raise ValueError(
                f"Missing required columns: {missing}\n"
                f"Expected aliases:\n{aliases_str}\n"
                f"Found columns: {sorted(self.df.columns)}"
            )
        
        return result
    
    def rename_columns(self, standard_names: List[str], inplace: bool = False) -> pd.DataFrame:
        """
        Rename DataFrame columns to standard names.
        
        Args:
            standard_names: List of standard column names to rename
            inplace: If True, modify the DataFrame in place
            
        Returns:
            DataFrame with renamed columns
        """
        column_mapping = {}
        for std_name in standard_names:
            actual_col = self.get_column(std_name, required=False)
            if actual_col:
                column_mapping[actual_col] = std_name
        
        if inplace:
            self.df.rename(columns=column_mapping, inplace=True)
            return self.df
        else:
            return self.df.rename(columns=column_mapping)


def standardize_gene_loading(df: pd.DataFrame) -> pd.DataFrame:
    """
    Standardize a gene loading DataFrame to have columns: Name, Score, program_id.
    
    Args:
        df: Gene loading DataFrame
        
    Returns:
        DataFrame with standardized column names
    """
    mapper = ColumnMapper(df)
    
    # Get the required columns
    cols = mapper.get_columns(['gene', 'score', 'program_id'], required=True)
    
    # Rename to standard names
    rename_map = {
        cols['gene']: 'Name',
        cols['score']: 'Score',
        cols['program_id']: 'program_id'
    }
    
    return df.rename(columns=rename_map)


def standardize_celltype_enrichment(df: pd.DataFrame) -> pd.DataFrame:
    """
    Standardize a cell-type enrichment DataFrame.
    Expected output columns: cell_type, program, log2_fc_in_vs_out, fdr
    
    Args:
        df: Cell-type enrichment DataFrame
        
    Returns:
        DataFrame with standardized column names
    """
    mapper = ColumnMapper(df)
    
    # Get the required columns
    cols = mapper.get_columns(['cell_type', 'program_id', 'log2_fc', 'fdr'], required=True)
    
    # Rename to standard names
    rename_map = {
        cols['cell_type']: 'cell_type',
        cols['program_id']: 'program',
        cols['log2_fc']: 'log2_fc_in_vs_out',
        cols['fdr']: 'fdr'
    }
    
    return df.rename(columns=rename_map)
