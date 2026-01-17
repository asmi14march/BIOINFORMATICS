"""
TCGA Colorectal Cancer Multi-Omics Analysis Pipeline

This module integrates RNA-seq expression data and somatic mutation data from TCGA
colorectal cancer samples to identify candidate driver genes.

Features:
- RNA-seq expression normalization (log2 TPM+1)
- Gene-wise nonsynonymous mutation burden calculation from MAF files
- Dataset merging by gene ID
- Pearson and Spearman correlation analysis
- Statistical significance testing (p<0.05)
- Scatter plots and heatmaps visualization
- Export of ranked candidate driver genes
"""

import pandas as pd
import numpy as np
from scipy import stats
import seaborn as sns
import matplotlib.pyplot as plt
from typing import Tuple, Dict, Optional
import warnings

warnings.filterwarnings('ignore')


class TCGAMultiOmicsPipeline:
    """
    Pipeline for integrating RNA-seq expression and somatic mutation data.
    """
    
    def __init__(self):
        self.expression_data = None
        self.mutation_data = None
        self.merged_data = None
        self.correlation_results = None
        
    def load_expression_data(self, filepath: str, gene_col: str = 'gene_id', 
                            tpm_col: str = 'TPM') -> pd.DataFrame:
        """
        Load RNA-seq expression data from file.
        
        Args:
            filepath: Path to expression data file (CSV/TSV)
            gene_col: Column name containing gene identifiers
            tpm_col: Column name containing TPM values
            
        Returns:
            DataFrame with expression data
        """
        try:
            # Try CSV first
            df = pd.read_csv(filepath)
        except:
            # Try TSV
            df = pd.read_csv(filepath, sep='\t')
        
        print(f"Loaded expression data: {df.shape[0]} genes")
        self.expression_data = df
        return df
    
    def normalize_expression(self, df: Optional[pd.DataFrame] = None, 
                            tpm_col: str = 'TPM') -> pd.DataFrame:
        """
        Normalize expression data using log2(TPM + 1) transformation.
        
        Args:
            df: DataFrame with expression data (uses self.expression_data if None)
            tpm_col: Column name containing TPM values
            
        Returns:
            DataFrame with normalized expression
        """
        if df is None:
            df = self.expression_data.copy()
        else:
            df = df.copy()
        
        # Apply log2(TPM + 1) normalization
        df['log2_TPM'] = np.log2(df[tpm_col] + 1)
        
        print(f"Normalized expression data: mean={df['log2_TPM'].mean():.2f}, "
              f"std={df['log2_TPM'].std():.2f}")
        
        self.expression_data = df
        return df
    
    def load_maf_data(self, filepath: str, gene_col: str = 'Hugo_Symbol',
                     variant_class_col: str = 'Variant_Classification') -> pd.DataFrame:
        """
        Load somatic mutation data from MAF file.
        
        Args:
            filepath: Path to MAF file
            gene_col: Column name containing gene symbols
            variant_class_col: Column name containing variant classifications
            
        Returns:
            DataFrame with mutation data
        """
        # MAF files are typically tab-delimited
        try:
            df = pd.read_csv(filepath, sep='\t', comment='#')
        except:
            # Try CSV as fallback
            try:
                df = pd.read_csv(filepath, comment='#')
            except Exception as e:
                raise ValueError(f"Could not read MAF file: {e}")
        
        print(f"Loaded MAF data: {df.shape[0]} mutations")
        self.mutation_data = df
        return df
    
    def calculate_mutation_burden(self, df: Optional[pd.DataFrame] = None,
                                 gene_col: str = 'Hugo_Symbol',
                                 variant_class_col: str = 'Variant_Classification') -> pd.DataFrame:
        """
        Calculate gene-wise nonsynonymous mutation burden.
        
        Nonsynonymous mutations include:
        - Missense_Mutation
        - Nonsense_Mutation
        - Frame_Shift_Del
        - Frame_Shift_Ins
        - In_Frame_Del
        - In_Frame_Ins
        - Nonstop_Mutation
        
        Args:
            df: DataFrame with MAF data (uses self.mutation_data if None)
            gene_col: Column name containing gene symbols
            variant_class_col: Column name containing variant classifications
            
        Returns:
            DataFrame with mutation burden per gene
        """
        if df is None:
            df = self.mutation_data.copy()
        else:
            df = df.copy()
        
        # Check if gene_col exists in dataframe
        if gene_col not in df.columns:
            raise ValueError(f"Column '{gene_col}' not found in MAF data. Available columns: {df.columns.tolist()}")
        
        # Define nonsynonymous mutation types
        nonsynonymous_types = [
            'Missense_Mutation',
            'Nonsense_Mutation',
            'Frame_Shift_Del',
            'Frame_Shift_Ins',
            'In_Frame_Del',
            'In_Frame_Ins',
            'Nonstop_Mutation'
        ]
        
        # Filter for nonsynonymous mutations
        if variant_class_col in df.columns:
            df_nonsyn = df[df[variant_class_col].isin(nonsynonymous_types)]
            print(f"Filtered {df_nonsyn.shape[0]} nonsynonymous mutations from {df.shape[0]} total mutations")
        else:
            # If no variant classification column, use all mutations
            df_nonsyn = df
            print(f"No variant classification column found, using all {df.shape[0]} mutations")
        
        # Remove rows with missing gene symbols
        df_nonsyn = df_nonsyn.dropna(subset=[gene_col])
        
        # Count mutations per gene
        mutation_burden = df_nonsyn.groupby(gene_col).size().reset_index()
        mutation_burden.columns = ['gene_id', 'mutation_count']
        
        print(f"Calculated mutation burden for {mutation_burden.shape[0]} genes")
        print(f"Mean mutation count: {mutation_burden['mutation_count'].mean():.2f}")
        
        self.mutation_data = mutation_burden
        return mutation_burden
    
    def merge_datasets(self, expression_df: Optional[pd.DataFrame] = None,
                      mutation_df: Optional[pd.DataFrame] = None,
                      gene_col: str = 'gene_id') -> pd.DataFrame:
        """
        Merge expression and mutation data by gene ID.
        
        Args:
            expression_df: DataFrame with expression data
            mutation_df: DataFrame with mutation burden
            gene_col: Column name for gene identifiers
            
        Returns:
            Merged DataFrame
        """
        if expression_df is None:
            expression_df = self.expression_data
        if mutation_df is None:
            mutation_df = self.mutation_data
        
        # Ensure gene_id column exists in both dataframes
        if gene_col not in expression_df.columns:
            # Try to find a similar column
            gene_cols = [col for col in expression_df.columns if 'gene' in col.lower()]
            if gene_cols:
                expression_df = expression_df.rename(columns={gene_cols[0]: gene_col})
        
        if gene_col not in mutation_df.columns:
            gene_cols = [col for col in mutation_df.columns if 'gene' in col.lower()]
            if gene_cols:
                mutation_df = mutation_df.rename(columns={gene_cols[0]: gene_col})
        
        # Merge datasets
        merged = pd.merge(expression_df, mutation_df, on=gene_col, how='inner')
        
        # Fill missing mutation counts with 0
        if 'mutation_count' in merged.columns:
            merged['mutation_count'] = merged['mutation_count'].fillna(0)
        
        print(f"Merged datasets: {merged.shape[0]} genes with both expression and mutation data")
        
        self.merged_data = merged
        return merged
    
    def calculate_correlations(self, df: Optional[pd.DataFrame] = None,
                              expression_col: str = 'log2_TPM',
                              mutation_col: str = 'mutation_count') -> pd.DataFrame:
        """
        Calculate Pearson and Spearman correlations between expression and mutation burden.
        
        Args:
            df: Merged DataFrame (uses self.merged_data if None)
            expression_col: Column name for expression values
            mutation_col: Column name for mutation counts
            
        Returns:
            DataFrame with correlation results per gene
        """
        if df is None:
            df = self.merged_data
        
        # For gene-wise analysis, we need samples as columns
        # If we have a single value per gene, calculate overall correlation
        if expression_col in df.columns and mutation_col in df.columns:
            # Calculate overall correlations
            pearson_r, pearson_p = stats.pearsonr(df[expression_col], df[mutation_col])
            spearman_r, spearman_p = stats.spearmanr(df[expression_col], df[mutation_col])
            
            # Create results dataframe
            results = df[[col for col in ['gene_id', expression_col, mutation_col] 
                         if col in df.columns]].copy()
            results['pearson_r'] = pearson_r
            results['pearson_p'] = pearson_p
            results['spearman_r'] = spearman_r
            results['spearman_p'] = spearman_p
            
            print(f"Overall Pearson correlation: r={pearson_r:.3f}, p={pearson_p:.3e}")
            print(f"Overall Spearman correlation: r={spearman_r:.3f}, p={spearman_p:.3e}")
        else:
            results = df.copy()
        
        self.correlation_results = results
        return results
    
    def identify_significant_genes(self, df: Optional[pd.DataFrame] = None,
                                   p_threshold: float = 0.05,
                                   correlation_col: str = 'pearson_p') -> pd.DataFrame:
        """
        Identify genes with significant correlations (p < threshold).
        
        Args:
            df: DataFrame with correlation results
            p_threshold: P-value threshold for significance
            correlation_col: Column name for p-values
            
        Returns:
            DataFrame with significant genes only
        """
        if df is None:
            df = self.correlation_results
        
        if correlation_col in df.columns:
            significant = df[df[correlation_col] < p_threshold].copy()
            print(f"Found {significant.shape[0]} significant genes (p < {p_threshold})")
        else:
            # If no p-value column, return all genes
            significant = df.copy()
            print(f"No p-value column found, returning all {significant.shape[0]} genes")
        
        return significant
    
    def create_scatter_plot(self, df: Optional[pd.DataFrame] = None,
                           expression_col: str = 'log2_TPM',
                           mutation_col: str = 'mutation_count',
                           output_path: str = 'expression_mutation_scatter.png',
                           figsize: Tuple[int, int] = (10, 8)) -> None:
        """
        Generate scatter plot of expression vs mutation burden.
        
        Args:
            df: DataFrame with merged data
            expression_col: Column name for expression values
            mutation_col: Column name for mutation counts
            output_path: Path to save the plot
            figsize: Figure size (width, height)
        """
        if df is None:
            df = self.merged_data
        
        plt.figure(figsize=figsize)
        
        # Create scatter plot
        plt.scatter(df[expression_col], df[mutation_col], alpha=0.5, s=50)
        
        # Add regression line
        z = np.polyfit(df[expression_col], df[mutation_col], 1)
        p = np.poly1d(z)
        plt.plot(df[expression_col], p(df[expression_col]), 
                "r--", alpha=0.8, linewidth=2, label='Linear fit')
        
        # Calculate and display correlation
        if hasattr(self, 'correlation_results') and self.correlation_results is not None:
            if 'pearson_r' in self.correlation_results.columns:
                r = self.correlation_results['pearson_r'].iloc[0]
                p = self.correlation_results['pearson_p'].iloc[0]
                plt.title(f'Expression vs Mutation Burden\nPearson r={r:.3f}, p={p:.3e}', 
                         fontsize=14, fontweight='bold')
        else:
            plt.title('Expression vs Mutation Burden', fontsize=14, fontweight='bold')
        
        plt.xlabel(f'{expression_col}', fontsize=12)
        plt.ylabel(f'{mutation_col}', fontsize=12)
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Scatter plot saved to: {output_path}")
    
    def create_heatmap(self, df: Optional[pd.DataFrame] = None,
                      top_n: int = 50,
                      expression_col: str = 'log2_TPM',
                      mutation_col: str = 'mutation_count',
                      output_path: str = 'top_genes_heatmap.png',
                      figsize: Tuple[int, int] = (12, 10)) -> None:
        """
        Generate heatmap for top genes based on mutation burden.
        
        Args:
            df: DataFrame with merged data
            top_n: Number of top genes to display
            expression_col: Column name for expression values
            mutation_col: Column name for mutation counts
            output_path: Path to save the plot
            figsize: Figure size (width, height)
        """
        if df is None:
            df = self.merged_data
        
        # Select top N genes by mutation burden
        top_genes = df.nlargest(top_n, mutation_col)
        
        # Prepare data for heatmap
        if 'gene_id' in top_genes.columns:
            heatmap_data = top_genes.set_index('gene_id')[[expression_col, mutation_col]]
        else:
            heatmap_data = top_genes[[expression_col, mutation_col]]
        
        # Normalize data for better visualization
        from sklearn.preprocessing import StandardScaler
        scaler = StandardScaler()
        heatmap_data_normalized = pd.DataFrame(
            scaler.fit_transform(heatmap_data),
            index=heatmap_data.index,
            columns=heatmap_data.columns
        )
        
        plt.figure(figsize=figsize)
        sns.heatmap(heatmap_data_normalized.T, cmap='RdYlBu_r', 
                   cbar_kws={'label': 'Z-score'},
                   yticklabels=['Expression', 'Mutation Burden'],
                   linewidths=0.5)
        plt.title(f'Top {top_n} Genes by Mutation Burden', fontsize=14, fontweight='bold')
        plt.xlabel('Genes', fontsize=12)
        plt.ylabel('Features', fontsize=12)
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Heatmap saved to: {output_path}")
    
    def export_candidate_drivers(self, df: Optional[pd.DataFrame] = None,
                                output_path: str = 'candidate_driver_genes.csv',
                                rank_by: str = 'mutation_count') -> pd.DataFrame:
        """
        Export ranked candidate driver genes to CSV file.
        
        Args:
            df: DataFrame with analysis results
            output_path: Path to save the CSV file
            rank_by: Column name to rank genes by
            
        Returns:
            DataFrame with ranked genes
        """
        if df is None:
            if self.correlation_results is not None:
                df = self.correlation_results
            else:
                df = self.merged_data
        
        # Rank genes
        if rank_by in df.columns:
            ranked_genes = df.sort_values(rank_by, ascending=False)
        else:
            ranked_genes = df
        
        # Add rank column
        ranked_genes = ranked_genes.reset_index(drop=True)
        ranked_genes.insert(0, 'rank', range(1, len(ranked_genes) + 1))
        
        # Export to CSV
        ranked_genes.to_csv(output_path, index=False)
        
        print(f"Exported {ranked_genes.shape[0]} candidate driver genes to: {output_path}")
        print(f"\nTop 10 candidate driver genes:")
        display_cols = [col for col in ['rank', 'gene_id', 'log2_TPM', 'mutation_count', 
                                        'pearson_r', 'pearson_p'] if col in ranked_genes.columns]
        print(ranked_genes[display_cols].head(10).to_string(index=False))
        
        return ranked_genes
    
    def run_pipeline(self, expression_file: str, maf_file: str,
                    output_dir: str = './results',
                    p_threshold: float = 0.05) -> Dict[str, pd.DataFrame]:
        """
        Run the complete multi-omics analysis pipeline.
        
        Args:
            expression_file: Path to expression data file
            maf_file: Path to MAF file
            output_dir: Directory to save results
            p_threshold: P-value threshold for significance
            
        Returns:
            Dictionary containing all analysis results
        """
        import os
        os.makedirs(output_dir, exist_ok=True)
        
        print("="*80)
        print("TCGA Colorectal Cancer Multi-Omics Analysis Pipeline")
        print("="*80)
        
        # Step 1: Load and normalize expression data
        print("\n[Step 1] Loading and normalizing expression data...")
        self.load_expression_data(expression_file)
        self.normalize_expression()
        
        # Step 2: Load and process mutation data
        print("\n[Step 2] Loading and processing mutation data...")
        self.load_maf_data(maf_file)
        self.calculate_mutation_burden()
        
        # Step 3: Merge datasets
        print("\n[Step 3] Merging datasets...")
        self.merge_datasets()
        
        # Step 4: Calculate correlations
        print("\n[Step 4] Calculating correlations...")
        self.calculate_correlations()
        
        # Step 5: Identify significant genes
        print("\n[Step 5] Identifying significant genes...")
        significant_genes = self.identify_significant_genes(p_threshold=p_threshold)
        
        # Step 6: Generate visualizations
        print("\n[Step 6] Generating visualizations...")
        self.create_scatter_plot(output_path=os.path.join(output_dir, 'expression_mutation_scatter.png'))
        self.create_heatmap(output_path=os.path.join(output_dir, 'top_genes_heatmap.png'))
        
        # Step 7: Export results
        print("\n[Step 7] Exporting results...")
        ranked_genes = self.export_candidate_drivers(
            output_path=os.path.join(output_dir, 'candidate_driver_genes.csv')
        )
        
        print("\n" + "="*80)
        print("Pipeline completed successfully!")
        print("="*80)
        
        return {
            'expression_data': self.expression_data,
            'mutation_data': self.mutation_data,
            'merged_data': self.merged_data,
            'correlation_results': self.correlation_results,
            'significant_genes': significant_genes,
            'ranked_genes': ranked_genes
        }


def main():
    """
    Example usage of the TCGA Multi-Omics Pipeline.
    """
    # Initialize pipeline
    pipeline = TCGAMultiOmicsPipeline()
    
    # Example: Run with sample data files
    # Note: Replace these paths with actual TCGA data files
    expression_file = 'tcga_crc_expression.csv'
    maf_file = 'tcga_crc_mutations.maf'
    
    try:
        results = pipeline.run_pipeline(
            expression_file=expression_file,
            maf_file=maf_file,
            output_dir='./results',
            p_threshold=0.05
        )
        print("\nAnalysis complete! Check the 'results' directory for outputs.")
    except FileNotFoundError as e:
        print(f"\nError: {e}")
        print("\nPlease provide TCGA colorectal cancer data files:")
        print("1. Expression data (CSV/TSV) with columns: gene_id, TPM")
        print("2. MAF file with somatic mutations")


if __name__ == '__main__':
    main()
