"""
Advanced Examples for TCGA Multi-Omics Pipeline

This script demonstrates advanced usage patterns and customization options
for the TCGA colorectal cancer multi-omics analysis pipeline.
"""

from tcga_multiomics_pipeline import TCGAMultiOmicsPipeline
import pandas as pd
import os


def example_1_basic_pipeline():
    """
    Example 1: Basic pipeline execution with default parameters.
    """
    print("\n" + "="*80)
    print("Example 1: Basic Pipeline Execution")
    print("="*80)
    
    pipeline = TCGAMultiOmicsPipeline()
    
    results = pipeline.run_pipeline(
        expression_file='./demo_data/sample_expression.csv',
        maf_file='./demo_data/sample_mutations.maf',
        output_dir='./example1_results',
        p_threshold=0.05
    )
    
    print(f"\n✓ Analysis complete!")
    print(f"  - Total genes: {len(results['merged_data'])}")
    print(f"  - Significant genes: {len(results['significant_genes'])}")


def example_2_step_by_step():
    """
    Example 2: Step-by-step analysis with custom parameters.
    """
    print("\n" + "="*80)
    print("Example 2: Step-by-Step Analysis with Custom Parameters")
    print("="*80)
    
    pipeline = TCGAMultiOmicsPipeline()
    
    # Step 1: Load and normalize expression data
    print("\n[Step 1] Loading expression data...")
    pipeline.load_expression_data('./demo_data/sample_expression.csv')
    pipeline.normalize_expression(tpm_col='TPM')
    
    # Step 2: Load and process mutation data
    print("\n[Step 2] Processing mutation data...")
    pipeline.load_maf_data('./demo_data/sample_mutations.maf')
    pipeline.calculate_mutation_burden(
        gene_col='Hugo_Symbol',
        variant_class_col='Variant_Classification'
    )
    
    # Step 3: Merge datasets
    print("\n[Step 3] Merging datasets...")
    merged_data = pipeline.merge_datasets()
    
    # Step 4: Calculate correlations
    print("\n[Step 4] Calculating correlations...")
    pipeline.calculate_correlations(
        expression_col='log2_TPM',
        mutation_col='mutation_count'
    )
    
    # Step 5: Filter for highly mutated genes
    print("\n[Step 5] Filtering for highly mutated genes...")
    highly_mutated = merged_data[merged_data['mutation_count'] >= 5]
    print(f"Found {len(highly_mutated)} genes with ≥5 mutations")
    
    # Step 6: Create visualizations
    print("\n[Step 6] Creating visualizations...")
    os.makedirs('./example2_results', exist_ok=True)
    pipeline.create_scatter_plot(
        df=highly_mutated,
        output_path='./example2_results/highly_mutated_scatter.png'
    )
    pipeline.create_heatmap(
        df=highly_mutated,
        top_n=30,
        output_path='./example2_results/top30_heatmap.png'
    )
    
    # Step 7: Export results
    print("\n[Step 7] Exporting results...")
    pipeline.export_candidate_drivers(
        df=highly_mutated,
        output_path='./example2_results/highly_mutated_genes.csv',
        rank_by='mutation_count'
    )
    
    print(f"\n✓ Custom analysis complete!")


def example_3_multiple_thresholds():
    """
    Example 3: Analyzing significance at multiple p-value thresholds.
    """
    print("\n" + "="*80)
    print("Example 3: Multiple Significance Thresholds")
    print("="*80)
    
    pipeline = TCGAMultiOmicsPipeline()
    
    # Run pipeline
    pipeline.load_expression_data('./demo_data/sample_expression.csv')
    pipeline.normalize_expression()
    pipeline.load_maf_data('./demo_data/sample_mutations.maf')
    pipeline.calculate_mutation_burden()
    pipeline.merge_datasets()
    pipeline.calculate_correlations()
    
    # Test different thresholds
    thresholds = [0.001, 0.01, 0.05, 0.10]
    
    print("\nSignificance threshold analysis:")
    print("-" * 40)
    for threshold in thresholds:
        sig_genes = pipeline.identify_significant_genes(p_threshold=threshold)
        print(f"p < {threshold:5.3f}: {len(sig_genes):4d} significant genes")
    
    print(f"\n✓ Threshold analysis complete!")


def example_4_custom_gene_list():
    """
    Example 4: Analyzing specific genes of interest (e.g., known cancer genes).
    """
    print("\n" + "="*80)
    print("Example 4: Analyzing Specific Genes of Interest")
    print("="*80)
    
    pipeline = TCGAMultiOmicsPipeline()
    
    # Run pipeline
    pipeline.load_expression_data('./demo_data/sample_expression.csv')
    pipeline.normalize_expression()
    pipeline.load_maf_data('./demo_data/sample_mutations.maf')
    pipeline.calculate_mutation_burden()
    merged_data = pipeline.merge_datasets()
    
    # Define genes of interest (in real analysis, these would be known cancer genes)
    # For demo, we'll use the top mutated genes
    genes_of_interest = merged_data.nlargest(20, 'mutation_count')['gene_id'].tolist()
    
    print(f"\nAnalyzing {len(genes_of_interest)} genes of interest:")
    
    # Filter for genes of interest
    filtered_data = merged_data[merged_data['gene_id'].isin(genes_of_interest)]
    
    print(f"Found data for {len(filtered_data)} genes")
    
    # Export focused results
    os.makedirs('./example4_results', exist_ok=True)
    pipeline.create_heatmap(
        df=filtered_data,
        top_n=len(filtered_data),
        output_path='./example4_results/genes_of_interest_heatmap.png'
    )
    
    filtered_data.to_csv('./example4_results/genes_of_interest.csv', index=False)
    
    print(f"\n✓ Gene-specific analysis complete!")


def example_5_comparative_analysis():
    """
    Example 5: Comparative statistics summary.
    """
    print("\n" + "="*80)
    print("Example 5: Comprehensive Statistical Summary")
    print("="*80)
    
    pipeline = TCGAMultiOmicsPipeline()
    
    # Run pipeline
    pipeline.load_expression_data('./demo_data/sample_expression.csv')
    pipeline.normalize_expression()
    pipeline.load_maf_data('./demo_data/sample_mutations.maf')
    pipeline.calculate_mutation_burden()
    merged_data = pipeline.merge_datasets()
    pipeline.calculate_correlations()
    
    # Print comprehensive statistics
    print("\nExpression Statistics:")
    print("-" * 40)
    print(f"  Mean log2(TPM+1): {merged_data['log2_TPM'].mean():.3f}")
    print(f"  Median log2(TPM+1): {merged_data['log2_TPM'].median():.3f}")
    print(f"  Std dev: {merged_data['log2_TPM'].std():.3f}")
    print(f"  Range: [{merged_data['log2_TPM'].min():.3f}, {merged_data['log2_TPM'].max():.3f}]")
    
    print("\nMutation Burden Statistics:")
    print("-" * 40)
    print(f"  Mean mutations/gene: {merged_data['mutation_count'].mean():.3f}")
    print(f"  Median mutations/gene: {merged_data['mutation_count'].median():.3f}")
    print(f"  Std dev: {merged_data['mutation_count'].std():.3f}")
    print(f"  Range: [{merged_data['mutation_count'].min():.0f}, {merged_data['mutation_count'].max():.0f}]")
    
    # Percentile analysis
    print("\nMutation Burden Percentiles:")
    print("-" * 40)
    percentiles = [25, 50, 75, 90, 95, 99]
    for p in percentiles:
        value = merged_data['mutation_count'].quantile(p/100)
        print(f"  {p:2d}th percentile: {value:.1f} mutations")
    
    print("\nCorrelation Results:")
    print("-" * 40)
    if 'pearson_r' in pipeline.correlation_results.columns:
        print(f"  Pearson correlation: r={pipeline.correlation_results['pearson_r'].iloc[0]:.4f}, "
              f"p={pipeline.correlation_results['pearson_p'].iloc[0]:.4e}")
        print(f"  Spearman correlation: r={pipeline.correlation_results['spearman_r'].iloc[0]:.4f}, "
              f"p={pipeline.correlation_results['spearman_p'].iloc[0]:.4e}")
    
    print(f"\n✓ Statistical summary complete!")


def main():
    """
    Run all examples sequentially.
    """
    print("\n" + "="*80)
    print("TCGA Multi-Omics Pipeline: Advanced Examples")
    print("="*80)
    print("\nThese examples demonstrate various usage patterns of the pipeline.")
    print("Make sure you have run demo_pipeline.py first to generate sample data!")
    
    # Check if demo data exists
    if not os.path.exists('./demo_data/sample_expression.csv'):
        print("\n⚠ Demo data not found. Running demo_pipeline.py first...")
        import demo_pipeline
        demo_pipeline.main()
    
    # Run examples
    try:
        example_1_basic_pipeline()
        example_2_step_by_step()
        example_3_multiple_thresholds()
        example_4_custom_gene_list()
        example_5_comparative_analysis()
        
        print("\n" + "="*80)
        print("All examples completed successfully!")
        print("="*80)
        print("\nGenerated output directories:")
        print("  - example1_results/: Basic pipeline output")
        print("  - example2_results/: Custom step-by-step analysis")
        print("  - example4_results/: Gene-specific analysis")
        
    except FileNotFoundError as e:
        print(f"\n✗ Error: {e}")
        print("Please ensure demo data exists by running demo_pipeline.py first.")


if __name__ == '__main__':
    main()
