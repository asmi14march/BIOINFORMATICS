"""
Demo script for TCGA Multi-Omics Pipeline with sample data generation.

This script generates synthetic TCGA colorectal cancer data for demonstration
purposes and runs the complete analysis pipeline.
"""

import pandas as pd
import numpy as np
import os
from tcga_multiomics_pipeline import TCGAMultiOmicsPipeline


def generate_sample_expression_data(n_genes: int = 1000, output_file: str = 'sample_expression.csv'):
    """
    Generate sample RNA-seq expression data.
    
    Args:
        n_genes: Number of genes to generate
        output_file: Output file path
    """
    print(f"Generating sample expression data for {n_genes} genes...")
    
    # Generate gene IDs
    gene_ids = [f"GENE_{i:04d}" for i in range(n_genes)]
    
    # Generate TPM values (log-normal distribution to simulate real TPM data)
    tpm_values = np.random.lognormal(mean=2, sigma=2, size=n_genes)
    tpm_values = np.clip(tpm_values, 0, 10000)  # Clip extreme values
    
    # Create DataFrame
    df = pd.DataFrame({
        'gene_id': gene_ids,
        'TPM': tpm_values,
        'gene_name': [f"Gene_{i}" for i in range(n_genes)]
    })
    
    # Save to file
    df.to_csv(output_file, index=False)
    print(f"Sample expression data saved to: {output_file}")
    print(f"TPM range: {df['TPM'].min():.2f} - {df['TPM'].max():.2f}")
    
    return df


def generate_sample_maf_data(n_mutations: int = 5000, n_genes: int = 1000, 
                             output_file: str = 'sample_mutations.maf'):
    """
    Generate sample MAF (Mutation Annotation Format) data.
    
    Args:
        n_mutations: Number of mutations to generate
        n_genes: Number of unique genes
        output_file: Output file path
    """
    print(f"Generating sample MAF data with {n_mutations} mutations...")
    
    # Generate gene symbols
    gene_symbols = [f"GENE_{i:04d}" for i in range(n_genes)]
    
    # Mutation types
    variant_classifications = [
        'Missense_Mutation',
        'Nonsense_Mutation',
        'Frame_Shift_Del',
        'Frame_Shift_Ins',
        'In_Frame_Del',
        'In_Frame_Ins',
        'Nonstop_Mutation',
        'Silent',
        'Splice_Site'
    ]
    
    # Nonsynonymous types are more common (70%)
    nonsynonymous_types = variant_classifications[:7]
    
    # Generate mutations with some genes having more mutations (driver genes)
    # Create a power-law distribution for mutation counts
    gene_mutation_counts = np.random.zipf(a=2, size=n_genes)
    gene_mutation_counts = np.clip(gene_mutation_counts, 1, 50)
    
    mutations = []
    for gene_idx, count in enumerate(gene_mutation_counts):
        if len(mutations) >= n_mutations:
            break
        gene = gene_symbols[gene_idx]
        for _ in range(int(count)):
            if len(mutations) >= n_mutations:
                break
            # 70% nonsynonymous, 30% other
            if np.random.random() < 0.7:
                variant = np.random.choice(nonsynonymous_types)
            else:
                variant = np.random.choice(variant_classifications[7:])
            
            mutations.append({
                'Hugo_Symbol': gene,
                'Variant_Classification': variant,
                'Tumor_Sample_Barcode': f"TCGA-XX-{np.random.randint(1000, 9999)}",
                'Chromosome': np.random.choice([str(i) for i in range(1, 23)] + ['X', 'Y']),
                'Start_Position': np.random.randint(1000000, 100000000),
                'Reference_Allele': np.random.choice(['A', 'C', 'G', 'T']),
                'Tumor_Seq_Allele2': np.random.choice(['A', 'C', 'G', 'T'])
            })
    
    # Create DataFrame
    df = pd.DataFrame(mutations)
    
    # Save to file
    df.to_csv(output_file, sep='\t', index=False)
    print(f"Sample MAF data saved to: {output_file}")
    print(f"Mutations per gene range: {df.groupby('Hugo_Symbol').size().min()} - {df.groupby('Hugo_Symbol').size().max()}")
    
    return df


def main():
    """
    Run the demo pipeline with sample data.
    """
    print("="*80)
    print("TCGA Multi-Omics Pipeline Demo")
    print("="*80)
    print("\nThis demo generates synthetic data and runs the complete analysis pipeline.\n")
    
    # Create data directory
    data_dir = './demo_data'
    os.makedirs(data_dir, exist_ok=True)
    
    # Generate sample data
    print("\n[Data Generation]")
    print("-"*80)
    expression_file = os.path.join(data_dir, 'sample_expression.csv')
    maf_file = os.path.join(data_dir, 'sample_mutations.maf')
    
    generate_sample_expression_data(n_genes=1000, output_file=expression_file)
    print()
    generate_sample_maf_data(n_mutations=5000, n_genes=1000, output_file=maf_file)
    
    # Run the pipeline
    print("\n[Pipeline Execution]")
    print("-"*80)
    pipeline = TCGAMultiOmicsPipeline()
    
    results = pipeline.run_pipeline(
        expression_file=expression_file,
        maf_file=maf_file,
        output_dir='./demo_results',
        p_threshold=0.05
    )
    
    # Display summary statistics
    print("\n[Summary Statistics]")
    print("-"*80)
    print(f"Total genes analyzed: {results['merged_data'].shape[0]}")
    print(f"Significant genes (p<0.05): {results['significant_genes'].shape[0]}")
    print(f"Mean expression (log2 TPM+1): {results['merged_data']['log2_TPM'].mean():.2f}")
    print(f"Mean mutation burden: {results['merged_data']['mutation_count'].mean():.2f}")
    
    print("\n[Output Files]")
    print("-"*80)
    print("Results saved in './demo_results/' directory:")
    print("  - candidate_driver_genes.csv: Ranked list of candidate driver genes")
    print("  - expression_mutation_scatter.png: Scatter plot of expression vs mutations")
    print("  - top_genes_heatmap.png: Heatmap of top mutated genes")
    
    print("\n" + "="*80)
    print("Demo completed successfully!")
    print("="*80)


if __name__ == '__main__':
    main()
