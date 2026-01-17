# BIOINFORMATICS
Multi-omics analysis pipeline for colorectal cancer using TCGA data.

## Overview

This repository contains a Python bioinformatics pipeline for integrating RNA-seq expression data and somatic mutation data from TCGA (The Cancer Genome Atlas) colorectal cancer samples. The pipeline identifies candidate driver genes by analyzing the relationship between gene expression and mutation burden.

## Features

- **RNA-seq Expression Normalization**: Log2(TPM+1) transformation
- **Mutation Burden Calculation**: Gene-wise nonsynonymous mutation counts from MAF files
- **Data Integration**: Merge datasets by gene identifiers
- **Statistical Analysis**: Pearson and Spearman correlation calculations
- **Significance Testing**: Identify genes with p-value < 0.05
- **Visualizations**: 
  - Scatter plots of expression vs mutation burden
  - Heatmaps of top mutated genes
- **Export Functionality**: Ranked candidate driver genes in CSV format

## Installation

### Requirements

- Python 3.8 or higher
- pip package manager

### Setup

1. Clone the repository:
```bash
git clone https://github.com/asmi14march/BIOINFORMATICS.git
cd BIOINFORMATICS
```

2. Install dependencies:
```bash
pip install -r requirements.txt
```

## Usage

### Quick Start with Demo Data

Run the demo script to see the pipeline in action with synthetic data:

```bash
python demo_pipeline.py
```

This will:
1. Generate sample expression and mutation data
2. Run the complete analysis pipeline
3. Create visualizations and export results to `./demo_results/`

### Using Real TCGA Data

#### Option 1: Using the Pipeline Class

```python
from tcga_multiomics_pipeline import TCGAMultiOmicsPipeline

# Initialize pipeline
pipeline = TCGAMultiOmicsPipeline()

# Run complete analysis
results = pipeline.run_pipeline(
    expression_file='path/to/tcga_expression.csv',
    maf_file='path/to/tcga_mutations.maf',
    output_dir='./results',
    p_threshold=0.05
)
```

#### Option 2: Step-by-Step Analysis

```python
from tcga_multiomics_pipeline import TCGAMultiOmicsPipeline

# Initialize pipeline
pipeline = TCGAMultiOmicsPipeline()

# Load and normalize expression data
pipeline.load_expression_data('expression_data.csv')
pipeline.normalize_expression()

# Load and process mutation data
pipeline.load_maf_data('mutations.maf')
pipeline.calculate_mutation_burden()

# Merge datasets
pipeline.merge_datasets()

# Calculate correlations
pipeline.calculate_correlations()

# Identify significant genes
significant_genes = pipeline.identify_significant_genes(p_threshold=0.05)

# Generate visualizations
pipeline.create_scatter_plot(output_path='scatter.png')
pipeline.create_heatmap(output_path='heatmap.png')

# Export results
ranked_genes = pipeline.export_candidate_drivers(output_path='drivers.csv')
```

## Input Data Format

### Expression Data (CSV/TSV)

Required columns:
- `gene_id`: Gene identifier (e.g., HUGO symbol or Ensembl ID)
- `TPM`: Transcripts Per Million values

Example:
```
gene_id,TPM,gene_name
TP53,125.5,Tumor protein p53
KRAS,89.3,KRAS proto-oncogene
APC,234.7,APC regulator of WNT signaling
```

### MAF File (Mutation Annotation Format)

Required columns:
- `Hugo_Symbol`: Gene symbol
- `Variant_Classification`: Type of mutation

Standard MAF format from TCGA is supported. The pipeline filters for nonsynonymous mutations:
- Missense_Mutation
- Nonsense_Mutation
- Frame_Shift_Del
- Frame_Shift_Ins
- In_Frame_Del
- In_Frame_Ins
- Nonstop_Mutation

## Output Files

The pipeline generates several output files in the specified output directory:

1. **candidate_driver_genes.csv**: Ranked list of genes with expression and mutation data
   - Columns: rank, gene_id, log2_TPM, mutation_count, pearson_r, pearson_p, spearman_r, spearman_p

2. **expression_mutation_scatter.png**: Scatter plot showing relationship between expression and mutation burden
   - Includes regression line and correlation statistics

3. **top_genes_heatmap.png**: Heatmap of top 50 genes by mutation burden
   - Shows normalized expression and mutation patterns

## Dependencies

- pandas >= 1.5.0
- numpy >= 1.23.0
- scipy >= 1.9.0
- seaborn >= 0.12.0
- matplotlib >= 3.6.0
- scikit-learn (for data normalization in heatmaps)

## Pipeline Workflow

```
1. Load Expression Data
   ↓
2. Normalize Expression (log2 TPM+1)
   ↓
3. Load MAF Data
   ↓
4. Calculate Mutation Burden
   ↓
5. Merge Datasets by Gene ID
   ↓
6. Calculate Correlations (Pearson & Spearman)
   ↓
7. Identify Significant Genes (p < 0.05)
   ↓
8. Generate Visualizations
   ↓
9. Export Ranked Candidate Drivers
```

## Examples

### Analyzing TCGA Colorectal Cancer Data

```python
from tcga_multiomics_pipeline import TCGAMultiOmicsPipeline

# Initialize
pipeline = TCGAMultiOmicsPipeline()

# Run analysis on TCGA-COAD (Colon Adenocarcinoma) data
results = pipeline.run_pipeline(
    expression_file='TCGA_COAD_expression.csv',
    maf_file='TCGA_COAD_mutations.maf',
    output_dir='./coad_results',
    p_threshold=0.05
)

# Access results
print(f"Analyzed {len(results['merged_data'])} genes")
print(f"Found {len(results['significant_genes'])} significant correlations")
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

See LICENSE file for details.

## Citation

If you use this pipeline in your research, please cite:
```
TCGA Multi-Omics Analysis Pipeline for Colorectal Cancer
https://github.com/asmi14march/BIOINFORMATICS
```

## Contact

For questions or issues, please open an issue on GitHub.
