import mygene
import pandas as pd
import gseapy as gp
from gseapy import barplot, dotplot
import argparse
import os


def load_and_filter_data(input_file: str,padj_threshold: float = 0.05,
                        log2fc_threshold: float = 1):
    """
    Load and filter the differential expression data.

    Args:
        padj_threshold (float): Adjusted p-value threshold
        log2fc_threshold (float): Absolute log2 fold change threshold
    """
    try:
        df = pd.read_csv(input_file)

        filtered_df = df[
            (df['padj'] < padj_threshold) &
            (df['log2FoldChange'].abs() > log2fc_threshold)
        ]

        # Clean gene IDs
        filtered_df['geneID'] = filtered_df['geneID'].str.split('.').str[0]
        ensembl_ids = filtered_df['geneID'].dropna().unique()

        # Convert to gene symbols
        mg = mygene.MyGeneInfo()
        results = mg.getgenes(ensembl_ids, fields='symbol')
        gene_symbols_list = [
            result['symbol'] for result in results if 'symbol' in result
        ]
        print(f"Number of genes: {len(gene_symbols_list)}")
        print("First 10 genes:", gene_symbols_list[:10])

        return gene_symbols_list

    except Exception as e:
        raise Exception(f"Error converting gene IDs: {str(e)}")


def perform_enrichment(input_file: str,output_dir: str):
        """
        Perform enrichment analysis.

        Args:
            gene_sets (List[str]): List of gene set databases to use
        """




        attempts = 5
        for i in range(attempts):
            print(f"Attempt {i + 1} of {attempts}")
            try:
                # Create output directory if it doesn't exist
                os.makedirs(output_dir, exist_ok=True)
                gene_list = load_and_filter_data(input_file)
                enr = gp.enrichr(
                    gene_list=gene_list,
                    gene_sets=['MSigDB_Hallmark_2020', 'KEGG_2, 'GO_Biological_Process_2021'],
                    organism='human',
                    outdir=None,
                )

                # Show results
                print("\nTop 5 enriched terms:")
                print(enr.results)
                enr.results.to_csv(os.path.join(output_dir, 'enrichment_results.csv'), index=False)

                # Create visualizations only if enrichment was successful
                # Create dotplot visualization
                dot_plot = dotplot(enr.results,
                                   column="Adjusted P-value",
                                   x='Gene_set',
                                   size=10,
                                   top_term=10,
                                   figsize=(3, 5),
                                   title="Results",
                                   xticklabels_rot=45,
                                   show_ring=True,
                                   marker='o')
                dot_plot_path = os.path.join(output_dir, 'dotplot.png')
                dot_plot.figure.savefig(dot_plot_path, bbox_inches='tight', dpi=300)

                # Create and save categorical barplot
                bar_plot = barplot(enr.results,
                                   column="Adjusted P-value",
                                   group='Gene_set',
                                   size=10,
                                   top_term=5,
                                   figsize=(3, 5),
                                   color={'KEGG_2021_Human': 'salmon',
                                          'MSigDB_Hallmark_2020': 'darkblue'})
                categorical_path = os.path.join(output_dir, 'barplot_categorical.png')
                bar_plot.figure.savefig(categorical_path, bbox_inches='tight', dpi=300)

                # If everything was successful, break the loop
                break

            except Exception as e:
                print(f"Error in enrichment analysis: {str(e)}")
                if i == attempts - 1:  # If this was the last attempt
                    raise Exception("Max attempts reached. Exiting.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Perform GSEA analysis on differential expression data.')

    parser.add_argument('--input', '-i',
                        type=str,
                        required=True,
                        help='Path to input CSV file containing DEG data')

    parser.add_argument('--output', '-o',
                        type=str,
                        required=True,
                        help='Directory to save output files')

    parser.add_argument('--padj', '-p',
                        type=float,
                        default=0.05,
                        help='Adjusted p-value threshold (default: 0.05)')

    args = parser.parse_args()

    # Use the parsed arguments
    perform_enrichment(
        input_file=args.input,
        output_dir=args.output,
        )

