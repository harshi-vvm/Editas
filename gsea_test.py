import mygene
import pandas as pd
import gseapy as gp
from gseapy import barplot, dotplot
import argparse
import os






df = pd.read_csv('/Users/harshini.muthukumar/Downloads/BCL11A-edited_vs_Mock_DEGs.csv')
filtered_df = df[(df['padj'] < 0.05) & (df['log2FoldChange'].abs() > 1)]

# Clean gene IDs and get symbols
filtered_df['geneID'] = filtered_df['geneID'].str.split('.').str[0]
ensembl_ids = filtered_df['geneID'].dropna().unique()

# Get gene symbols
mg = mygene.MyGeneInfo()
results = mg.getgenes(ensembl_ids, fields='symbol')
gene_symbols_list = [result['symbol'] for result in results if 'symbol' in result]

print(f"Number of genes: {len(gene_symbols_list)}")
print("First 10 genes:", gene_symbols_list[:10])

# Perform enrichment analysis
print("\nPerforming enrichment analysis...")
enr = gp.enrichr(
    gene_list=gene_symbols_list,
    gene_sets=['MSigDB_Hallmark_2020', 'KEGG_2021_Human','GO'],
    organism='human',
    outdir=None
)

# Show results
print("\nTop 5 enriched terms:")
print(enr.results.head())