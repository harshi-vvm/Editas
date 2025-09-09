import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load DESeq2 results
df = pd.read_excel("/Users/harshini.muthukumar/Downloads/target_Pair5_vs_Mock.DGEs.xlsx")
df = df.dropna(subset=["pvalue"])

# Set thresholds
logfc_thresh = 1
pval_thresh = 0.01

# Filter up/down
up = df[(df["log2FoldChange"] >= logfc_thresh) & (df["pvalue"] < pval_thresh)]
down = df[(df["log2FoldChange"] <= -logfc_thresh) & (df["pvalue"] < pval_thresh)]

# Top 5 genes
top_up = up.sort_values("log2FoldChange", ascending=False).head(5)
top_down = down.sort_values("log2FoldChange").head(5)

# Combine
top_genes = pd.concat([top_down, top_up]).copy()

# Compute -log10(pvalue)
top_genes["neg_log10_pvalue"] = -np.log10(top_genes["pvalue"])

# Sort gene names for plot order
top_genes["geneSymbol"] = pd.Categorical(
    top_genes["geneSymbol"],
    categories=top_genes.sort_values("log2FoldChange")["geneSymbol"],
    ordered=True
)

# Set up color mapping
norm = plt.Normalize(top_genes["neg_log10_pvalue"].min(), top_genes["neg_log10_pvalue"].max())
cmap = plt.cm.viridis
colors = cmap(norm(top_genes["neg_log10_pvalue"].values))
color_list = [tuple(c) for c in colors]  # convert to list of tuples (acceptable for seaborn)

# Create plot and axis
fig, ax = plt.subplots(figsize=(10, 8))

# Barplot
sns.barplot(
    data=top_genes,
    x="log2FoldChange",
    y="geneSymbol",
    palette=color_list,
    ax=ax
)

# Reference line
ax.axvline(0, color="black", linestyle="--")

# Customize labels
ax.set_xlabel("Log2 Fold Change", fontsize=14)
ax.set_ylabel("Gene Symbol", fontsize=14, fontweight="bold")
ax.set_title("Top 5 Up/Down Genes Colored by -log10(p-value)", fontsize=16, fontweight="bold")

# Style gene labels
for label in ax.get_yticklabels():
    label.set_fontweight('bold')
    label.set_fontsize(14)

# Colorbar
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = fig.colorbar(sm, ax=ax)
cbar.set_label("-log10(p-value)", fontsize=12)

plt.tight_layout()
plt.show()
