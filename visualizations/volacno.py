import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Load data
df = pd.read_excel("/Users/harshini.muthukumar/Downloads/target_Pair5_vs_Mock.DGEs.xlsx")
df = df.dropna(subset=["pvalue", "log2FoldChange"])

# Compute -log10(pvalue)
df["neg_log10_pvalue"] = -np.log10(df["pvalue"])

# Set thresholds
logfc_thresh = 1
pval_thresh = 0.01

# Classify regulation
def get_regulation(row):
    if row["log2FoldChange"] >= logfc_thresh and row["pvalue"] < pval_thresh:
        return "Up"
    elif row["log2FoldChange"] <= -logfc_thresh and row["pvalue"] < pval_thresh:
        return "Down"
    else:
        return "Not Significant"

df["regulation"] = df.apply(get_regulation, axis=1)

# Get top annotated genes
top_up = df[df["regulation"] == "Up"].nlargest(5, "neg_log10_pvalue")
top_down = df[df["regulation"] == "Down"].nlargest(5, "neg_log10_pvalue")
top_annotated = pd.concat([top_up, top_down])

# Plot volcano
plt.figure(figsize=(10, 8))
sns.scatterplot(
    data=df,
    x="log2FoldChange",
    y="neg_log10_pvalue",
    hue="regulation",
    palette={"Up": "red", "Down": "blue", "Not Significant": "gray"},
    alpha=0.7,
    edgecolor=None,
    s=50  # standard size for all dots
)

# Add gene symbol annotations with larger, bold font
for _, row in top_annotated.iterrows():
    plt.text(
        row["log2FoldChange"],
        row["neg_log10_pvalue"] + 0.2,
        row["geneSymbol"],
        fontsize=12,          # larger font size
        fontweight='bold',    # bold font
        ha='center'
    )

# Add threshold lines
plt.axhline(-np.log10(pval_thresh), color='black', linestyle='--', linewidth=0.7)
plt.axvline(logfc_thresh, color='black', linestyle='--', linewidth=0.7)
plt.axvline(-logfc_thresh, color='black', linestyle='--', linewidth=0.7)

# Customize axis and title
plt.xlabel("Log2 Fold Change", fontsize=14)
plt.ylabel("-Log10(p-value)", fontsize=14)
plt.title("Volcano Plot with Gene Annotations", fontsize=16, fontweight='bold')

# Customize legend
legend = plt.legend(
    title="Regulation",
    title_fontsize=14,
    fontsize=12,
    loc="upper right",
    frameon=True,
    framealpha=0.9
)
plt.setp(legend.get_title(), weight='bold')

plt.tight_layout()
plt.show()
