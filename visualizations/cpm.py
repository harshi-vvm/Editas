import pandas as pd
import matplotlib.pyplot as plt



# Load the CPM file
cpm = pd.read_csv("/Users/harshini.muthukumar/Downloads/all_samples.merged_cpm.csv")
cpm.set_index("Geneid",inplace = True)# Contains gene_id, sample, and CPM
print(cpm.head())

# First get the data
target_reads = cpm.loc["ENSG00000130164.13"]

# Now assign column names
target_reads.columns = ["sample_name", "reads"]

print(target_reads)

ax = target_reads.plot(kind='scatter',
            figsize=(10, 8),  # Increased height from 6 to 8
            c='#4575b4',
            edgecolor='blue',
            linewidth=0.5,
            x = 'sample_name',
            y = 'reads',
            s = 100)

plt.ylabel('CPM', fontsize=28)
plt.xlabel('Sample', fontsize=28)
plt.xticks(rotation=45, ha='right', fontsize=14, weight = 'bold')
plt.yticks(fontsize=14, weight = 'bold')

plt.tight_layout()

plt.savefig('/Users/harshini.muthukumar/Desktop/ppt/target_scatter.png')
plt.show()
