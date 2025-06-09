import pandas as pd
import matplotlib.pyplot as plt


df = pd.read_excel("/Users/harshini.muthukumar/Desktop/ppt/uditas_1.xlsx")
print(df)

#df.set_index("guide_pair",inplace=True)




# Split x-axis labels into two lines by replacing "/" with "/\n"
df['guide_pair'] = df['guide_pair'].str.replace('/', '/\n')

df.set_index("guide_pair", inplace=True)

color = {
    'deletion': 'steelblue',
    'inversion': 'powderblue',
    'wildtype': 'lightgray'
}

ax = df.plot(kind='bar',
             figsize=(10, 8),
             stacked=True,
             width=0.6,
             edgecolor='black',
             linewidth=0.5,
             color=color)

for c in ax.containers:
    ax.bar_label(c, label_type='center', fmt='%g',  # using %g instead of %.1f
                 fontsize=18, color='black', weight='bold')


ax.set_xbound(-0.5, len(df) - 0.5)

plt.xlabel(" ")
plt.ylabel('% of DNA reads', fontsize=28)
plt.xticks(rotation=0, ha='center', fontsize=14, weight='bold')
plt.yticks(fontsize=14, weight='bold')
ax.legend_.remove()

# Remove legend
# plt.legend(title='Transcript Type', bbox_to_anchor=(1.05, 1), loc='upper left').get_title().set_weight('bold')

plt.tight_layout()
plt.savefig('/Users/harshini.muthukumar/Desktop/ppt/ldlr_efficiency.png')
plt.show()




