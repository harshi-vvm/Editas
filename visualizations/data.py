import matplotlib.pyplot as plt
import numpy as np

# Data
categories = ['Input Regions', 'Successfully Converted', 'Failed to Convert']
values = [242254, 29019, 426470]

# Calculate percentages
percentages = [100, (29019/242254)*100, (426470/242254)*100]

# Create bar plot
plt.figure(figsize=(10, 6))
bars = plt.bar(categories, values, color=['blue', 'green', 'red'])

# Add value labels on top of each bar
for bar in bars:
    height = bar.get_height()
    plt.text(bar.get_x() + bar.get_width()/2., height,
             f'{int(height):,}\n({height/values[0]*100:.1f}%)',
             ha='center', va='bottom')

# Customize the plot
plt.title('LiftOver Results: Human to Macaque Genome Conversion', pad=20)
plt.ylabel('Number of Regions')
plt.grid(axis='y', linestyle='--', alpha=0.7)

# Adjust layout to prevent label cutoff
plt.tight_layout()
plt.show()