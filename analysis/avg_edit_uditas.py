import pandas as pd

# Read the source Excel file
# Replace 'input.xlsx' with your source Excel file name
df = pd.read_excel('/Users/harshini.muthukumar/Desktop/excel/REQ4509-001_big_results.xlsx')

# Create a new DataFrame with the specified columns
columns=['name', 'guide_id_1', 'guide_id_2', 'WT', 'Large Deletion', 'Inversion']

uditas = df[columns]


# If you need to copy data from the source file to the new file,
# you can map the columns accordingly
# For example:
# new_df['name'] = df['source_name_column']
# new_df['guide_id_1'] = df['source_guide1_column']
# etc.

# Save the new DataFrame to an Excel file
uditas.to_excel('/Users/harshini.muthukumar/Desktop/excel/uditas.xlsx', index=False)
