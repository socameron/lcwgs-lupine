# sum_depths.py (simplified example)

import sys
import pandas as pd

# Read in each depth file and store in a list
depth_dfs = [pd.read_csv(depth_file, sep='\t', header=None, names=['chrom', 'pos', f'depth_{i}'])
             for i, depth_file in enumerate(sys.argv[1:], 1)]

# Merge all dataframes on 'chrom' and 'pos'
merged_df = depth_dfs[0]
for df in depth_dfs[1:]:
    merged_df = pd.merge(merged_df, df, on=['chrom', 'pos'])

# Sum the depths across all samples
merged_df['total_depth'] = merged_df.filter(like='depth_').sum(axis=1)

# Output the aggregated depth information
merged_df.to_csv(sys.stdout, sep='\t', index=False, columns=['chrom', 'pos', 'total_depth'])
