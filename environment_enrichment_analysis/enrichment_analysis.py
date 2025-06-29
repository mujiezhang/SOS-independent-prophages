"""
Function: Calculate environmental distribution enrichment and perform multiple test correction
Input: env-distribution2.tsv (environmental distribution data)
Output: env-distribution-enrichment-adjust.tsv (corrected enrichment results)
"""

import scipy.stats as stats
import numpy as np
from statsmodels.stats.multitest import multipletests

# Define file paths
input_file = r'E:\课题\4.tem-induce\new-round-refseq\figure\figure2\4\chiplot\env-distribution2.tsv'
output_file = r'E:\课题\4.tem-induce\new-round-refseq\figure\figure2\4\chiplot\env-distribution-enrichment-adjust.tsv'

# Global statistical parameters
total = 49346
si_total = 11806
sd_total = 37540

# Calculate global proportions
p_si = si_total / total
p_sd = sd_total / total
print(f"Global proportions - SI: {p_si:.4f}, SD: {p_sd:.4f}")

# Store results and p-values
rows = []
p_values_list = []

# Read input file and calculate enrichment
with open(input_file, 'r') as f1:
    # Skip possible header lines
    next(f1)
    
    for line in f1:
        data = line.strip().split('\t')
        if not data[0]:  # Skip empty lines
            continue
            
        env = data[0]
        sd_count = int(data[1])
        si_count = int(data[2])
        group_total = sd_count + si_count
        
        # Calculate enrichment factors
        ef1 = (sd_count / group_total) / p_sd
        ef2 = (si_count / group_total) / p_si
        
        # Hypergeometric test for overrepresentation
        p_sd_test = stats.hypergeom.sf(sd_count, total, sd_total, group_total)
        p_si_test = stats.hypergeom.sf(si_count, total, si_total, group_total)
        
        # Store results and p-values
        rows.append((env, ef1, ef2, p_sd_test, p_si_test))
        p_values_list.extend([p_sd_test, p_si_test])

# Convert to NumPy array for multiple test correction
p_values_array = np.array(p_values_list)
reject, corrected_p_values, _, _ = multipletests(p_values_array, alpha=0.05, method='fdr_bh')

# Write results to output file
with open(output_file, 'w') as f2:
    # Write header
    f2.write("Environment\tEF_SD\tEF_SI\tP_SD\tP_SI\tAdjP_SD\tAdjP_SI\n")
    
    # Write corrected results
    for idx, row in enumerate(rows):
        env, ef1, ef2, p1, p2 = row
        adj_p1 = corrected_p_values[idx*2]
        adj_p2 = corrected_p_values[idx*2 + 1]
        
        # Format p-values in scientific notation for small values
        f2.write(f"{env}\t{ef1:.6f}\t{ef2:.6f}\t{p1:.6e}\t{p2:.6e}\t{adj_p1:.6e}\t{adj_p2:.6e}\n")

print("Analysis completed. Results saved to:", output_file)
