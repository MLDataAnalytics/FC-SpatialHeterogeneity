import pandas as pd
from scipy.stats import mannwhitneyu

# Function to compare two phenotype groups using Mann-Whitney U Test

def compare_phenotypes(group1, group2):
    stat, p_value = mannwhitneyu(group1, group2)
    return stat, p_value

# Example usage
if __name__ == '__main__':
    # Sample data
    phenotype_A = [1.1, 2.0, 3.5, 2.5, 3.3]
    phenotype_B = [1.9, 2.1, 1.5, 2.8, 3.6]

    statistic, p = compare_phenotypes(phenotype_A, phenotype_B)

    print(f'Statistic: {statistic}, P-value: {p}')