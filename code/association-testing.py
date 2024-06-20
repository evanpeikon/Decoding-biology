import pandas as pd
from scipy.stats import chi2_contingency

variants = pd.read_csv('variants_table.txt', sep='\t', header=None, names=['CHROM', 'POS', 'REF', 'ALT'])

phenotypes = pd.read_csv('phenotypes.csv')

merged_data = pd.merge(variants, phenotypes, on='SampleID')

contingency_table = pd.crosstab(merged_data['ALT'], merged_data['Phenotype'])

chi2, p = chi2_contingency(contingency_table)
print(f'Chi-square: {chi2}, p-value: {p}')
