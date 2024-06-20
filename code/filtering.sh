bcftools filter -i 'QUAL > 20' && DP > 10' variants.vcf > filtered_variants.vcf
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' filtered_variants.vcf > variants_table.txt
