# Variant calling 
samtools sort -o sorted_reads.bam aligned_reads.bam
samtools index sorted_reads.bam
samtools mpileup -uf reference_genome.fa sorted_reads.bam > output.pileup
bcftools call -mv -Ov -o variants.vcf output.pileup
