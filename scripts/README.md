# COSMIC structural variant (SV) breakpoints

Convert COSMIC breakpoints to VCF notation for [delly-sansa](https://github.com/dellytools/sansa) annotation

`python convertCosmic.py -v hg38.vcf.gz -b CosmicBreakpointsExport.tsv.gz -o out.vcf`

Sort, compress and index:

`bcftools sort -O b -o cosmic.bcf out.vcf && bcftools index cosmic.bcf`

Some COSMIC SV type encodings are not entirely clear to me. Hence, please use `-n` if you use cosmic data as the SV annotation database.

`sansa annotate -n -d cosmic.bcf input.vcf.gz`

# COSMIC copy-number alterations (CNAs)

Convert COSMIC CNAs to VCF

`python convertCosmic.py -v hg19.vcf.gz -b CosmicCompleteCNA.tsv.gz -o out.vcf`
