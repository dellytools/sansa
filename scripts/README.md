# COSMIC SVs

Convert COSMIC breakpoints to delly VCF notation

`python convertCosmic.py -v hg38.vcf.gz -b CosmicBreakpointsExport.tsv.gz -o out.vcf`

Sort, compress and index:

`bcftools sort -O b -o cosmic.bcf out.vcf && bcftools index cosmic.bcf`

Some COSMIC SV type encodings are not entirely clear to me. Hence, please use `-n` if you use cosmic data as the SV annotation database.

`sansa annotate -n -d cosmic.bcf input.vcf.gz`
