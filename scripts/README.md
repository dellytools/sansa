# COSMIC SVs

Convert COSMIC breakpoints to delly VCF notation

`python convertCosmic.py -v hg38.vcf.gz -b CosmicBreakpointsExport.tsv.gz -o out.vcf`

Sort, compress and index:

`vcf-sort out.vcf | bgzip > cosmic.vcf.gz && tabix cosmic.vcf.gz`

Some COSMIC SV type encodings are not entirely clear to me. Hence, please use `-n` if you use cosmic data as the SV annotation database.

`sansa annotate -n -d cosmic.vcf.gz input.vcf.gz`
