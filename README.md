# sansa
Structural Variant Database Annotation

# Installation

`git clone --recursive https://github.com/dellytools/sansa.git`

`cd sansa/`

`make all`

# Annotation

Download an annotation database. Examples are [gnomAD-SV](https://gnomad.broadinstitute.org/) and then run the annotation.

`sansa annotate -d gnomad_v2.1_sv.sites.vcf.gz input.vcf.gz`

The method generates two output files: `anno.bcf` with annotation SVs augmented by a unique ID (INFO/ANNOID) and `query.txt.gz` with query SVs matched to annotation IDs.

[bcftools](https://github.com/samtools/bcftools) can be used to extract all INFO fields you want as annotation. For instance, let's annotate with the VCF ID and EUR_AF for the European allele frequency in gnomad-SV. Always include INFO/ANNOID as the first column.

`bcftools query -H -f "%INFO/ANNOID %ID %INFO/EUR_AF\n" anno.bcf | sed -e 's/^# //' > anno.txt`

Last is a simple join of query SVs with matched database SVs based on the first column (ANNOID).

`join <(zcat query.txt.gz | sort -k 1b,1) anno.txt > results.txt`
