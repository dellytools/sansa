# sansa
Structural Variant Database Annotation

# Installation

`git clone --recursive https://github.com/dellytools/sansa.git`

`cd sansa/`

`make all`

# Annotation

Download an annotation database. Examples are [gnomAD-SV](https://gnomad.broadinstitute.org/) and then run the annotation.

`sansa annotate -d gnomad_v2.1_sv.sites.vcf.gz input.vcf`

The method generates two output files: `query.tsv.gz` with query SVs and `anno.tsv.gz` with annotation SVs.


