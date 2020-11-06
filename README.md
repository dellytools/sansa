# sansa

Structural variant (SV) annotation using database SVs

# Installation

`git clone --recursive https://github.com/dellytools/sansa.git`

`cd sansa/`

`make all`

# SV annotation

Download an annotation database. Examples are [gnomAD-SV](https://gnomad.broadinstitute.org/) and then run the annotation.

`sansa annotate -d gnomad_v2.1_sv.sites.vcf.gz input.vcf.gz`

The method generates two output files: `anno.bcf` with annotation SVs augmented by a unique ID (INFO/ANNOID) and `query.tsv.gz` with query SVs matched to annotation IDs.

[bcftools](https://github.com/samtools/bcftools) can be used to extract all INFO fields you want as annotation. For instance, let's annotate with the VCF ID and EUR_AF for the European allele frequency in gnomad-SV. Always include INFO/ANNOID as the first column.

`bcftools query -H -f "%INFO/ANNOID\t%ID\t%INFO/EUR_AF\n" anno.bcf | sed -e 's/^# //' > anno.tsv`

Last is a simple join of query SVs with matched database SVs based on the first column (ANNOID).

`join anno.tsv <(zcat query.tsv.gz | sort -k 1b,1) > results.tsv`

# Parameters

[Sansa](https://github.com/dellytools/sansa) matches SVs based on the absolute difference in breakpoint locations (`-b`) and the size ratio (`-r`) of the smaller SV compared to the larger SV. By default, the SVs need to have their start and end breakpoint within 50bp and differ in size by less than 20% (`-r 0.8`).

`sansa annotate -b 50 -r 0.8 -d gnomad_v2.1_sv.sites.vcf.gz input.vcf.gz`

By default, [sansa](https://github.com/dellytools/sansa) only reports the best matching SVs. You can change the matching strategy to `all` using `-s`.

`sansa annotate -s all -d gnomad_v2.1_sv.sites.vcf.gz input.vcf.gz`

You can also include unmatched query SVs in the output using `-m`.

`sansa annotate -m -d gnomad_v2.1_sv.sites.vcf.gz input.vcf.gz`

By default, SVs are only compared within the same SV type (DELs with DELs, INVs with INVs, and so on). For [delly](https://github.com/dellytools/delly) this comparison is INFO/CT aware. You can deactivate this SV type check using `-n`.

`sansa annotate -n -d gnomad_v2.1_sv.sites.vcf.gz input.vcf.gz`

# Citation

Tobias Rausch, Thomas Zichner, Andreas Schlattl, Adrian M. Stuetz, Vladimir Benes, Jan O. Korbel.      
DELLY: structural variant discovery by integrated paired-end and split-read analysis.     
Bioinformatics. 2012 Sep 15;28(18):i333-i339.       
[https://doi.org/10.1093/bioinformatics/bts378](https://doi.org/10.1093/bioinformatics/bts378)

# License

Sansa is distributed under the BSD 3-Clause license. Consult the accompanying [LICENSE](https://github.com/dellytools/sansa/blob/master/LICENSE) file for more details.


