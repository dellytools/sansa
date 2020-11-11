[![C/C++ CI](https://github.com/dellytools/sansa/workflows/C/C++%20CI/badge.svg)](https://github.com/dellytools/sansa/actions)
[![Docker CI](https://github.com/dellytools/sansa/workflows/Docker%20CI/badge.svg)](https://hub.docker.com/r/dellytools/sansa/)
[![GitHub license](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://github.com/dellytools/sansa/blob/master/LICENSE)
[![GitHub Releases](https://img.shields.io/github/release/dellytools/sansa.svg)](https://github.com/dellytools/sansa/releases)

# Sansa

Structural variant (SV) annotation.

## Installation

`git clone --recursive https://github.com/dellytools/sansa.git`

`cd sansa/`

`make all`

## SV annotation

Download an annotation database. Examples are [gnomAD-SV](https://gnomad.broadinstitute.org/) or [1000 Genomes phase 3](https://www.internationalgenome.org/phase-3-structural-variant-dataset) and then run the annotation.

`sansa annotate -d gnomad_v2.1_sv.sites.vcf.gz input.vcf.gz`

The method generates two output files: `anno.bcf` with annotation SVs augmented by a unique ID (INFO/ANNOID) and `query.tsv.gz` with query SVs matched to annotation IDs.

[bcftools](https://github.com/samtools/bcftools) can be used to extract all INFO fields you want as annotation. For instance, let's annotate with the VCF ID and EUR_AF for the European allele frequency in gnomad-SV. Always include INFO/ANNOID as the first column.

`bcftools query -H -f "%INFO/ANNOID\t%ID\t%INFO/EUR_AF\n" anno.bcf | sed -e 's/^# //' > anno.tsv`

Last is a simple join of query SVs with matched database SVs based on the first column (ANNOID).

`join anno.tsv <(zcat query.tsv.gz | sort -k 1b,1) > results.tsv`

## SV annotation parameters

[Sansa](https://github.com/dellytools/sansa) matches SVs based on the absolute difference in breakpoint locations (`-b`) and the size ratio (`-r`) of the smaller SV compared to the larger SV. By default, the SVs need to have their start and end breakpoint within 50bp and differ in size by less than 20% (`-r 0.8`).

`sansa annotate -b 50 -r 0.8 -d gnomad_v2.1_sv.sites.vcf.gz input.vcf.gz`

By default, [sansa](https://github.com/dellytools/sansa) only reports the best matching SVs. You can change the matching strategy to `all` using `-s`.

`sansa annotate -s all -d gnomad_v2.1_sv.sites.vcf.gz input.vcf.gz`

You can also include unmatched query SVs in the output using `-m`.

`sansa annotate -m -d gnomad_v2.1_sv.sites.vcf.gz input.vcf.gz`

By default, SVs are only compared within the same SV type (DELs with DELs, INVs with INVs, and so on). For [delly](https://github.com/dellytools/delly) this comparison is INFO/CT aware. You can deactivate this SV type check using `-n`.

`sansa annotate -n -d gnomad_v2.1_sv.sites.vcf.gz input.vcf.gz`

## Feature/Gene annotation

Based on a distance cutoff (`-t`) [sansa](https://github.com/dellytools/sansa) matches SVs to nearby genes. The gene annotation file can be in [gtf/gff2](https://en.wikipedia.org/wiki/General_feature_format) or [gff3](https://en.wikipedia.org/wiki/General_feature_format) format.

`sansa annotate -g Homo_sapiens.GRCh37.87.gtf.gz input.vcf.gz`

`sansa annotate -i Name -g Homo_sapiens.GRCh37.87.gff3.gz input.vcf.gz`

The output has 2 columns for genes near the SV start breakpoint and genes near the SV end breakpoint. For each gene, the output lists the gene name and in paranthesis the distance (negative values: before SV breakpoint, 0: SV breakpoint within gene, positive values: after SV breakpoint) and the strand of the gene (+/-/*).

You can also use the Ensembl gene id or annotate exons instead of genes.

`sansa annotate -i gene_id -g Homo_sapiens.GRCh37.87.gff3.gz input.vcf.gz`

`sansa annotate -f exon -i exon_id -g Homo_sapiens.GRCh37.87.gff3.gz input.vcf.gz`

Gene and SV annotation can be run in a single command.

`sansa annotate -g Homo_sapiens.GRCh37.87.gtf.gz -d gnomad_v2.1_sv.sites.vcf.gz input.vcf.gz`

## Discovering gene fusion candidates

Using [delly](https://github.com/dellytools/delly) and the `INFO/CT` values one can identify gene fusion candidates. Here is the mapping from gene strand to CT values with classical cancer genomics examples (GRCh37 coordinates).

| chr  | start    | chr2 | end       | svtype | ct   | startfeature | endfeature   |
|------|----------|------|-----------|--------|------|--------------|--------------|
| chrA | posStart | chrA | posEnd    | DEL    | 3to5 | geneA(0;+)   | geneB(0;+)   |
| chrA | posStart | chrA | posEnd    | DEL    | 3to5 | geneC(0;-)   | geneD(0;-)   |
| 21   | 39887792 | 21   | 42869743  | DEL    | 3to5 | ERG(0;-)     | TMPRSS2(0;-) |
|------|----------|------|-----------|--------|------|--------------|--------------|
| chrA | posA     | chrB | posB      | BND    | 3to3 | geneA(0;+)   | geneB(0;-)   |
| chrA | posA     | chrB | posB      | BND    | 3to3 | geneC(0;-)   | geneD(0;+)   |
| chrA | posA     | chrB | posB      | BND    | 3to5 | geneA(0;+)   | geneB(0;+)   |
| chrA | posA     | chrB | posB      | BND    | 3to5 | geneC(0;-)   | geneD(0;-)   |
| 21   | 42867595 | 7    | 14027003  | BND    | 3to5 | TMPRSS2(0;-) | ETV1(0;-)    |
| chrA | posA     | chrB | posB      | BND    | 5to5 | geneA(0;+)   | geneB(0;-)   |
| chrA | posA     | chrB | posB      | BND    | 5to5 | geneC(0;-)   | geneD(0;+)   |
| chrA | posA     | chrB | posB      | BND    | 5to3 | geneA(0;+)   | geneB(0;+)   |
| chrA | posA     | chrB | posB      | BND    | 5to3 | geneC(0;-)   | geneD(0;-)   |
| 21   | 39826990 | 1    | 205637229 | BND    | 5to3 | ERG(0;-)     | SLC45A3(0;-) |


## Citation

Tobias Rausch, Thomas Zichner, Andreas Schlattl, Adrian M. Stuetz, Vladimir Benes, Jan O. Korbel.      
DELLY: structural variant discovery by integrated paired-end and split-read analysis.     
Bioinformatics. 2012 Sep 15;28(18):i333-i339.       
[https://doi.org/10.1093/bioinformatics/bts378](https://doi.org/10.1093/bioinformatics/bts378)

## License

Sansa is distributed under the BSD 3-Clause license. Consult the accompanying [LICENSE](https://github.com/dellytools/sansa/blob/master/LICENSE) file for more details.


