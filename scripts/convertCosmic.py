#!/usr/bin/env python

from __future__ import print_function
import argparse
import csv
import sys
import cyvcf2
import gzip


# Parse command line
parser = argparse.ArgumentParser(description='Convert cosmic SVs')
parser.add_argument('-v', '--vcf', metavar='input.vcf.gz', required=True, dest='vcf', help='input VCF file (required)')
parser.add_argument('-b', '--break', metavar='breaks.tsv.gz', required=True, dest='breaks', help='COSMIC breakpoint data (required)')
parser.add_argument('-o', '--out', metavar='out.vcf', required=True, dest='out', help='output VCF file (required)')
args = parser.parse_args()

# Parse header
seq = list()
vcf = cyvcf2.VCF(args.vcf)
for rec in vcf.header_iter():
    d = rec.info()
    if "HeaderType" in d.keys():
        if d['HeaderType'] == "CONTIG":
            seq.append(d['ID'])

# Set VCF INFO fields
reservedkey = set(["GRCh", "Chrom From", "Location From min", "Location From max", "Strand From", "Chrom To", "Location To min", "Location To max", "Strand To", "Mutation Type", "Chromosome:G_Start..G_Stop"])
keydict = dict()
with gzip.open(args.breaks) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        for k in row.keys():
            if k in reservedkey:
                continue
            key = k.replace(" ", "").upper()
            if key in keydict.values():
                print("Duplicate key: ", key, file=sys.stderr)
                quit()
            else:
                keydict[k] = key
                vcf.add_info_to_header({'ID': key, 'Description': k,'Type': 'String', 'Number': '1'})
        w = cyvcf2.Writer(args.out, vcf)
        break
w.close()
vcf.close()

# Append VCF records
with open(args.out, "a") as outf:
    with gzip.open(args.breaks) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            svtype = None
            ct = None
            chr1 = None
            svStart = None
            chr2 = None
            svEnd = None
            if "Chromosome:G_Start..G_Stop" in row.keys():
                # Assume CNA format
                fields = row["Chromosome:G_Start..G_Stop"].split(':')
                chr1 = fields[0]
                if row["GRCh"] == "38":
                    chr1 = "chr" + chr1
                if chr1 in seq:
                    chr2 = chr1
                    coord = fields[1].replace("..","-")
                    pos = coord.split("-")
                    svStart = int(pos[0])
                    svEnd = int(pos[1])
                    svtype = "SCNA"
                    if svStart == 0:
                        svStart = 1
                    if svStart > svEnd:
                        tmpPos = svStart
                        svStart = svEnd
                        svEnd = tmpPos
            else:
                # Assume breakpoint format
                chr1 = row["Chrom From"]
                svStart = int((int(row["Location From min"]) + int(row["Location From max"])) / 2)
                chr2 = row["Chrom To"]
                svEnd = int((int(row["Location To min"]) + int(row["Location To max"])) / 2)
                if row["GRCh"] == "38":
                    chr1 = "chr" + chr1
                    chr2 = "chr" + chr2
                if chr1 in seq:
                    if chr2 in seq:
                        idx1 = seq.index(chr1)
                        idx2 = seq.index(chr2)
                        if ((idx1 == idx2) and (svStart <= svEnd)) or (idx1 > idx2):
                            strand = row["Strand From"] + "to" + row["Strand To"]
                        else:
                            tmpChr = chr1
                            chr1 = chr2
                            chr2 = tmpChr
                            tmpPos = svStart
                            svStart = svEnd
                            svEnd = tmpPos
                            strand = row["Strand To"] + "to" + row["Strand From"]
                        mut = row["Mutation Type"]
                        if idx1 == idx2:                        
                            if mut == "intrachromosomal deletion":
                                svtype = "DEL"
                                ct = "3to5"
                            elif mut == "intrachromosomal insertion":
                                svtype = "INS"
                                ct = "NtoN"
                            elif mut == "intrachromosomal tandem duplication":
                                svtype = "DUP"
                                ct = "5to3"
                            elif mut == "intrachromosomal inversion":
                                svtype = "INV"
                                if (strand == "+to+") or (strand == "+to-"):
                                    ct = "3to3"
                                else:
                                    ct = "5to5"
                            elif mut == "intrachromosomal with inverted orientation":
                                svtype = "INV"
                                if (strand == "+to+") or (strand == "+to-"):
                                    ct = "3to3"
                                else:
                                    ct = "5to5"
                            elif mut == "intrachromosomal with non-inverted orientation":
                                svtype = "DEL"
                                ct = "3to5"
                            elif mut == "Intrachromosomal unknown type":
                                svtype = "DEL"
                                ct = "3to5"
                            else:
                                print("Warning: Unknown intrachromosomal mutation type ", mut, file=sys.stderr)
                        else:
                            if mut == "Interchromosomal unknown type":
                                if strand == "+to+":
                                    ct = "3to5"
                                elif strand == "-to-":
                                    ct = "5to3"
                                elif strand == "+to-":
                                    ct = "3to3"
                                elif strand == "-to+":
                                    ct = "5to5"
                                svtype = "BND"
                            else:
                                print("Warning: Unknown interchromosomal mutation type ", mut, file=sys.stderr)
            if svtype is not None:
                info = "CHR2=" + chr2
                info += ";POS2=" + str(svEnd)
                info += ";SVTYPE=" + svtype
                if ct is not None:
                    info += ";CT=" + ct
                for k in row.keys():
                    if k in reservedkey:
                        continue
                    val = row[k].replace(" ", "")
                    if len(val):
                        info += ";" + keydict[k] + "=" + val
                print(chr1, svStart, ".", "N", "<" + svtype + ">", ".", "PASS", info, sep="\t", file=outf)
