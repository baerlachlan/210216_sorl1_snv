#!/bin/bash

## This script follows the recommended gatk workflow for RNA-seq short variant discovery (SNPs + Indels):
## - https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels-

## Minor adjustments have been made as follows:
## - The workflow assumes raw data is in the unmapped BAM format
##   - As the raw data for this dataset is in FASTQ format, the initial SamToFastq method was skipped
##   - Later in the workflow, the MergeBamAlignment step was also skipped
##     - Unmapped BAM files contain metadata and this step allows for merging of that metadata into the aligned BAMS
##     - This was not needed as raw data was FASTQ which does not contain metadata
## - A bootstrapping approach was taken to define a set of known variants
##   - Our zebrafish are most likely quite different to the reference models due to years of breeding
##   - It is better to therefore generate a set of known variants from our data than to use pre-defined ones in a database (dbSNP , Ensembl)
##   - The reasoning and process for this is located at https://gatk.broadinstitute.org/hc/en-us/articles/360035890531?id=44 (Section 3)

SAMPLES = [
    "SRR11951228", "SRR11951229", "SRR11951230", "SRR11951231", "SRR11951232", "SRR11951233",
    "SRR11951234", "SRR11951235", "SRR11951236", "SRR11951237", "SRR11951238", "SRR11951239",
    "SRR11951240", "SRR11951241", "SRR11951242", "SRR11951243", "SRR11951244", "SRR11951245",
    "SRR11951246", "SRR11951247", "SRR11951248", "SRR11951249", "SRR11951250", "SRR11951251"
]
REF_EXT = ["dict", "fa.fai"]
FQC_DIRS = ["00_rawData", "01_trim", "02_align", "07_recalBases"]
FQC_EXT = ["zip", "html"]
VCF_EXT = ["vcf.gz", "vcf.gz.tbi"]

## Set variables that may change between datasets
REFDIR = "/hpcfs/users/a1647910/refs/ensembl-release-101/danio_rerio/"
READ_LEN = 75

rule all:
	input:
		expand("{DIR}/FastQC/{SAMPLE}_fastqc.{EXT}", DIR = FQC_DIRS, SAMPLE = SAMPLES, EXT = FQC_EXT),
        "02_align/featureCounts/genes.out",
        expand("08_callSnvs/selected/{SAMPLE}.vcf.gz", SAMPLE = SAMPLES),
        expand("11_geneiase/ase/{SAMPLE}.static.pval.tsv", SAMPLE = SAMPLES)

include: "smk/modules/refs.smk"
include: "smk/modules/fastqc_raw.smk"
include: "smk/modules/trim.smk"
include: "smk/modules/align.smk"
include: "smk/modules/featureCounts.smk"
include: "smk/modules/markDuplicates.smk"
include: "smk/modules/splitNCigar.smk"
include: "smk/modules/addRG.smk"
include: "smk/modules/knownSnvs.smk"
include: "smk/modules/recalBases.smk"
include: "smk/modules/callSnvs.smk"
include: "smk/modules/wasp.smk"
include: "smk/modules/aseReadCounter.smk"
include: "smk/modules/geneiase.smk"
