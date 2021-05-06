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
FQC_EXT = ["zip", "html"]
VCF_EXT = ["vcf.gz", "vcf.gz.tbi"]

## Set variables that may change between datasets
REFDIR = "/hpcfs/users/a1647910/refs/ensembl-release-101/danio_rerio/"
READ_LEN = 75

rule all:
	input:
		expand("00_rawData/FastQC/{SAMPLE}_fastqc.{EXT}", SAMPLE = SAMPLES, EXT = FQC_EXT),
		expand("01_trimmedData/FastQC/{SAMPLE}_fastqc.{EXT}", SAMPLE = SAMPLES, EXT = FQC_EXT),
		expand("02_alignedData/FastQC/{SAMPLE}Aligned.sortedByCoord.out_fastqc.{EXT}", SAMPLE = SAMPLES, EXT = FQC_EXT),
		"02_alignedData/featureCounts/genes.out",
		expand("03_markDuplicates/FastQC/{SAMPLE}_fastqc.{EXT}", SAMPLE = SAMPLES, EXT = FQC_EXT),
		expand("12_selectVariants/mergedVcf/mergedVcf.{EXT}", EXT = VCF_EXT),
		expand("13_countAlleles/counts/{SAMPLE}.alleleCounts.tsv", SAMPLE = SAMPLES),
		expand("14_aseReadCounts/{SAMPLE}.tsv", SAMPLE = SAMPLES)

###########################
## Build reference files ##
###########################

rule unzip_refFa:
	input:
		REFDIR + "Danio_rerio.GRCz11.dna.primary_assembly.fa.gz"
	output:
		temp(REFDIR + "Danio_rerio.GRCz11.dna.primary_assembly.fa")
	conda:
		"snakemake/envs/default.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 4000,
		hours = 0,
		mins = 30
	shell:
		"gunzip -c {input} > {output}"

rule star_index:
	input:
		refFa = REFDIR + "Danio_rerio.GRCz11.dna.primary_assembly.fa",
		gtf = REFDIR + "Danio_rerio.GRCz11.101.chr.gtf.gz"
	output:
		temp(directory(REFDIR + "star/"))
	params:
		overhang = READ_LEN-1
	conda:
		"snakemake/envs/default.yaml"
	resources:
		cpu = 16,
		ntasks = 1,
		mem_mb = 32000,
		hours = 0,
		mins = 30
	shell:
		"""
		zcat {input.gtf} > temp.gtf

		STAR \
			--runThreadN {resources.cpu} \
			--runMode genomeGenerate \
			--genomeDir {output} \
			--genomeFastaFiles {input.refFa} \
			--sjdbGTFfile temp.gtf \
			--sjdbOverhang {params.overhang}

		rm temp.gtf
		"""

## Reference dictionary and index needs to be created as described in:
## https://gatk.broadinstitute.org/hc/en-us/articles/360035531652-FASTA-Reference-genome-format
rule ref_dict:
	input:
		REFDIR + "Danio_rerio.GRCz11.dna.primary_assembly.fa"
	output:
		temp(REFDIR + "Danio_rerio.GRCz11.dna.primary_assembly.dict")
	conda:
		"snakemake/envs/default.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 4000,
		hours = 0,
		mins = 10
	shell:
		"gatk CreateSequenceDictionary -R {input}"

rule ref_index:
	input:
		REFDIR + "Danio_rerio.GRCz11.dna.primary_assembly.fa"
	output:
		temp(REFDIR + "Danio_rerio.GRCz11.dna.primary_assembly.fa.fai")
	conda:
		"snakemake/envs/default.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 200,
		hours = 0,
		mins = 10
	shell:
		"samtools faidx {input}"

################################################
## Workflow for counting alleles at SNV sites ##
################################################

rule fastqc_raw:
	input:
		"00_rawData/fastq/{SAMPLE}.fastq.gz"
	output:
		"00_rawData/FastQC/{SAMPLE}_fastqc.zip",
		"00_rawData/FastQC/{SAMPLE}_fastqc.html"
	params:
		outDir = "00_rawData/FastQC/"
	conda:
		"snakemake/envs/default.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 2000,
		hours = 1,
		mins = 0
	shell:
		"fastqc -t {resources.cpu} -o {params.outDir} --noextract {input}"

rule trim:
	input:
		R1 = "00_rawData/fastq/{SAMPLE}.fastq.gz"
	output:
		R1 = temp("01_trimmedData/fastq/{SAMPLE}.fastq.gz"),
		html = "01_trimmedData/log/{SAMPLE}.html"
	conda:
		"snakemake/envs/default.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 2000,
		hours = 2,
		mins = 0
	shell:
		"""
		fastp \
			-i {input.R1}  \
        	-o {output.R1} \
			--qualified_quality_phred 20 \
			--length_required 35 \
			--unqualified_percent_limit 100 \
			--complexity_threshold 100 \
			--cut_front \
			--cut_tail \
			--trim_poly_g \
			--thread 1 \
			--html {output.html} \
			--json /dev/null \
		"""

rule fastqc_trim:
	input:
		"01_trimmedData/fastq/{SAMPLE}.fastq.gz"
	output:
		"01_trimmedData/FastQC/{SAMPLE}_fastqc.zip",
		"01_trimmedData/FastQC/{SAMPLE}_fastqc.html"
	params:
		outDir = "01_trimmedData/FastQC/"
	conda:
		"snakemake/envs/default.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 2000,
		hours = 1,
		mins = 0
	shell:
		"fastqc -t {resources.cpu} -o {params.outDir} --noextract {input}"

rule align:
	input:
		R1 = "01_trimmedData/fastq/{SAMPLE}.fastq.gz",
		starIndex = REFDIR + "star/"
	output:
		bam = temp("02_alignedData/bam/{SAMPLE}Aligned.sortedByCoord.out.bam"),
		bamIndex = temp("02_alignedData/bam/{SAMPLE}Aligned.sortedByCoord.out.bam.bai")
	params:
		overhang = READ_LEN-1,
		bname = "02_alignedData/bam/{SAMPLE}"
	conda:
		"snakemake/envs/default.yaml"
	resources:
		cpu = 16,
		ntasks = 1,
		mem_mb = 32000,
		hours = 5,
		mins = 0
	shell:
		"""
		STAR \
			--genomeDir {input.starIndex}\
			--runThreadN {resources.cpu} \
			--readFilesIn {input.R1} \
			--readFilesCommand "gunzip -c" \
			--sjdbOverhang {params.overhang} \
			--outSAMtype BAM SortedByCoordinate \
			--twopassMode Basic \
			--outFileNamePrefix {params.bname}

		mkdir -p 02_alignedData/log
		mv {params.bname}*out 02_alignedData/log
		mv {params.bname}*tab 02_alignedData/log

		samtools index {output.bam}
		"""

rule fastqc_align:
	input:
		"02_alignedData/bam/{SAMPLE}Aligned.sortedByCoord.out.bam"
	output:
		"02_alignedData/FastQC/{SAMPLE}Aligned.sortedByCoord.out_fastqc.zip",
		"02_alignedData/FastQC/{SAMPLE}Aligned.sortedByCoord.out_fastqc.html"
	params:
		outDir = "02_alignedData/FastQC/"
	conda:
		"snakemake/envs/default.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 2000,
		hours = 1,
		mins = 0
	shell:
		"fastqc -t {resources.cpu} -o {params.outDir} --noextract {input}"

rule featureCounts:
	input:
		bam = expand("02_alignedData/bam/{SAMPLE}Aligned.sortedByCoord.out.bam", SAMPLE = SAMPLES),
		gtf = REFDIR + "Danio_rerio.GRCz11.101.chr.gtf.gz"
	output:
		counts = "02_alignedData/featureCounts/counts.out",
		summary = "02_alignedData/featureCounts/counts.out.summary",
		genes = "02_alignedData/featureCounts/genes.out"
	conda:
		"snakemake/envs/default.yaml"
	resources:
		cpu = 2,
		ntasks = 2,
		mem_mb = 4000,
		hours = 2,
		mins = 0
	shell:
		"""
		featureCounts \
			-Q 10 \
			-s 0 \
			-T 4 \
			-a {input.gtf} \
			-o {output.counts} {input.bam}

		## Storing the output in a single file
		cut -f1,7- {output.counts} | \
		sed 1d > {output.genes}
		"""

rule mark_duplicates:
	input:
		"02_alignedData/bam/{SAMPLE}Aligned.sortedByCoord.out.bam"
	output:
		bam = temp("03_markDuplicates/bam/{SAMPLE}.bam"),
		bamIndex = temp("03_markDuplicates/bam/{SAMPLE}.bai"),
		metrics = "03_markDuplicates/log/{SAMPLE}.metrics"
	conda:
		"snakemake/envs/default.yaml"
	resources:
		cpu = 2,
		ntasks = 1,
		mem_mb = 4000,
		hours = 4,
		mins = 0
	shell:
		"""
		gatk \
 	    	MarkDuplicates \
 	        --INPUT {input} \
 	        --OUTPUT {output.bam}  \
 	        --CREATE_INDEX true \
 	        --VALIDATION_STRINGENCY SILENT \
 	        --METRICS_FILE {output.metrics}
		"""

rule fastqc_duplicates:
	input:
		"03_markDuplicates/bam/{SAMPLE}.bam"
	output:
		"03_markDuplicates/FastQC/{SAMPLE}_fastqc.zip",
		"03_markDuplicates/FastQC/{SAMPLE}_fastqc.html"
	params:
		outDir = "03_markDuplicates/FastQC/"
	conda:
		"snakemake/envs/default.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 2000,
		hours = 1,
		mins = 0
	shell:
		"fastqc -t {resources.cpu} -o {params.outDir} --noextract {input}"

rule splitNCigar:
	input:
		bam = "03_markDuplicates/bam/{SAMPLE}.bam",
		bamIndex = "03_markDuplicates/bam/{SAMPLE}.bai",
		refFa = REFDIR + "Danio_rerio.GRCz11.dna.primary_assembly.fa",
		refIndex = REFDIR + "Danio_rerio.GRCz11.dna.primary_assembly.fa.fai",
		refDict = REFDIR + "Danio_rerio.GRCz11.dna.primary_assembly.dict"
	output:
		bam = temp("04_splitNCigar/bam/{SAMPLE}.bam"),
		bamIndex = temp("04_splitNCigar/bam/{SAMPLE}.bai")
	conda:
		"snakemake/envs/default.yaml"
	resources:
		cpu = 8,
		ntasks = 1,
		mem_mb = 32000,
		hours = 4,
		mins = 0
	shell:
		"""
		gatk \
            SplitNCigarReads \
            -R {input.refFa} \
            -I {input.bam} \
            -O {output.bam}
		"""

## gatk HaplotypeCaller requires read group (RG) tags
## An explanation of this is found at https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups
rule addRG:
	input:
		bam = "04_splitNCigar/bam/{SAMPLE}.bam",
		bamIndex = "04_splitNCigar/bam/{SAMPLE}.bai"
	output:
		bam = temp("05_addRG/bam/{SAMPLE}.bam"),
		bamIndex = temp("05_addRG/bam/{SAMPLE}.bai")
	conda:
		"snakemake/envs/default.yaml"
	resources:
		cpu = 1,
		ntasks = 2,
		mem_mb = 4000,
		hours = 2,
		mins = 0
	shell:
		"""
		gatk \
			AddOrReplaceReadGroups \
    		-I {input.bam} \
   			-O {output.bam} \
    		-SORT_ORDER coordinate \
    		-RGID default \
    		-RGLB default \
    		-RGSM {wildcards.SAMPLE} \
			-RGPU default \
    		-RGPL ILLUMINA \
    		-CREATE_INDEX True
		"""

rule callVariants_noRecal:
	input:
		bam = "05_addRG/bam/{SAMPLE}.bam",
		bamIndex = "05_addRG/bam/{SAMPLE}.bai",
		refFa = REFDIR + "Danio_rerio.GRCz11.dna.primary_assembly.fa",
		refIndex = REFDIR + "Danio_rerio.GRCz11.dna.primary_assembly.fa.fai",
		refDict = REFDIR + "Danio_rerio.GRCz11.dna.primary_assembly.dict"
	output:
		vcf = temp("06_callVariants_noRecal/vcf/{SAMPLE}.vcf.gz"),
		vcfIndex = temp("06_callVariants_noRecal/vcf/{SAMPLE}.vcf.gz.tbi")
	conda:
		"snakemake/envs/default.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 8000,
		hours = 72,
		mins = 0
	shell:
		"""
		gatk --java-options "-Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
		HaplotypeCaller \
		-R {input.refFa} \
		-I {input.bam} \
		-O {output.vcf} \
		-dont-use-soft-clipped-bases \
		--standard-min-confidence-threshold-for-calling 20
		"""

rule knownVariants:
	input:
		vcf = "06_callVariants_noRecal/vcf/{SAMPLE}.vcf.gz",
		vcfIndex = "06_callVariants_noRecal/vcf/{SAMPLE}.vcf.gz.tbi",
		refFa = REFDIR + "Danio_rerio.GRCz11.dna.primary_assembly.fa",
		refIndex = REFDIR + "Danio_rerio.GRCz11.dna.primary_assembly.fa.fai",
		refDict = REFDIR + "Danio_rerio.GRCz11.dna.primary_assembly.dict"
	output:
		vcf = temp("07_knownVariants/vcf/{SAMPLE}.vcf.gz"),
		vcfIndex = temp("07_knownVariants/vcf/{SAMPLE}.vcf.gz.tbi")
	conda:
		"snakemake/envs/default.yaml"
	resources:
		cpu = 1,
		ntasks = 2,
		mem_mb = 1000,
		hours = 2,
		mins = 0
	shell:
		"""
		gatk \
		    VariantFiltration \
			--R {input.refFa} \
			--V {input.vcf} \
			--window 35 \
			--cluster 3 \
			--filter-name "FS" \
			--filter "FS > 30.0" \
			--filter-name "QD" \
			--filter "QD < 2.0" \
			-O {output.vcf}
		"""

## We now follow the GATK workflow as normal
rule baseRecal:
	input:
		bam = "05_addRG/bam/{SAMPLE}.bam",
		bamIndex = "05_addRG/bam/{SAMPLE}.bai",
		refFa = REFDIR + "Danio_rerio.GRCz11.dna.primary_assembly.fa",
		refIndex = REFDIR + "Danio_rerio.GRCz11.dna.primary_assembly.fa.fai",
		refDict = REFDIR + "Danio_rerio.GRCz11.dna.primary_assembly.dict",
		dbsnp = "07_knownVariants/vcf/{SAMPLE}.vcf.gz",
		dbsnpIndex = "07_knownVariants/vcf/{SAMPLE}.vcf.gz.tbi"
	output:
		temp("08_baseRecal/{SAMPLE}_recal")
	conda:
		"snakemake/envs/default.yaml"
	resources:
		cpu = 1,
		ntasks = 2,
		mem_mb = 4000,
		hours = 2,
		mins = 0
	shell:
		"""
		gatk --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
            -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
            -Xloggc:gc_log.log -Xms4000m" \
            BaseRecalibrator \
            -R {input.refFa} \
            -I {input.bam} \
            --use-original-qualities \
            -O {output} \
            -known-sites {input.dbsnp}
		"""

rule applyRecal:
	input:
		bam = "05_addRG/bam/{SAMPLE}.bam",
		bamIndex = "05_addRG/bam/{SAMPLE}.bai",
		refFa = REFDIR + "Danio_rerio.GRCz11.dna.primary_assembly.fa",
		refIndex = REFDIR + "Danio_rerio.GRCz11.dna.primary_assembly.fa.fai",
		refDict = REFDIR + "Danio_rerio.GRCz11.dna.primary_assembly.dict",
		recal = "08_baseRecal/{SAMPLE}_recal"
	output:
		bam = temp("09_applyRecal/bam/{SAMPLE}.bam"),
		bamIndex = temp("09_applyRecal/bam/{SAMPLE}.bai")
	conda:
		"snakemake/envs/default.yaml"
	resources:
		cpu = 1,
		ntasks = 2,
		mem_mb = 4000,
		hours = 8,
		mins = 0
	shell:
		"""
		gatk \
            --java-options "-XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps \
            -XX:+PrintGCDetails -Xloggc:gc_log.log \
            -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms3000m" \
            ApplyBQSR \
            --add-output-sam-program-record \
            -R {input.refFa} \
            -I {input.bam} \
            --use-original-qualities \
            -O {output.bam} \
            --bqsr-recal-file {input.recal}
		"""

rule callVariants:
	input:
		bam = "09_applyRecal/bam/{SAMPLE}.bam",
		bamIndex = "09_applyRecal/bam/{SAMPLE}.bai",
		refFa = REFDIR + "Danio_rerio.GRCz11.dna.primary_assembly.fa",
		refIndex = REFDIR + "Danio_rerio.GRCz11.dna.primary_assembly.fa.fai",
		refDict = REFDIR + "Danio_rerio.GRCz11.dna.primary_assembly.dict"
	output:
		vcf = temp("10_callVariants/vcf/{SAMPLE}.vcf.gz"),
		vcfIndex = temp("10_callVariants/vcf/{SAMPLE}.vcf.gz.tbi")
	conda:
		"snakemake/envs/default.yaml"
	resources:
		cpu = 1,
		ntasks = 2,
		mem_mb = 8000,
		hours = 72,
		mins = 0
	shell:
		"""
		gatk --java-options "-Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
			HaplotypeCaller \
			-R {input.refFa} \
			-I {input.bam} \
			-O {output.vcf} \
			-dont-use-soft-clipped-bases \
			--standard-min-confidence-threshold-for-calling 20
		"""

rule filterVariants:
	input:
		vcf = "10_callVariants/vcf/{SAMPLE}.vcf.gz",
		vcfIndex = "10_callVariants/vcf/{SAMPLE}.vcf.gz.tbi",
		refFa = REFDIR + "Danio_rerio.GRCz11.dna.primary_assembly.fa",
		refIndex = REFDIR + "Danio_rerio.GRCz11.dna.primary_assembly.fa.fai",
		refDict = REFDIR + "Danio_rerio.GRCz11.dna.primary_assembly.dict"
	output:
		vcf = temp("11_filterVariants/vcf/{SAMPLE}.vcf.gz"),
		vcfIndex = temp("11_filterVariants/vcf/{SAMPLE}.vcf.gz.tbi")
	conda:
		"snakemake/envs/default.yaml"
	resources:
		cpu = 1,
		ntasks = 2,
		mem_mb = 1000,
		hours = 0,
		mins = 30
	shell:
		"""
		gatk \
		    VariantFiltration \
			--R {input.refFa} \
			--V {input.vcf} \
			--window 35 \
			--cluster 3 \
			--filter-name "FS" \
			--filter "FS > 30.0" \
			--filter-name "QD" \
			--filter "QD < 2.0" \
			-O {output.vcf}
		"""

rule selectVariants:
	input:
		vcf = "11_filterVariants/vcf/{SAMPLE}.vcf.gz",
		refFa = REFDIR + "Danio_rerio.GRCz11.dna.primary_assembly.fa"
	output:
		vcf = temp("12_selectVariants/vcf/{SAMPLE}.vcf.gz"),
		vcfIndex = temp("12_selectVariants/vcf/{SAMPLE}.vcf.gz.tbi")
	conda:
		"snakemake/envs/default.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 2000,
		hours = 0,
		mins = 30
	shell:
		"""
		gatk \
			SelectVariants \
			-R {input.refFa} \
			-V {input.vcf} \
			--select-type-to-include SNP \
			-O {output.vcf}
		"""

rule mergeSelected:
	input:
		vcf = expand("12_selectVariants/vcf/{SAMPLE}.vcf.gz", SAMPLE = SAMPLES),
		vcfIndex = expand("12_selectVariants/vcf/{SAMPLE}.vcf.gz.tbi", SAMPLE = SAMPLES)
	output:
		vcf = "12_selectVariants/mergedVcf/mergedVcf.vcf.gz",
		vcfIndex = "12_selectVariants/mergedVcf/mergedVcf.vcf.gz.tbi"
	conda:
		"snakemake/envs/default.yaml"
	resources:
		cpu = 4,
		ntasks = 2,
		mem_mb = 2000,
		hours = 2,
		mins = 0
	shell:
		"""
		bcftools merge --threads {resources.cpu} -o {output.vcf} -O z {input.vcf}
		bcftools index --threads {resources.cpu} -t {output.vcf}
		"""

rule countAlleles:
	input:
		bam = "09_applyRecal/bam/{SAMPLE}.bam",
		bamIndex = "09_applyRecal/bam/{SAMPLE}.bai",
		refFa = REFDIR + "Danio_rerio.GRCz11.dna.primary_assembly.fa",
		intervals = "13_countAlleles/intervals/snvs.intervals"
	output:
		counts = "13_countAlleles/counts/{SAMPLE}.alleleCounts.tsv"
	conda:
		"snakemake/envs/default.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 8000,
		hours = 2,
		mins = 0
	shell:
		"""
		gatk --java-options "-Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
			CollectAllelicCounts \
			-I {input.bam} \
			-R {input.refFa} \
			-L {input.intervals} \
			-O {output.counts} \
			--read-filter NotSecondaryAlignmentReadFilter \
			--read-filter GoodCigarReadFilter \
			--read-filter PassesVendorQualityCheckReadFilter \
			--read-filter MappingQualityAvailableReadFilter

		## Remove header
		sed -i '/^@/d' {output.counts}
		"""

rule aseReadCounts:
	input:
		bam = "09_applyRecal/bam/{SAMPLE}.bam",
		bamIndex = "09_applyRecal/bam/{SAMPLE}.bai",
		vcf = "12_selectVariants/vcf/{SAMPLE}.vcf.gz",
		refFa = REFDIR + "Danio_rerio.GRCz11.dna.primary_assembly.fa"
	output:
		tsv = "14_aseReadCounts/{SAMPLE}.tsv"
	conda:
		"snakemake/envs/default.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 4000,
		hours = 4,
		mins = 0
	shell:
		"""
		gatk \
			ASEReadCounter \
			-I {input.bam} \
			-V {input.vcf} \
			-R {input.refFa} \
			-O {output.tsv}
		"""