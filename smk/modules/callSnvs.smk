rule callVariants:
	input:
		bam = "07_recalBases/bam/{SAMPLE}.bam",
		bamIndex = "07_recalBases/bam/{SAMPLE}.bai",
		refFa = "refs/Danio_rerio.GRCz11.dna.primary_assembly.fa",
		refIndex = "refs/Danio_rerio.GRCz11.dna.primary_assembly.fa.fai",
		refDict = "refs/Danio_rerio.GRCz11.dna.primary_assembly.dict",
		intervals = "refs/exons.intervals",
		snvs = "06_dbsnp/4_selected/{SAMPLE}_snvs.vcf.gz",
		indels = "06_dbsnp/4_selected/{SAMPLE}_indels.vcf.gz"
	output:
		vcf = temp("08_callSnvs/1_called/{SAMPLE}.vcf.gz"),
		vcfIndex = temp("08_callSnvs/1_called/{SAMPLE}.vcf.gz.tbi"),
		detailMetrics = "08_callSnvs/1_called/log/{SAMPLE}.variant_calling_detail_metrics",
		summaryMetrics = "08_callSnvs/1_called/log/{SAMPLE}.variant_calling_summary_metrics"
	params:
		metricsBname = "08_callSnvs/1_called/log/{SAMPLE}"
	conda:
		"../envs/ase.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 8000,
		time = "00-04:00:00"
	shell:
		"""
		gatk --java-options "-Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
			HaplotypeCaller \
			-R {input.refFa} \
			-I {input.bam} \
			-L {input.intervals} \
			-O {output.vcf} \
			-dont-use-soft-clipped-bases \
			--standard-min-confidence-threshold-for-calling 20

		gatk \
			CollectVariantCallingMetrics \
			--DBSNP {input.snvs} \
			--INPUT {output.vcf} \
			--OUTPUT {params.metricsBname}
		"""

## We are only interested in Single Nucleotide Variants for ASE analysis
rule extractSnvs:
	input:
		vcf = "08_callSnvs/1_called/{SAMPLE}.vcf.gz",
		vcfIndex = "08_callSnvs/1_called/{SAMPLE}.vcf.gz.tbi",
		refFa = "refs/Danio_rerio.GRCz11.dna.primary_assembly.fa",
		refIndex = "refs/Danio_rerio.GRCz11.dna.primary_assembly.fa.fai",
		refDict = "refs/Danio_rerio.GRCz11.dna.primary_assembly.dict",
		snvs = "06_dbsnp/4_selected/{SAMPLE}_snvs.vcf.gz"
	output:
		vcf = temp("08_callSnvs/2_extracted/{SAMPLE}.vcf.gz"),
		vcfIndex = temp("08_callSnvs/2_extracted/{SAMPLE}.vcf.gz.tbi"),
		detailMetrics = "08_callSnvs/2_extracted/log/{SAMPLE}.variant_calling_detail_metrics",
		summaryMetrics = "08_callSnvs/2_extracted/log/{SAMPLE}.variant_calling_summary_metrics"
	params:
		metricsBname = "08_callSnvs/2_extracted/log/{SAMPLE}"
	conda:
		"../envs/ase.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 8000,
		time = "00-12:00:00"
	shell:
		"""
		gatk \
			SelectVariants \
			-R {input.refFa} \
			-V {input.vcf} \
			--select-type-to-include SNP \
			--restrict-alleles-to BIALLELIC \
			-O {output.vcf}

		gatk \
			CollectVariantCallingMetrics \
			--DBSNP {input.snvs} \
			--INPUT {output.vcf} \
			--OUTPUT {params.metricsBname}
		"""

rule filterSnvs:
	input:
		vcf = "08_callSnvs/2_extracted/{SAMPLE}.vcf.gz",
		vcfIndex = "08_callSnvs/2_extracted/{SAMPLE}.vcf.gz.tbi",
		refFa = "refs/Danio_rerio.GRCz11.dna.primary_assembly.fa",
		refIndex = "refs/Danio_rerio.GRCz11.dna.primary_assembly.fa.fai",
		refDict = "refs/Danio_rerio.GRCz11.dna.primary_assembly.dict",
		snvs = "06_dbsnp/4_selected/{SAMPLE}_snvs.vcf.gz"
	output:
		vcf = temp("08_callSnvs/3_filtered/{SAMPLE}.vcf.gz"),
		vcfIndex = temp("08_callSnvs/3_filtered/{SAMPLE}.vcf.gz.tbi"),
		detailMetrics = "08_callSnvs/3_filtered/log/{SAMPLE}.variant_calling_detail_metrics",
		summaryMetrics = "08_callSnvs/3_filtered/log/{SAMPLE}.variant_calling_summary_metrics"
	params:
		metricsBname = "08_callSnvs/3_filtered/log/{SAMPLE}"
	conda:
		"../envs/ase.yaml"
	resources:
		cpu = 1,
		ntasks = 2,
		mem_mb = 2000,
		time = "00-00:30:00"
	shell:
		"""
		gatk \
		    VariantFiltration \
			--R {input.refFa} \
			--V {input.vcf} \
			--window 35 \
			--cluster 3 \
			--filter-name "FS" --filter "FS > 60.0" \
			--filter-name "QD" --filter "QD < 2.0" \
			--filter-name "MQ" --filter "MQ < 40.0" \
			--filter-name "SOR" --filter "SOR > 4.0" \
			--filter-name "MQRankSum" --filter "MQRankSum < -12.5" \
			--filter-name "ReadPosRankSum" --filter "ReadPosRankSum < -8.0" \
			-O {output.vcf}

		gatk \
			CollectVariantCallingMetrics \
			--DBSNP {input.snvs} \
			--INPUT {output.vcf} \
			--OUTPUT {params.metricsBname}
		"""

rule selectSnvs:
	input:
		vcf = "08_callSnvs/3_filtered/{SAMPLE}.vcf.gz",
		vcfIndex = "08_callSnvs/3_filtered/{SAMPLE}.vcf.gz.tbi",
		refFa = "refs/Danio_rerio.GRCz11.dna.primary_assembly.fa",
		refIndex = "refs/Danio_rerio.GRCz11.dna.primary_assembly.fa.fai",
		refDict = "refs/Danio_rerio.GRCz11.dna.primary_assembly.dict",
		snvs = "06_dbsnp/4_selected/{SAMPLE}_snvs.vcf.gz"
	output:
		vcf = "08_callSnvs/4_selected/{SAMPLE}.vcf.gz",
		vcfIndex = "08_callSnvs/4_selected/{SAMPLE}.vcf.gz.tbi",
		detailMetrics = "08_callSnvs/4_selected/log/{SAMPLE}.variant_calling_detail_metrics",
		summaryMetrics = "08_callSnvs/4_selected/log/{SAMPLE}.variant_calling_summary_metrics"
	params:
		metricsBname = "08_callSnvs/4_selected/log/{SAMPLE}"
	conda:
		"../envs/ase.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 2000,
		time = "00-00:30:00"
	shell:
		"""
		gatk \
			SelectVariants \
			--exclude-filtered \
			-R {input.refFa} \
			-V {input.vcf} \
			-O {output.vcf}

		gatk \
			CollectVariantCallingMetrics \
			--DBSNP {input.snvs} \
			--INPUT {output.vcf} \
			--OUTPUT {params.metricsBname}
		"""