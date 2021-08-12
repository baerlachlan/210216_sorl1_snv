## Define a known set of variants based on the data, opposed to using a pre-existing database like dbSNP
## First we detect all potenital SNVs
rule callDbsnp:
	input:
		bam = "05_addRG/bam/{SAMPLE}.bam",
		bamIndex = "05_addRG/bam/{SAMPLE}.bai",
		refFa = "refs/Danio_rerio.GRCz11.dna.primary_assembly.fa",
		refIndex = "refs/Danio_rerio.GRCz11.dna.primary_assembly.fa.fai",
		refDict = "refs/Danio_rerio.GRCz11.dna.primary_assembly.dict"
	output:
		vcf = temp("06_dbsnp/1_called/{SAMPLE}.vcf.gz"),
		vcfIndex = temp("06_dbsnp/1_called/{SAMPLE}.vcf.gz.tbi")
	conda:
		"../envs/ase.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 8000,
		time = "00-12:00:00"
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

rule extractDbsnp:
	input:
		vcf = "06_dbsnp/1_called/{SAMPLE}.vcf.gz",
		vcfIndex = "06_dbsnp/1_called/{SAMPLE}.vcf.gz.tbi",
		refFa = "refs/Danio_rerio.GRCz11.dna.primary_assembly.fa",
		refIndex = "refs/Danio_rerio.GRCz11.dna.primary_assembly.fa.fai",
		refDict = "refs/Danio_rerio.GRCz11.dna.primary_assembly.dict"
	output:
		snvs = temp("06_dbsnp/2_extracted/{SAMPLE}_snvs.vcf.gz"),
		indels = temp("06_dbsnp/2_extracted/{SAMPLE}_indels.vcf.gz")
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
			-O {output.snvs}

		gatk \
			SelectVariants \
			-R {input.refFa} \
			-V {input.vcf} \
			--select-type-to-include INDEL \
			-O {output.indels}
		"""

## Here we filter the variants based on GATK's recommendations: https://gencore.bio.nyu.edu/variant-calling-pipeline-gatk4/
## These SNVs will be used for measuring base recalibration
rule filterDbsnp:
	input:
		snvs = "06_dbsnp/2_extracted/{SAMPLE}_snvs.vcf.gz",
		indels = "06_dbsnp/2_extracted/{SAMPLE}_indels.vcf.gz",
		refFa = "refs/Danio_rerio.GRCz11.dna.primary_assembly.fa",
		refIndex = "refs/Danio_rerio.GRCz11.dna.primary_assembly.fa.fai",
		refDict = "refs/Danio_rerio.GRCz11.dna.primary_assembly.dict"
	output:
		snvs = temp("06_dbsnp/3_filtered/{SAMPLE}_snvs.vcf.gz"),
		indels = temp("06_dbsnp/3_filtered/{SAMPLE}_indels.vcf.gz")
	conda:
		"../envs/ase.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 1000,
		time = "00-01:00:00"
	shell:
		"""
		gatk \
		    VariantFiltration \
			--R {input.refFa} \
			--V {input.snvs} \
			--window 35 \
			--cluster 3 \
			--filter-name "FS" --filter "FS > 60.0" \
			--filter-name "QD" --filter "QD < 2.0" \
			--filter-name "MQ" --filter "MQ < 40.0" \
			--filter-name "SOR" --filter "SOR > 4.0" \
			--filter-name "MQRankSum" --filter "MQRankSum < -12.5" \
			--filter-name "ReadPosRankSum" --filter "ReadPosRankSum < -8.0" \
			-O {output.snvs}

		gatk \
		    VariantFiltration \
			--R {input.refFa} \
			--V {input.indels} \
			--window 35 \
			--cluster 3 \
			--filter-name "FS" --filter "FS > 200.0" \
			--filter-name "QD" --filter "QD < 2.0" \
			--filter-name "SOR" --filter "SOR > 10.0" \
			-O {output.indels}
		"""

rule selectDbsnp:
	input:
		snvs = "06_dbsnp/3_filtered/{SAMPLE}_snvs.vcf.gz",
		indels = "06_dbsnp/3_filtered/{SAMPLE}_indels.vcf.gz",
		refFa = "refs/Danio_rerio.GRCz11.dna.primary_assembly.fa",
		refIndex = "refs/Danio_rerio.GRCz11.dna.primary_assembly.fa.fai",
		refDict = "refs/Danio_rerio.GRCz11.dna.primary_assembly.dict"
	output:
		snvs = temp("06_dbsnp/4_selected/{SAMPLE}_snvs.vcf.gz"),
		indels = temp("06_dbsnp/4_selected/{SAMPLE}_indels.vcf.gz")
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
			-V {input.snvs} \
			-O {output.snvs}

		gatk \
			SelectVariants \
			--exclude-filtered \
			-R {input.refFa} \
			-V {input.indels} \
			-O {output.indels}
		"""