rule aseRC:
	input:
		bam = "09_wasp/5_merge/{SAMPLE}.keep.merge.sort.bam",
		bamIndex = "09_wasp/5_merge/{SAMPLE}.keep.merge.sort.bam.bai",
		vcf = "08_callSnvs/4_selected/{SAMPLE}.vcf.gz",
		refFa = "refs/Danio_rerio.GRCz11.dna.primary_assembly.fa",
		intervals = "08_callSnvs/intervals/exons.intervals"
	output:
		tsv = "10_aseReadCounter/wasp/{SAMPLE}.tsv"
	conda:
		"../envs/ase.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 4000,
		time = "00-01:00:00"
	shell:
		"""
		gatk \
			ASEReadCounter \
			-I {input.bam} \
			-V {input.vcf} \
			-R {input.refFa} \
			-L {input.intervals} \
			-O {output.tsv} \
			--min-mapping-quality 10 \
			--min-base-quality 20
		"""

rule aseRC_nowasp:
	input:
		bam = "07_recalBases/bam/{SAMPLE}.bam",
		bamIndex = "07_recalBases/bam/{SAMPLE}.bai",
		vcf = "08_callSnvs/4_selected/{SAMPLE}.vcf.gz",
		refFa = "refs/Danio_rerio.GRCz11.dna.primary_assembly.fa",
		intervals = "08_callSnvs/intervals/exons.intervals"
	output:
		tsv = "10_aseReadCounter/nowasp/{SAMPLE}.tsv"
	conda:
		"../envs/ase.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 4000,
		time = "00-01:00:00"
	shell:
		"""
		gatk \
			ASEReadCounter \
			-I {input.bam} \
			-V {input.vcf} \
			-R {input.refFa} \
			-L {input.intervals} \
			-O {output.tsv} \
			--min-mapping-quality 10 \
			--min-base-quality 20
		"""