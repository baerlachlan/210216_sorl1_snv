rule splitNCigar:
	input:
		bam = "03_markDuplicates/bam/{SAMPLE}.bam",
		bamIndex = "03_markDuplicates/bam/{SAMPLE}.bai",
		refFa = "refs/Danio_rerio.GRCz11.dna.primary_assembly.fa",
		refIndex = "refs/Danio_rerio.GRCz11.dna.primary_assembly.fa.fai",
		refDict = "refs/Danio_rerio.GRCz11.dna.primary_assembly.dict"
	output:
		bam = temp("04_splitNCigar/bam/{SAMPLE}.bam"),
		bamIndex = temp("04_splitNCigar/bam/{SAMPLE}.bai"),
		samstats = "04_splitNCigar/samstats/{SAMPLE}.tsv"
	conda:
		"../envs/ase.yaml"
	resources:
		cpu = 4,
		ntasks = 1,
		mem_mb = 16000,
		time = "00-04:00:00"
	shell:
		"""
		gatk \
            SplitNCigarReads \
            -R {input.refFa} \
            -I {input.bam} \
            -O {output.bam}

		samtools stats -d {output.bam} > {output.samstats}
		"""