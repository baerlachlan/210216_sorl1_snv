## When MarkDuplicates (Picard) is run on coordinate sorted BAM files, unmapped mates of mapped records and supplementary/secondary alignments are excluded from the duplication test
## For variant analysis with GATK this is not a problem because HaplotypeCaller filters unmapped reads and secondary alignments before analysing
rule markDuplicates:
	input:
		bam = "02_align/bam/{SAMPLE}.bam"
	output:
		bam = temp("03_markDuplicates/bam/{SAMPLE}.bam"),
		bamIndex = temp("03_markDuplicates/bam/{SAMPLE}.bai"),
		metrics = "03_markDuplicates/metrics/{SAMPLE}.tsv",
		samstats = "03_markDuplicates/samstats/{SAMPLE}.tsv"
	conda:
		"../envs/ase.yaml"
	resources:
		cpu = 1,
		ntasks = 2,
		mem_mb = 4000,
		time = "00-01:00:00"
	shell:
		"""
		gatk \
			MarkDuplicates \
			--INPUT {input} \
			--OUTPUT {output.bam}  \
			--DUPLICATE_SCORING_STRATEGY RANDOM \
			--CREATE_INDEX true \
			--METRICS_FILE {output.metrics}

		samtools stats {output.bam} > {output.samstats}
		"""