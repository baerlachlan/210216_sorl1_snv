rule rmdup:
	input:
		bam = "02_align/bam/{SAMPLE}Aligned.sortedByCoord.out.bam"
	output:
		bam = temp("03_rmdup/bam/{SAMPLE}Aligned.sortedByCoord.out.bam"),
		bamIndex = temp("03_rmdup/bam/{SAMPLE}Aligned.sortedByCoord.out.bam.bai")
	conda:
		"../envs/wasp.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 32000,
		time = "00-04:00:00"
	shell:
		"""
		python ../packages/WASP/mapping/rmdup.py {input.bam} {output.bam}

		samtools index {output.bam}
		"""