rule groupUmis:
	input:
		bam = rules.addRG.output.bam,
		bamIndex = rules.addRG.output.bamIndex,
	output:
		bam = temp(os.path.join("results", groupUmis_dir, "bam", "{SAMPLE}{MERGETAG, .*}.bam")),
		bamIndex = temp(os.path.join("results", groupUmis_dir, "bam", "{SAMPLE}{MERGETAG, .*}.bam.bai")),
	conda:
		"../envs/gatk.yml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 8000,
		time = "00-05:00:00",
	shell:
		"""
		umi_tools group \
			-I {input.bam} \
			-S {output.bam} \
			--temp-dir=. \
			--output-bam \
			--method=unique \
			--extract-umi-method=read_id \
			--umi-separator=":" \
			--paired

		samtools index {output.bam}
		"""