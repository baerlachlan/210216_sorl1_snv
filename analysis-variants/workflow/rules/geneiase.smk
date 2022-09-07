rule geneiase:
	input:
		tsv = "13_geneiase/1_counts/{SAMPLE}.static.tsv"
	output:
		tsv = "13_geneiase/2_ase/{SAMPLE}.static.pval.tsv"
	conda:
		"../envs/geneiase.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 4000,
		time = "01-00:00:00"
	shell:
		"""
		../../packages/geneiase-1.0.1/bin/geneiase \
			-t static \
			-m 1 \
			-i {input.tsv} \
			-o {output.tsv}
		"""