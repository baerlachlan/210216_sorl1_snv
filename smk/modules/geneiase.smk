rule geneiase:
	input:
		"12_geneiase/counts/{SAMPLE}.static.tsv"
	output:
		"12_geneiase/ase/{SAMPLE}.static.pval.tsv"
	conda:
		"../envs/geneiase.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 4000,
		time = "01-00:00:00"
	shell:
		"""
		../packages/geneiase-1.0.1/bin/geneiase \
			-t static \
			-i {input} \
			-o {output}
		"""