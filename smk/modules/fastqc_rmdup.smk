rule fastqc_dedup:
	input:
		"03_rmdup/bam/{SAMPLE}.bam"
	output:
		"03_rmdup/FastQC/{SAMPLE}_fastqc.zip",
		"03_rmdup/FastQC/{SAMPLE}_fastqc.html"
	params:
		outDir = "03_rmdup/FastQC/"
	conda:
		"../envs/ase.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 2000,
		time = "00-02:00:00"
	shell:
		"fastqc -t {resources.cpu} -o {params.outDir} --noextract {input}"