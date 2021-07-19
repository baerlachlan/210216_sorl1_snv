rule trim:
	input:
		R1 = "00_rawData/fastq/{SAMPLE}.fastq.gz"
	output:
		R1 = temp("01_trim/fastq/{SAMPLE}.fastq.gz"),
		html = "01_trim/log/{SAMPLE}.html"
	conda:
		"../envs/ase.yaml"
	resources:
		cpu = 1,
		ntasks = 2,
		mem_mb = 2000,
		time = "00-02:00:00"
	shell:
		"""
		fastp \
			-i {input.R1}  \
        	-o {output.R1} \
			--qualified_quality_phred 20 \
			--length_required 35 \
			--trim_poly_g \
			--thread 1 \
			--html {output.html} \
			--json /dev/null \
		"""

rule fastqc_trim:
	input:
		"01_trim/fastq/{SAMPLE}.fastq.gz"
	output:
		"01_trim/FastQC/{SAMPLE}_fastqc.zip",
		"01_trim/FastQC/{SAMPLE}_fastqc.html"
	params:
		outDir = "01_trim/FastQC/"
	conda:
		"../envs/ase.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 2000,
		time = "00-01:00:00"
	shell:
		"fastqc -t {resources.cpu} -o {params.outDir} --noextract {input}"