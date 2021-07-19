rule markDuplicates:
	input:
		bam = "02_align/bam/{SAMPLE}Aligned.sortedByCoord.out.bam"
	output:
		bam = temp("03_markDuplicates/bam/{SAMPLE}.bam"),
		bamIndex = temp("03_markDuplicates/bam/{SAMPLE}.bai"),
		metrics = "03_markDuplicates/log/{SAMPLE}.tsv"
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
		"""

rule fastqc_markDuplicates:
	input:
		"03_markDuplicates/bam/{SAMPLE}.bam"
	output:
		"03_markDuplicates/FastQC/{SAMPLE}_fastqc.zip",
		"03_markDuplicates/FastQC/{SAMPLE}_fastqc.html"
	params:
		outDir = "03_markDuplicates/FastQC/"
	conda:
		"../envs/ase.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 2000,
		time = "00-02:00:00"
	shell:
		"fastqc -t {resources.cpu} -o {params.outDir} --noextract {input}"