## gatk HaplotypeCaller requires read group (RG) tags
## An explanation of this is found at https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups
rule addRG:
	input:
		bam = "04_splitNCigar/bam/{SAMPLE}.bam",
		bamIndex = "04_splitNCigar/bam/{SAMPLE}.bai"
	output:
		bam = temp("05_addRG/bam/{SAMPLE}.bam"),
		bamIndex = temp("05_addRG/bam/{SAMPLE}.bai")
	conda:
		"../envs/ase.yaml"
	resources:
		cpu = 1,
		ntasks = 2,
		mem_mb = 4000,
		time = "00-02:00:00"
	shell:
		"""
		gatk \
			AddOrReplaceReadGroups \
    		-I {input.bam} \
   			-O {output.bam} \
    		-SORT_ORDER coordinate \
    		-RGID default \
    		-RGLB default \
    		-RGSM {wildcards.SAMPLE} \
			-RGPU default \
    		-RGPL ILLUMINA \
    		-CREATE_INDEX True
		"""

rule addRG_metrics:
	input:
		bam = "05_addRG/bam/{SAMPLE}.bam"
	output:
		metrics = "05_addRG/metrics/{SAMPLE}.tsv"
	conda:
		"../envs/ase.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 4000,
		time = "00-00:10:00"
	shell:
		"""
		samtools stats -d {input.bam} > {output.metrics}
		"""

# ## This rule provides metrics of base coverage, which is important to inspect for base recalibration purposes
# ## Base recalibration requires a coverage of at least 100 million bases to be effective
# rule baseCoverage:
# 	input:
# 		bam = "05_addRG/bam/{SAMPLE}.bam",
# 		bed = "08_callSnvs/intervals/exons.bed"
# 		# bam = "05_addRG/bam/{SAMPLE}.bam",
# 		# bamIndex = "05_addRG/bam/{SAMPLE}.bai"
# 	output:
# 		cov = "05_addRG/log/{SAMPLE}.baseCoverage"
# 	conda:
# 		"../envs/ase.yaml"
# 	resources:
# 		cpu = 1,
# 		ntasks = 1,
# 		mem_mb = 2000,
# 		time = "00-00:10:00"
# 	shell:
# 		"""
# 		samtools depth {input.bam} | \
# 		awk '{{a+=$3}}END{{print a}}' > \
# 		{output.cov}
# 		"""