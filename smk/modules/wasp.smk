rule findIntersecting:
	input:
		bam = "07_recalBases/bam/{SAMPLE}.bam",
		snpDir = "09_wasp/1_snvs/{SAMPLE}"
	output:
		fq = temp("09_wasp/2_findIntersecting/{SAMPLE}/{SAMPLE}.remap.fq.gz"),
		to_remap = temp("09_wasp/2_findIntersecting/{SAMPLE}/{SAMPLE}.to.remap.bam"),
		keep_intersect = temp("09_wasp/2_findIntersecting/{SAMPLE}/{SAMPLE}.keep.bam")
	params:
		outDir = temp(directory("09_wasp/2_findIntersecting/{SAMPLE}"))
	conda:
		"../envs/wasp.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 2000,
		time = "00-01:00:00"
	shell:
		"""
		python ../packages/WASP/mapping/find_intersecting_snps.py \
    	    --is_sorted \
			--output_dir {params.outDir} \
    	    --snp_dir {input.snpDir} \
    	    {input.bam}
		"""

rule remap:
	input:
		fq = "09_wasp/2_findIntersecting/{SAMPLE}/{SAMPLE}.remap.fq.gz",
		starIndex = "refs/star/"
	output:
		remapped_unsorted = temp("09_wasp/3_remap/{SAMPLE}Aligned.out.bam"),
		remapped_sorted = temp("09_wasp/3_remap/{SAMPLE}sorted.out.bam"),
		index = temp("09_wasp/3_remap/{SAMPLE}sorted.out.bam.bai"),
		STARgenome = temp(directory("09_wasp/3_remap/{SAMPLE}_STARgenome")),
		STARpass1 = temp(directory("09_wasp/3_remap/{SAMPLE}_STARpass1"))
	params:
		overhang = READ_LEN-1,
		bname = "09_wasp/3_remap/{SAMPLE}"
	conda:
		"../envs/wasp.yaml"
	resources:
		cpu = 16,
		ntasks = 1,
		mem_mb = 32000,
		time = "00-01:00:00"
	shell:
		"""
		STAR \
			--genomeDir {input.starIndex}\
			--runThreadN {resources.cpu} \
			--readFilesIn {input.fq} \
			--readFilesCommand "gunzip -c" \
			--sjdbOverhang {params.overhang} \
			--outSAMtype BAM Unsorted \
			--twopassMode Basic \
			--outFileNamePrefix {params.bname}

		mkdir -p 09_wasp/3_remap/log
		mv {params.bname}*out 09_wasp/3_remap/log
		mv {params.bname}*tab 09_wasp/3_remap/log

		samtools sort -o {output.remapped_sorted} {output.remapped_unsorted}
		samtools index {output.remapped_sorted}
		"""

rule filterRemapped:
	input:
		to_remap = "09_wasp/2_findIntersecting/{SAMPLE}/{SAMPLE}.to.remap.bam",
		remapped_unsorted = "09_wasp/3_remap/{SAMPLE}Aligned.out.bam",
		remapped_sorted = "09_wasp/3_remap/{SAMPLE}sorted.out.bam"
	output:
		keep_filter = temp("09_wasp/4_filterRemapped/{SAMPLE}.keep.bam")
	conda:
		"../envs/wasp.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 24000,
		time = "00-00:30:00"
	shell:
		"""
		python ../packages/WASP/mapping/filter_remapped_reads.py \
			{input.to_remap} \
			{input.remapped_sorted} \
			{output.keep_filter}
		"""

rule merge:
	input:
		keep_filter = "09_wasp/4_filterRemapped/{SAMPLE}.keep.bam",
		keep_intersect = "09_wasp/2_findIntersecting/{SAMPLE}/{SAMPLE}.keep.bam"
	output:
		keep_merged = temp("09_wasp/5_merge/{SAMPLE}.keep.merge.bam"),
		keep_sorted = temp("09_wasp/5_merge/{SAMPLE}.keep.merge.sort.bam"),
		keep_sortedIndex = temp("09_wasp/5_merge/{SAMPLE}.keep.merge.sort.bam.bai")
	conda:
		"../envs/wasp.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 4000,
		time = "00-00:30:00"
	shell:
		"""
		samtools merge {output.keep_merged} \
			{input.keep_filter} \
			{input.keep_intersect}
		samtools sort -o {output.keep_sorted} \
			{output.keep_merged}
		samtools index {output.keep_sorted}
		"""
