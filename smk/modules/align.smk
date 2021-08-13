rule align:
	input:
		R1 = "01_trim/fastq/{SAMPLE}.fastq.gz",
		starIndex = "refs/star/"
	output:
		bamRenamed = temp("02_align/bam/{SAMPLE}.bam"),
		bamIndex = temp("02_align/bam/{SAMPLE}.bam.bai"),
		STARgenome = temp(directory("02_align/bam/{SAMPLE}_STARgenome")),
		STARpass1 = temp(directory("02_align/bam/{SAMPLE}_STARpass1"))
	params:
		overhang = READ_LEN-1,
		bname = "02_align/bam/{SAMPLE}",
		bamBeforeRename = "02_align/bam/{SAMPLE}Aligned.sortedByCoord.out.bam"
	conda:
		"../envs/ase.yaml"
	resources:
		cpu = 16,
		ntasks = 1,
		mem_mb = 32000,
		time = "00-05:00:00"
	shell:
		"""
		STAR \
			--genomeDir {input.starIndex}\
			--runThreadN {resources.cpu} \
			--readFilesIn {input.R1} \
			--readFilesCommand "gunzip -c" \
			--sjdbOverhang {params.overhang} \
			--outSAMtype BAM SortedByCoordinate \
			--twopassMode Basic \
			--outFileNamePrefix {params.bname}

		mv {params.bamBeforeRename} {output.bamRenamed}

		mkdir -p 02_align/log
		mv {params.bname}*out 02_align/log
		mv {params.bname}*tab 02_align/log

		samtools index {output.bamRenamed}
		"""

rule fastqc_align:
	input:
		"02_align/bam/{SAMPLE}.bam"
	output:
		"02_align/FastQC/{SAMPLE}_fastqc.zip",
		"02_align/FastQC/{SAMPLE}_fastqc.html"
	params:
		outDir = "02_align/FastQC/"
	conda:
		"../envs/ase.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 2000,
		time = "00-01:00:00"
	shell:
		"fastqc -t {resources.cpu} -o {params.outDir} --noextract {input}"