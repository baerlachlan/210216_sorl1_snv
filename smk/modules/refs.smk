rule unzip_refFa:
	input:
		REFDIR + "Danio_rerio.GRCz11.dna.primary_assembly.fa.gz"
	output:
		temp("refs/Danio_rerio.GRCz11.dna.primary_assembly.fa")
	conda:
		"../envs/ase.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 4000,
		time = "00-00:30:00"
	shell:
		"gunzip -c {input} > {output}"

rule star_index:
	input:
		refFa = "refs/Danio_rerio.GRCz11.dna.primary_assembly.fa",
		gtf = REFDIR + "Danio_rerio.GRCz11.101.chr.gtf.gz"
	output:
		temp(directory("refs/star/"))
	params:
		overhang = READ_LEN-1
	conda:
		"../envs/ase.yaml"
	resources:
		cpu = 16,
		ntasks = 1,
		mem_mb = 32000,
		time = "00-00:30:00"
	shell:
		"""
		zcat {input.gtf} > temp.gtf

		STAR \
			--runThreadN {resources.cpu} \
			--runMode genomeGenerate \
			--genomeDir {output} \
			--genomeFastaFiles {input.refFa} \
			--sjdbGTFfile temp.gtf \
			--sjdbOverhang {params.overhang}

		rm temp.gtf
		"""

## Reference dictionary and index needs to be created as described in:
## https://gatk.broadinstitute.org/hc/en-us/articles/360035531652-FASTA-Reference-genome-format
rule ref_dict:
	input:
		"refs/Danio_rerio.GRCz11.dna.primary_assembly.fa"
	output:
		temp("refs/Danio_rerio.GRCz11.dna.primary_assembly.dict")
	conda:
		"../envs/ase.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 4000,
		time = "00-00:10:00"
	shell:
		"gatk CreateSequenceDictionary -R {input}"

rule ref_index:
	input:
		"refs/Danio_rerio.GRCz11.dna.primary_assembly.fa"
	output:
		temp("refs/Danio_rerio.GRCz11.dna.primary_assembly.fa.fai")
	conda:
		"../envs/ase.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 200,
		time = "00-00:10:00"
	shell:
		"samtools faidx {input}"