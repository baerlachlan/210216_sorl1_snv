rule recalBases:
	input:
		bam = "05_addRG/bam/{SAMPLE}.bam",
		bamIndex = "05_addRG/bam/{SAMPLE}.bai",
		refFa = "refs/Danio_rerio.GRCz11.dna.primary_assembly.fa",
		refIndex = "refs/Danio_rerio.GRCz11.dna.primary_assembly.fa.fai",
		refDict = "refs/Danio_rerio.GRCz11.dna.primary_assembly.dict",
		dbsnp = "06_knownSnvs/filtered/{SAMPLE}.vcf.gz",
		dbsnpIndex = "06_knownSnvs/filtered/{SAMPLE}.vcf.gz.tbi"
	output:
		temp("07_recalBases/recal/{SAMPLE}_recal")
	conda:
		"../envs/ase.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 4000,
		time = "00-01:00:00"
	shell:
		"""
		gatk --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
            -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
            -Xloggc:gc_log.log -Xms4000m" \
            BaseRecalibrator \
            -R {input.refFa} \
            -I {input.bam} \
            --use-original-qualities \
            -O {output} \
            -known-sites {input.dbsnp}
		"""

rule applyRecal:
	input:
		bam = "05_addRG/bam/{SAMPLE}.bam",
		bamIndex = "05_addRG/bam/{SAMPLE}.bai",
		refFa = "refs/Danio_rerio.GRCz11.dna.primary_assembly.fa",
		refIndex = "refs/Danio_rerio.GRCz11.dna.primary_assembly.fa.fai",
		refDict = "refs/Danio_rerio.GRCz11.dna.primary_assembly.dict",
		recal = "07_recalBases/recal/{SAMPLE}_recal"
	output:
		bam = "07_recalBases/bam/{SAMPLE}.bam",
		bamIndex = "07_recalBases/bam/{SAMPLE}.bai"
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
            --java-options "-XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps \
            -XX:+PrintGCDetails -Xloggc:gc_log.log \
            -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms3000m" \
            ApplyBQSR \
            --add-output-sam-program-record \
            -R {input.refFa} \
            -I {input.bam} \
            --use-original-qualities \
            -O {output.bam} \
            --bqsr-recal-file {input.recal}
		"""

rule fastqc_recalBases:
	input:
		"07_recalBases/bam/{SAMPLE}.bam"
	output:
		"07_recalBases/FastQC/{SAMPLE}_fastqc.zip",
		"07_recalBases/FastQC/{SAMPLE}_fastqc.html"
	params:
		outDir = "07_recalBases/FastQC/"
	conda:
		"../envs/ase.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 2000,
		time = "00-01:00:00"
	shell:
		"fastqc -t {resources.cpu} -o {params.outDir} --noextract {input}"