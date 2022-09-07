rule splitNCigar:
    input:
        bam = rules.markDuplicates.output.bam,
        bamIndex = rules.markDuplicates.output.bamIndex,
        refFa = rules.refs_downloadFa.output,
        refIndex = rules.refs_refIndex.output,
        refDict = rules.refs_refDict.output,
    output:
        bam = temp(os.path.join("results", splitNCigar_dir, "bam", "{SAMPLE}.bam")),
        bamIndex = temp(os.path.join("results", splitNCigar_dir, "bam", "{SAMPLE}.bai")),
        samstats = os.path.join("results", splitNCigar_dir, "samstats", "{SAMPLE}.tsv"),
    conda:
        "../envs/gatk.yml"
    resources:
        cpu = 8,
        ntasks = 1,
        mem_mb = 64000,
        time = "00-08:00:00",
    shell:
        """
        gatk --java-options "-Xms4g -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
            SplitNCigarReads \
            -R {input.refFa} \
            -I {input.bam} \
            -O {output.bam}

        samtools stats -d {output.bam} > {output.samstats}
        """