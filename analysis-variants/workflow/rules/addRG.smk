rule addRG:
    input:
        bam = rules.align.output.bam,
        bamIndex = rules.align.output.bamIndex,
    output:
        bam = temp(os.path.join("results", addRG_dir, "bam", "{SAMPLE}{MERGETAG, .*}.bam")),
        bamIndex = temp(os.path.join("results", addRG_dir, "bam", "{SAMPLE}{MERGETAG, .*}.bam.bai")),
        samstats = os.path.join("results", addRG_dir, "samstats", "{SAMPLE}{MERGETAG, .*}.tsv"),
    # params:
    #     RGID = lambda wildcard: analysis.RGID[wildcard.SAMPLE],
    #     RGSM = lambda wildcard: analysis.RGSM[wildcard.SAMPLE],
    conda:
        "../envs/gatk.yml"
    resources:
        cpu = 1,
        ntasks = 2,
        mem_mb = 4000,
        time = "00-02:00:00",
    shell:
        """
        gatk \
            AddOrReplaceReadGroups \
            -I {input.bam} \
            -O {output.bam} \
            -SORT_ORDER coordinate \
            -RGID {wildcards.SAMPLE}{wildcards.MERGETAG} \
            -RGPU null \
            -RGSM {wildcards.SAMPLE} \
            -RGPL ILLUMINA \
            -RGLB null \
            -CREATE_INDEX False

        samtools index {output.bam}
        samtools stats -d {output.bam} > {output.samstats}
        """