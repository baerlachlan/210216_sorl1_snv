rule mergeSamples:
    input:
        unpack(mergeSamples_inputs),
    output:
        bam = temp(os.path.join("results", mergeSamples_dir, "bam", "{SAMPLE}.bam")),
        bamIndex = temp(os.path.join("results", mergeSamples_dir, "bam", "{SAMPLE}.bam.bai")),
    conda:
        "../envs/gatk.yml"
    resources:
        cpu = 1,
        ntasks = 1,
        mem_mb = 8000,
        time = "00-05:00:00",
    shell:
        """
        samtools merge {output.bam} {input.bam}
        samtools index {output.bam}
        """