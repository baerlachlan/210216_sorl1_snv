rule aseRC:
    input:
        bam = rules.wasp_merge.output.keep_sorted,
        bamIndex = rules.wasp_merge.output.keep_sortedIndex,
        vcf = rules.variants_select.output.vcf,
        refFa = rules.refs_downloadFa.output,
        intervals = os.path.join("results", aseRC_dir, "intervals", "{SAMPLE}.intervals")
    output:
        tsv = os.path.join("results", aseRC_dir, "wasp", "{SAMPLE}.tsv"),
    conda:
        "../envs/gatk.yml"
    resources:
        cpu = 1,
        ntasks = 1,
        mem_mb = 4000,
        time = "00-04:00:00",
    shell:
        """
        gatk \
            ASEReadCounter \
            -I {input.bam} \
            -V {input.vcf} \
            -R {input.refFa} \
            -L {input.intervals} \
            -O {output.tsv} \
            --min-mapping-quality 10 \
            --min-base-quality 20
        """

rule aseRC_nowasp:
    input:
        bam = rules.bqsr_apply.output.bam,
        bamIndex = rules.bqsr_apply.output.bamIndex,
        vcf = rules.variants_select.output.vcf,
        refFa = rules.refs_downloadFa.output,
        intervals = os.path.join("results", aseRC_dir, "intervals", "{SAMPLE}.intervals")
    output:
        tsv = os.path.join("results", aseRC_dir, "no_wasp", "{SAMPLE}.tsv"),
    conda:
        "../envs/gatk.yml"
    resources:
        cpu = 1,
        ntasks = 1,
        mem_mb = 4000,
        time = "00-04:00:00",
    shell:
        """
        gatk \
            ASEReadCounter \
            -I {input.bam} \
            -V {input.vcf} \
            -R {input.refFa} \
            -L {input.intervals} \
            -O {output.tsv} \
            --min-mapping-quality 10 \
            --min-base-quality 20
        """