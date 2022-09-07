rule bqsr_firstPass:
    input:
        unpack(knownVariants_files),
        bam = rules.splitNCigar.output.bam,
        bamIndex = rules.splitNCigar.output.bamIndex,
        refFa = rules.refs_downloadFa.output,
        refIndex = rules.refs_refIndex.output,
        refDict = rules.refs_refDict.output,
    output:
        temp(os.path.join("results", bqsr_dir, "recal", "{SAMPLE}.firstPass.table")),
    conda:
        "../envs/gatk.yml"
    resources:
        cpu = 1,
        ntasks = 1,
        mem_mb = 8000,
        time = "00-04:00:00",
    shell:
        """
        gatk --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
            -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
            -Xloggc:gc_log.log -Xms4000m" \
            BaseRecalibrator \
            -R {input.refFa} \
            -I {input.bam} \
            -O {output} \
            --known-sites {input.knownVariants}
        """

rule bqsr_apply:
    input:
        bam = rules.splitNCigar.output.bam,
        bamIndex = rules.splitNCigar.output.bamIndex,
        refFa = rules.refs_downloadFa.output,
        refIndex = rules.refs_refIndex.output,
        refDict = rules.refs_refDict.output,
        recal = rules.bqsr_firstPass.output,
    output:
        bam = os.path.join("results", bqsr_dir, "bam", "{SAMPLE}.bam"),
        bamIndex = os.path.join("results", bqsr_dir, "bam", "{SAMPLE}.bai"),
        metrics = os.path.join("results", bqsr_dir, "metrics", "{SAMPLE}.tsv"),
    conda:
        "../envs/gatk.yml"
    resources:
        cpu = 1,
        ntasks = 1,
        mem_mb = 8000,
        time = "00-04:00:00",
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
            -O {output.bam} \
            --bqsr-recal-file {input.recal}

        samtools stats -d {output.bam} > {output.metrics}
        """

rule bqsr_secondPass:
    input:
        unpack(knownVariants_files),
        bam = rules.bqsr_apply.output.bam,
        bamIndex = rules.bqsr_apply.output.bamIndex,
        refFa = rules.refs_downloadFa.output,
        refIndex = rules.refs_refIndex.output,
        refDict = rules.refs_refDict.output,
    output:
        temp(os.path.join("results", bqsr_dir, "recal", "{SAMPLE}.secondPass.table")),
    conda:
        "../envs/gatk.yml"
    resources:
        cpu = 1,
        ntasks = 1,
        mem_mb = 8000,
        time = "00-04:00:00",
    shell:
        """
        gatk --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
            -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
            -Xloggc:gc_log.log -Xms4000m" \
            BaseRecalibrator \
            -R {input.refFa} \
            -I {input.bam} \
            -O {output} \
            --known-sites {input.knownVariants}
        """

rule bqsr_analyzeCovariates:
    input:
        firstPass = rules.bqsr_firstPass.output,
        secondPass = rules.bqsr_secondPass.output,
    output:
        os.path.join("results", bqsr_dir, "recal", "{SAMPLE}.analyzeCovariates.csv"),
    conda:
        "../envs/gatk.yml"
    resources:
        cpu = 1,
        ntasks = 1,
        mem_mb = 4000,
        time = "00-01:00:00",
    shell:
        """
        gatk AnalyzeCovariates \
            -before {input.firstPass} \
             -after {input.secondPass} \
             -csv {output}
        """