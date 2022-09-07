rule variants_gvcf:
    input:
        unpack(knownVariants_files),
        bam = rules.bqsr_apply.output.bam,
        bamIndex = rules.bqsr_apply.output.bamIndex,
        refFa = rules.refs_downloadFa.output,
        refIndex = rules.refs_refIndex.output,
        refDict = rules.refs_refDict.output,
        intervals = rules.intervals.output,
    output:
        gvcf = temp(os.path.join("results", variants_dir, "1_gvcf", "{SAMPLE}.g.vcf.gz")),
        gvcfIndex = temp(os.path.join("results", variants_dir, "1_gvcf", "{SAMPLE}.g.vcf.gz.tbi")),
        detailMetrics = os.path.join("results", variants_dir, "1_gvcf", "log", "{SAMPLE}.variant_calling_detail_metrics"),
        summaryMetrics = os.path.join("results", variants_dir, "1_gvcf", "log", "{SAMPLE}.variant_calling_summary_metrics"),
    params:
        metricsBname = os.path.join("results", variants_dir, "1_gvcf", "log", "{SAMPLE}"),
    conda:
        "../envs/gatk.yml"
    resources:
        cpu = 1,
        ntasks = 1,
        mem_mb = 8000,
        time = "03-00:00:00",
    shell:
        """
        gatk --java-options "-Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
            HaplotypeCaller \
            -R {input.refFa} \
            -I {input.bam} \
            -L {input.intervals} \
            -O {output.gvcf} \
            -dont-use-soft-clipped-bases \
            --standard-min-confidence-threshold-for-calling 20 \
            --emit-ref-confidence GVCF

        gatk \
            CollectVariantCallingMetrics \
            --GVCF_INPUT true \
            --DBSNP {input.knownVariants} \
            --INPUT {output.gvcf} \
            --OUTPUT {params.metricsBname}
        """

rule variants_sampleMap:
    input:
        expand(rules.variants_gvcf.output.gvcf, SAMPLE=samples),
    output:
        os.path.join("results", variants_dir, "sample_map.tsv"),
    params:
        variants_dir = variants_dir,
    conda:
        "../envs/r.yml"
    resources:
        cpu = 1,
        ntasks = 1,
        mem_mb = 2000,
        time = "00-00:10:00",
    script:
        "../scripts/sample_map.R"

rule variants_genomicsDB:
    input:
        gvcf = expand(rules.variants_gvcf.output.gvcf, SAMPLE=samples),
        gvcfIndex = expand(rules.variants_gvcf.output.gvcfIndex, SAMPLE=samples),
        sample_map = rules.variants_sampleMap.output,
        intervals = rules.intervals.output,
    output:
        genDB_dir = temp(directory(os.path.join("results", variants_dir, "2_genomicsDB"))),
    conda:
        "../envs/gatk.yml"
    resources:
        cpu = 1,
        ntasks = 1,
        mem_mb = 8000,
        time = "01-00:00:00",
    shell:
        """
        gatk --java-options "-Xmx4g -Xms4g" \
            GenomicsDBImport \
                --genomicsdb-workspace-path {output.genDB_dir} \
                --intervals {input.intervals} \
                --sample-name-map {input.sample_map} \
                --tmp-dir . \
                --merge-input-intervals
        """

rule variants_genotype:
    input:
        unpack(knownVariants_files),
        genDB_dir = rules.variants_genomicsDB.output.genDB_dir,
        refFa = rules.refs_downloadFa.output,
        refIndex = rules.refs_refIndex.output,
        refDict = rules.refs_refDict.output,
    output:
        vcf = temp(os.path.join("results", variants_dir, "3_jointGenotype", "all_samples.vcf.gz")),
        vcfIndex = temp(os.path.join("results", variants_dir, "3_jointGenotype", "all_samples.vcf.gz.tbi")),
        detailMetrics = os.path.join("results", variants_dir, "3_jointGenotype", "log", "all_samples.variant_calling_detail_metrics"),
        summaryMetrics = os.path.join("results", variants_dir, "3_jointGenotype", "log", "all_samples.variant_calling_summary_metrics"),
    params:
        gendb = os.path.join("gendb://" + "results", variants_dir, "2_genomicsDB"),
        metricsBname = os.path.join("results", variants_dir, "3_jointGenotype", "log", "all_samples"),
    conda:
        "../envs/gatk.yml"
    resources:
        cpu = 1,
        ntasks = 2,
        mem_mb = 32000,
        time = "01-00:00:00",
    shell:
        """
        gatk --java-options "-Xmx16g" \
            GenotypeGVCFs \
                -R {input.refFa} \
                -V {params.gendb} \
                -O {output.vcf}

        gatk \
            CollectVariantCallingMetrics \
            --DBSNP {input.knownVariants} \
            --INPUT {output.vcf} \
            --OUTPUT {params.metricsBname}
        """

rule variants_extract:
    input:
        unpack(knownVariants_files),
        vcf = rules.variants_genotype.output.vcf,
        vcfIndex = rules.variants_genotype.output.vcfIndex,
        refFa = rules.refs_downloadFa.output,
        refIndex = rules.refs_refIndex.output,
        refDict = rules.refs_refDict.output,
    output:
        vcf = os.path.join("results", variants_dir, "4_extract", "all_samples.vcf.gz"),
        vcfIndex = os.path.join("results", variants_dir, "4_extract", "all_samples.vcf.gz.tbi"),
        detailMetrics = os.path.join("results", variants_dir, "4_extract", "log", "all_samples.variant_calling_detail_metrics"),
        summaryMetrics = os.path.join("results", variants_dir, "4_extract", "log", "all_samples.variant_calling_summary_metrics"),
    params:
        metricsBname = os.path.join("results", variants_dir, "4_extract", "log", "all_samples"),
    conda:
        "../envs/gatk.yml"
    resources:
        cpu = 1,
        ntasks = 1,
        mem_mb = 8000,
        time = "00-01:00:00",
    shell:
        """
        gatk \
            SelectVariants \
            -R {input.refFa} \
            -V {input.vcf} \
            --select-type-to-include SNP \
            -O {output.vcf}

        gatk \
            CollectVariantCallingMetrics \
            --DBSNP {input.knownVariants} \
            --INPUT {output.vcf} \
            --OUTPUT {params.metricsBname}
        """

rule variants_filter:
    input:
        unpack(knownVariants_files),
        vcf = rules.variants_extract.output.vcf,
        vcfIndex = rules.variants_extract.output.vcfIndex,
        refFa = rules.refs_downloadFa.output,
        refIndex = rules.refs_refIndex.output,
        refDict = rules.refs_refDict.output,
    output:
        vcf = temp(os.path.join("results", variants_dir, "5_filter", "all_samples.vcf.gz")),
        vcfIndex = temp(os.path.join("results", variants_dir, "5_filter", "all_samples.vcf.gz.tbi")),
        detailMetrics = os.path.join("results", variants_dir, "5_filter", "log", "all_samples.variant_calling_detail_metrics"),
        summaryMetrics = os.path.join("results", variants_dir, "5_filter", "log", "all_samples.variant_calling_summary_metrics"),
    params:
        metricsBname = os.path.join("results", variants_dir, "5_filter", "log", "all_samples"),
    conda:
        "../envs/gatk.yml"
    resources:
        cpu = 1,
        ntasks = 2,
        mem_mb = 8000,
        time = "00-03:00:00",
    shell:
        """
        gatk \
            VariantFiltration \
            --R {input.refFa} \
            --V {input.vcf} \
            --window 35 \
            --cluster 3 \
            --filter-name "FS" --filter "FS > 60.0" \
            --filter-name "QD" --filter "QD < 2.0" \
            --filter-name "MQ" --filter "MQ < 40.0" \
            --filter-name "SOR" --filter "SOR > 4.0" \
            --filter-name "MQRankSum" --filter "MQRankSum < -12.5" \
            --filter-name "ReadPosRankSum" --filter "ReadPosRankSum < -8.0" \
            -O {output.vcf}

        gatk \
            CollectVariantCallingMetrics \
            --DBSNP {input.knownVariants} \
            --INPUT {output.vcf} \
            --OUTPUT {params.metricsBname}
        """

rule variants_select:
    input:
        unpack(knownVariants_files),
        vcf = rules.variants_filter.output.vcf,
        vcfIndex = rules.variants_filter.output.vcfIndex,
        refFa = rules.refs_downloadFa.output,
        refIndex = rules.refs_refIndex.output,
        refDict = rules.refs_refDict.output,
    output:
        vcf = os.path.join("results", variants_dir, "6_select", "all_samples.vcf.gz"),
        vcfIndex = os.path.join("results", variants_dir, "6_select", "all_samples.vcf.gz.tbi"),
        detailMetrics = os.path.join("results", variants_dir, "6_select", "log", "all_samples.variant_calling_detail_metrics"),
        summaryMetrics = os.path.join("results", variants_dir, "6_select", "log", "all_samples.variant_calling_summary_metrics"),
    params:
        metricsBname = os.path.join("results", variants_dir, "6_select", "log", "all_samples"),
    conda:
        "../envs/gatk.yml"
    resources:
        cpu = 1,
        ntasks = 1,
        mem_mb = 8000,
        time = "00-03:00:00",
    shell:
        """
        gatk \
            SelectVariants \
            --exclude-filtered \
            -R {input.refFa} \
            -V {input.vcf} \
            -O {output.vcf}

        gatk \
            CollectVariantCallingMetrics \
            --DBSNP {input.knownVariants} \
            --INPUT {output.vcf} \
            --OUTPUT {params.metricsBname}
        """