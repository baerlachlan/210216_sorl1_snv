rule knownVariants_gvcf:
    input:
        bam = rules.splitNCigar.output.bam,
        bamIndex = rules.splitNCigar.output.bamIndex,
        refFa = rules.refs_downloadFa.output,
        refIndex = rules.refs_refIndex.output,
        refDict = rules.refs_refDict.output,
        intervals = rules.intervals.output,
    output:
        gvcf = temp(os.path.join("results", knownVariants_dir, "1_gvcf", "{SAMPLE}.g.vcf.gz")),
        gvcfIndex = temp(os.path.join("results", knownVariants_dir, "1_gvcf", "{SAMPLE}.g.vcf.gz.tbi")),
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
        """

rule knownVariants_sampleMap:
    input:
        expand(rules.knownVariants_gvcf.output.gvcf, SAMPLE=samples),
    output:
        os.path.join("results", knownVariants_dir, "sample_map.tsv"),
    params:
        rule_dir = os.path.join("results", knownVariants_dir),
    conda:
        "../envs/r.yml"
    resources:
        cpu = 1,
        ntasks = 1,
        mem_mb = 2000,
        time = "00-00:10:00",
    script:
        "../scripts/variants_sampleMap.R"

rule knownVariants_genomicsDB:
    input:
        gvcf = expand(rules.knownVariants_gvcf.output.gvcf, SAMPLE=samples),
        sample_map = rules.knownVariants_sampleMap.output,
        intervals = rules.intervals.output,
    output:
        genDB_dir = temp(dir(os.path.join("results", knownVariants_dir, "2_genomicsDB"))),
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

rule knownVariants_genotype:
    input:
        genDB_dir = rules.knownVariants_genomicsDB.output.genDB_dir,
        refFa = rules.refs_downloadFa.output,
        refIndex = rules.refs_refIndex.output,
        refDict = rules.refs_refDict.output,
    output:
        vcf = temp(os.path.join("results", knownVariants_dir, "3_jointGenotype", "known_variants.vcf.gz")),
        vcfIndex = temp(os.path.join("results", knownVariants_dir, "3_jointGenotype", "known_variants.vcf.gz.tbi")),
    params:
        gendb = os.path.join("gendb://" + "results", knownVariants_dir, "2_genomicsDB"),
    conda:
        "../envs/gatk.yml"
    resources:
        cpu = 1,
        ntasks = 2,
        mem_mb = 32000,
        time = "01-00:00:00",
    shell:
        """
        gatk --java-options "-Xmx16g" GenotypeGVCFs \
            -R {input.refFa} \
            -V {params.gendb} \
            -O {output.vcf}
        """

rule knownVariants_extract:
    input:
        vcf = rules.knownVariants_genotype.output.vcf,
        vcfIndex = rules.knownVariants_genotype.output.vcfIndex,
        refFa = rules.refs_downloadFa.output,
        refIndex = rules.refs_refIndex.output,
        refDict = rules.refs_refDict.output,
    output:
        vcf = temp(os.path.join("results", knownVariants_dir, "4_extract", "known_variants.vcf.gz")),
        vcfIndex = temp(os.path.join("results", knownVariants_dir, "4_extract", "known_variants.vcf.gz.tbi")),
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
        """

rule knownVariants_filter:
    input:
        vcf = rules.knownVariants_extract.output.vcf,
        vcfIndex = rules.knownVariants_extract.output.vcfIndex,
        refFa = rules.refs_downloadFa.output,
        refIndex = rules.refs_refIndex.output,
        refDict = rules.refs_refDict.output,
    output:
        vcf = temp(os.path.join("results", knownVariants_dir, "5_filter", "known_variants.vcf.gz")),
        vcfIndex = temp(os.path.join("results", knownVariants_dir, "5_filter", "known_variants.vcf.gz.tbi")),
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
        """

rule knownVariants_select:
    input:
        vcf = rules.knownVariants_filter.output.vcf,
        vcfIndex = rules.knownVariants_filter.output.vcfIndex,
        refFa = rules.refs_downloadFa.output,
        refIndex = rules.refs_refIndex.output,
        refDict = rules.refs_refDict.output,
    output:
        vcf = os.path.join("results", knownVariants_dir, "6_select", "known_variants.vcf.gz"),
        vcfIndex = os.path.join("results", knownVariants_dir, "6_select", "known_variants.vcf.gz.tbi"),
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
        """