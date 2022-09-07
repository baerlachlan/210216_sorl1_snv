rule wasp_findIntersecting:
    input:
        bam = rules.bqsr_apply.output.bam,
        snvDir = os.path.join("results", wasp_dir, "1_snvs", "{SAMPLE}"),
    output:
        fq1 = temp(os.path.join("results", wasp_dir, "2_findIntersecting", "{SAMPLE}", "{SAMPLE}.remap.fq1.gz")),
        fq2 = temp(os.path.join("results", wasp_dir, "2_findIntersecting", "{SAMPLE}", "{SAMPLE}.remap.fq2.gz")),
        single = temp(os.path.join("results", wasp_dir, "2_findIntersecting", "{SAMPLE}", "{SAMPLE}.remap.single.fq.gz")),
        to_remap = temp(os.path.join("results", wasp_dir, "2_findIntersecting", "{SAMPLE}", "{SAMPLE}.to.remap.bam")),
        keep_intersect = temp(os.path.join("results", wasp_dir, "2_findIntersecting", "{SAMPLE}", "{SAMPLE}.keep.bam")),
    params:
        outDir = temp(directory(os.path.join("results", wasp_dir, "2_findIntersecting", "{SAMPLE}"))),
        wasp_scripts = "/hpcfs/users/a1647910/packages/WASP/mapping",
    conda:
        "../envs/ase.yml"
    resources:
        cpu = 1,
        ntasks = 1,
        mem_mb = 2000,
        time = "00-08:00:00",
    shell:
        """
        python {params.wasp_scripts}/find_intersecting_snps.py \
            --is_paired_end \
            --is_sorted \
            --output_dir {params.outDir} \
            --snp_dir {input.snvDir} \
            {input.bam}
        """

rule wasp_remap:
    input:
        fq1 = rules.wasp_findIntersecting.output.fq1,
        fq2 = rules.wasp_findIntersecting.output.fq2,
        starIndex = rules.refs_starIndex.output,
    output:
        remapped_unsorted = temp(os.path.join("results", wasp_dir, "3_remap", "{SAMPLE}Aligned.out.bam")),
        remapped_sorted = temp(os.path.join("results", wasp_dir, "3_remap", "{SAMPLE}sorted.out.bam")),
        index = temp(os.path.join("results", wasp_dir, "3_remap", "{SAMPLE}sorted.out.bam.bai")),
        STARgenome = temp(directory(os.path.join("results", wasp_dir, "3_remap", "{SAMPLE}_STARgenome"))),
        STARpass1 = temp(directory(os.path.join("results", wasp_dir, "3_remap", "{SAMPLE}_STARpass1"))),
    params:
        overhang = config["align"]["read_length"] - 1,
        bname = os.path.join("results", wasp_dir, "3_remap", "{SAMPLE}"),
        wasp_dir = os.path.join("results", wasp_dir)
    conda:
        "../envs/ase.yml"
    resources:
        cpu = 16,
        ntasks = 1,
        mem_mb = 32000,
        time = "00-04:00:00",
    shell:
        """
        STAR \
            --genomeDir {input.starIndex}\
            --runThreadN {resources.cpu} \
            --readFilesIn {input.fq1} {input.fq2} \
            --readFilesCommand "gunzip -c" \
            --sjdbOverhang {params.overhang} \
            --outSAMtype BAM Unsorted \
            --twopassMode Basic \
            --outFileNamePrefix {params.bname}

        mkdir -p 08_wasp/3_remap/log
        mv {params.bname}*out {params.wasp_dir}/log
        mv {params.bname}*tab {params.wasp_dir}/log

        samtools sort -o {output.remapped_sorted} {output.remapped_unsorted}
        samtools index {output.remapped_sorted}
        """

rule wasp_filterRemapped:
    input:
        to_remap = rules.wasp_findIntersecting.output.to_remap,
        remapped_unsorted = rules.wasp_remap.output.remapped_unsorted,
        remapped_sorted = rules.wasp_remap.output.remapped_sorted,
    output:
        keep_filter = temp(os.path.join("results", wasp_dir, "4_filterRemapped", "{SAMPLE}.keep.bam")),
    params:
        wasp_scripts = "/hpcfs/users/a1647910/packages/WASP/mapping",
    conda:
        "../envs/ase.yml"
    resources:
        cpu = 1,
        ntasks = 1,
        mem_mb = 24000,
        time = "00-02:00:00",
    shell:
        """
        python {params.wasp_scripts}/filter_remapped_reads.py \
            {input.to_remap} \
            {input.remapped_sorted} \
            {output.keep_filter}
        """

rule wasp_merge:
    input:
        keep_filter = rules.wasp_filterRemapped.output.keep_filter,
        keep_intersect = rules.wasp_findIntersecting.output.keep_intersect,
    output:
        keep_sorted = temp(os.path.join("results", wasp_dir, "5_merge", "{SAMPLE}.keep.merge.sort.bam")),
        keep_sortedIndex = temp(os.path.join("results", wasp_dir, "5_merge", "{SAMPLE}.keep.merge.sort.bam.bai")),
    params:
        keep_merged = temp(os.path.join("results", wasp_dir, "5_merge", "{SAMPLE}.keep.merge.bam")),
    conda:
        "../envs/ase.yml"
    resources:
        cpu = 1,
        ntasks = 1,
        mem_mb = 4000,
        time = "00-02:00:00",
    shell:
        """
        samtools merge {params.keep_merged} \
            {input.keep_filter} \
            {input.keep_intersect}
        samtools sort -o {output.keep_sorted} \
            {params.keep_merged}
        samtools index {output.keep_sorted}
        rm {params.keep_merged}
        """
