rule ase_setup:
    input:
        bams = expand(os.path.join("results", bqsr_dir, "bam", "{SAMPLE}.bam"), SAMPLE=samples),
        vcf = rules.variants_select.output.vcf
    output:
        snvDir = temp(directory(expand(os.path.join("results", wasp_dir, "1_snvs", "{SAMPLE}"), SAMPLE=samples))),
        intervals = temp(expand(os.path.join("results", aseRC_dir, "intervals", "{SAMPLE}.intervals"), SAMPLE=samples)),
    conda:
        "../envs/ase.yml"
    params:
        proj_root = os.getcwd(),
        variants_dir = variants_dir,
        wasp_dir = wasp_dir,
        aseRC_dir = aseRC_dir,
    resources:
        cpu = 1,
        ntasks = 1,
        mem_mb = 4000,
        time = "00-02:00:00",
    script:
        "../scripts/ase_setup.R"

