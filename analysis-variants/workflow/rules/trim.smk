rule trim:
    input:
        unpack(trim_inputs_se),
    output:
        R1 = temp(os.path.join("results", trim_dir, "fastq", "{SAMPLE}{MERGETAG, .*}" + config["fastq_ext"])),
        html = os.path.join("results", trim_dir, "log", "{SAMPLE}{MERGETAG, .*}.html"),
    conda:
        "../envs/gatk.yml"
    resources:
        cpu = 1,
        ntasks = 2,
        mem_mb = 2000,
        time = "00-02:00:00",
    shell:
        """
        fastp \
            -i {input.R1}  \
            -o {output.R1} \
            --qualified_quality_phred 20 \
            --length_required 35 \
            --trim_poly_g \
            --thread 1 \
            --html {output.html} \
            --json /dev/null \
        """