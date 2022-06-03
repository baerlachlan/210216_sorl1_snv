rule merge_fastq:
    input:
        fq = expand(
            os.path.join(
                "results", raw_dir, "fastq", "{{SAMPLE}}{MERGETAG}{{PAIRTAG}}" + config["fastq_ext"]
            ),
            MERGETAG=config["merge_fastq"]["tags"]
        ),
    output:
        fq = temp(os.path.join(
            "results", merge_dir, "fastq", "{SAMPLE}{PAIRTAG}" + config["fastq_ext"]
        )),
    resources:
        cpu = 1,
        ntasks = 1,
        mem_mb = 8000,
        time = "00-00:30:00",
    shell:
        """
        cat {input.fq} > {output.fq}
        """