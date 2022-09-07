rule intervals:
    output:
        os.path.join("resources", "exons.intervals"),
    params:
        species = config["ref"]["species"],
        ensembl_release = config["ref"]["ensembl_release"],
    conda:
        "../envs/r.yml"
    resources:
        cpu = 1,
        ntasks = 1,
        mem_mb = 4000,
        time = "00-00:10:00",
    script:
        "../scripts/intervals.R"