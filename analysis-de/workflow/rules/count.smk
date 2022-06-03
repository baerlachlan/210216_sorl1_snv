rule count:
    input:
        bam = expand(os.path.join("results", align_dir, "bam", "{SAMPLE}.bam"), SAMPLE = samples),
        gtf = rules.refs_downloadGtf.output,
    output:
        counts = os.path.join("results", counts_dir, "counts.out"),
        counts_summary = os.path.join("results", counts_dir, "counts.out.summary"),
    params:
        Q = config["count"]["min_qual"],
        s = config["count"]["strandedness"],
        min_overlap = config["count"]["min_overlap"],
        frac_overlap = config["count"]["frac_overlap"],
        extra = config["count"]["extra"],
    conda:
        "../envs/count.yml"
    resources:
        cpu = 4,
        ntasks = 1,
        mem_mb = 8000,
        time = "00-02:00:00",
    shell:
        """
        featureCounts \
            {params.extra} \
            -Q {params.Q} \
            -s {params.s} \
            --minOverlap {params.min_overlap} \
            --fracOverlap {params.frac_overlap} \
            -T {resources.cpu} \
            -a {input.gtf} \
            -o {output.counts} \
            {input.bam}
        """