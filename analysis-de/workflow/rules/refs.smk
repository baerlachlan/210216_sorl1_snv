rule refs_downloadFa:
    output:
        temp(refFa_path),
    params:
        ref_url = ref_url,
        refGz_path = refGz_path
    threads:
        1
    shell:
        """
        rsync -avP {params.ref_url} {params.refGz_path}
        gzip -d {params.refGz_path}
        """

rule refs_downloadGtf:
    output:
        temp(gtf_path),
    params:
        gtf_url = gtf_url,
    threads:
        1
    shell:
        """
        rsync -avP {params.gtf_url} {output}
        """

rule refs_starIndex:
    input:
        ref_fa = rules.refs_downloadFa.output,
        gtf = rules.refs_downloadGtf.output,
    output:
        temp(directory(os.path.join("resources", "star"))),
    params:
        overhang = config["align"]["read_length"] - 1,
    conda:
        "../envs/align.yml"
    resources:
        cpu = 16,
        ntasks = 1,
        mem_mb = 32000,
        time = "00-01:30:00",
    shell:
        """
        zcat {input.gtf} > temp.gtf

        STAR \
            --runThreadN {resources.cpu} \
            --runMode genomeGenerate \
            --genomeDir {output} \
            --genomeFastaFiles {input.ref_fa} \
            --sjdbGTFfile temp.gtf \
            --sjdbOverhang {params.overhang}

        rm temp.gtf
        """