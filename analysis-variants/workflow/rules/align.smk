rule align:
    input:
        R1 = rules.trim.output.R1,
        starIndex = rules.refs_starIndex.output,
    output:
        bam = temp(os.path.join("results", align_dir, "bam", "{SAMPLE}{MERGETAG, .*}.bam")),
        bamIndex = temp(os.path.join("results", align_dir, "bam", "{SAMPLE}{MERGETAG, .*}.bam.bai")),
        STARgenome = temp(directory(os.path.join("results", align_dir, "bam", "{SAMPLE}{MERGETAG, .*}_STARgenome"))),
        STARpass1 = temp(directory(os.path.join("results", align_dir, "bam", "{SAMPLE}{MERGETAG, .*}_STARpass1"))),
    params:
        overhang = config["align"]["read_length"] - 1,
        bname = os.path.join("results", align_dir, "bam", "{SAMPLE}{MERGETAG, .*}"),
        bamUnsorted = os.path.join("results", align_dir, "bam", "{SAMPLE}{MERGETAG, .*}Aligned.out.bam"),
        align_dir = os.path.join("results", align_dir),
    conda:
        "../envs/gatk.yml"
    resources:
        cpu = 16,
        ntasks = 1,
        mem_mb = 32000,
        time = "00-02:00:00",
    shell:
        """
        STAR \
            --genomeDir {input.starIndex}\
            --runThreadN {resources.cpu} \
            --readFilesIn {input.R1} \
            --readFilesCommand "gunzip -c" \
            --sjdbOverhang {params.overhang} \
            --outSAMtype BAM Unsorted \
            --twopassMode Basic \
            --outFileNamePrefix {params.bname}

        samtools sort {params.bamUnsorted} > {output.bam}
        samtools index {output.bam}
        rm {params.bamUnsorted}

        mkdir -p {params.align_dir}/log
        mv {params.bname}*out {params.align_dir}/log
        mv {params.bname}*tab {params.align_dir}/log
        """