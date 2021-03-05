from snakemake.shell import shell

inputs = "".join("--INPUT {} ".format(f) for f in snakemake.input.vcf)

shell(
    """
    gatk \
    --java-options "-Xms2000m" \
    MergeVcfs \
    {inputs}\
    --OUTPUT {snakemake.output.vcf}
    """
)
