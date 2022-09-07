####
## Manage config
####

if not config["umi"]["activate"]:
    config["umi"]["add_header"]["activate"] = False
if not config["umi"]["add_header"]["activate"]:
    config["umi"]["add_header"]["tag"] = [""]
if not config["paired_end"]["activate"]:
    config["paired_end"]["tags"] = [""]
if not config["merge_samples"]["activate"]:
    config["merge_samples"]["tags"] = [""]

####
## Directory structure
####

## Define list of all possible directories in workflow
dir_list = [
    "rawData", "addUmis", "trim", "align", "addRG", "groupUmis", "mergeSamples",
    "markDuplicates", "splitNCigar", "knownVariants", "bqsr", "variants",
    "wasp", "aseRC",
]
## Remove directories not required by workflow config
if not config["umi"]["add_header"]["activate"]:
    dir_list.remove("addUmis")
if not config["umi"]["activate"]:
    dir_list.remove("groupUmis")
if not config["merge_samples"]["activate"]:
    dir_list.remove("mergeSamples")
if not config["bootstrap_known_variants"]["activate"]:
    dir_list.remove("knownVariants")
if not config["ase_counts"]["activate"]:
    dir_list.remove("wasp")
    dir_list.remove("aseRC")
## Generate the directory ordering structure from list indices
glob_vars = vars()
for dir_name in dir_list:
    glob_vars[dir_name + "_dir"] = str(dir_list.index(dir_name)).zfill(2) + "_" + dir_name

####
## Sample names
####

if config["samples_tsv"]:
    samples_df = pd.read_csv(config["samples_tsv"], sep="\t")
    samples = list(dict.fromkeys(samples_df["sample"])) ## Remove duplicates
else:
    samples = os.listdir(os.path.join("results", rawData_dir, "fastq"))
    samples = [
        sample.replace(config["fastq_ext"], "") for sample in samples
    ]
    for tag in config["paired_end"]["tags"]:
        samples = [sample.replace(tag, "") for sample in samples]
    if config["merge_samples"]["activate"]:
        for tag in config["merge_samples"]["tags"]:
            samples = [sample.replace(tag, "") for sample in samples]
    if config["umi"]["add_header"]["activate"]:
        samples = [sample.replace(config["umi"]["add_header"]["tag"], "") for sample in samples]
    samples = list(dict.fromkeys(samples))  ## Remove duplicates

####
## Reference filenames, paths and urls for downloading
####

refGz = os.path.join(
    ".".join([
        config["ref"]["species"].capitalize(),
        config["ref"]["assembly"],
        "dna.primary_assembly.fa.gz"
    ])
)
refGz_path = os.path.join("resources", refGz)
refFa = refGz.rstrip(".gz")
refFa_path = os.path.join("resources", refFa)
ref_url = os.path.join(
    "rsync://ftp.ensembl.org/ensembl/pub",
    "release-" + str(config["ref"]["ensembl_release"]),
    "fasta",
    config["ref"]["species"],
    "dna",
    refGz
)
gtf = os.path.join(
    ".".join([
        config["ref"]["species"].capitalize(),
        config["ref"]["assembly"],
        str(config["ref"]["ensembl_release"]),
        "chr.gtf.gz"
    ])
)
gtf_path = os.path.join("resources", gtf)
gtf_url = os.path.join(
    "rsync://ftp.ensembl.org/ensembl/pub",
    "release-" + str(config["ref"]["ensembl_release"]),
    "gtf",
    config["ref"]["species"],
    gtf
)
knownVariants = os.path.join(".".join([config["ref"]["species"], "vcf.gz"]))
knownVariants_path = os.path.join("resources", knownVariants)
knownVariants_url = os.path.join(
    "rsync://ftp.ensembl.org/ensembl/pub",
    "release-" + str(config["ref"]["ensembl_release"]),
    "variation/vcf",
    config["ref"]["species"],
    knownVariants
)

####
## Input functions for rules affected by optional rules
####

def trim_inputs(wildcards):
    return {
        "R1": os.path.join(
            "results",
            addUmis_dir if config["umi"]["add_header"]["activate"] else rawData_dir,
            "fastq",
            wildcards.SAMPLE + wildcards.MERGETAG + config["paired_end"]["tags"][0] + config["fastq_ext"]
        ),
        "R2": os.path.join(
            "results",
            addUmis_dir if config["umi"]["add_header"]["activate"] else rawData_dir,
            "fastq",
            wildcards.SAMPLE + wildcards.MERGETAG + config["paired_end"]["tags"][1] + config["fastq_ext"]
        ),
    }

def trim_inputs_se(wildcards):
    return {
        "R1": os.path.join(
            "results",
            addUmis_dir if config["umi"]["add_header"]["activate"] else rawData_dir,
            "fastq",
            wildcards.SAMPLE + wildcards.MERGETAG + config["fastq_ext"]
        ),
    }

def mergeSamples_inputs(wildcards):
    return {
        "bam" : expand(os.path.join(
            "results",
            groupUmis_dir if config["umi"]["activate"] else addRG_dir,
            "bam",
            wildcards.SAMPLE + "{MERGETAG}" + ".bam"
        ), MERGETAG = config["merge_samples"]["tags"]),
        "bamIndex" : expand(os.path.join(
            "results",
            groupUmis_dir if config["umi"]["activate"] else addRG_dir,
            "bam",
            wildcards.SAMPLE + "{MERGETAG}" + ".bam.bai"
        ), MERGETAG = config["merge_samples"]["tags"])
    }

def markDuplicates_inputs(wildcards):
    if config["merge_samples"]["activate"]:
        input_dir = mergeSamples_dir
    elif config["umi"]["activate"]:
        input_dir = groupUmis_dir
    else:
        input_dir = addRG_dir
    return {
        "bam" : os.path.join(
            "results",
            input_dir,
            "bam",
            wildcards.SAMPLE + ".bam"
        ),
        "bamIndex" : os.path.join(
            "results",
            input_dir,
            "bam",
            wildcards.SAMPLE + ".bam.bai"
        ),
    }

def knownVariants_files(wildcards):
    if config["bootstrap_known_variants"]["activate"]:
        return {
            "knownVariants" : os.path.join(
                "results",
                knownVariants_dir,
                "6_select",
                "known_variants.vcf.gz"
            ),
            "knownVariantsIndex" : os.path.join(
                "results",
                knownVariants_dir,
                "6_select",
                "known_variants.vcf.gz.tbi"
            )
        }
    else:
        return {
            "knownVariants" : knownVariants_path,
            "knownVariantsIndex" : knownVariants_path + ".tbi",
            }

####
## All file endpoints
####

def workflow_outputs():

    outputs = []

    ## FastQC reports
    fqc_raw = expand(
        os.path.join("results", rawData_dir, "FastQC/{SAMPLE}{MERGETAG}{PAIRTAG}_fastqc.{EXT}"),
        SAMPLE=samples,
        MERGETAG=config["merge_samples"]["tags"],
        PAIRTAG=config["paired_end"]["tags"],
        EXT=["html", "zip"]
    )
    outputs.extend(fqc_raw)
    fqc_trim = expand(
        os.path.join("results", trim_dir, "FastQC/{SAMPLE}{MERGETAG}{PAIRTAG}_fastqc.{EXT}"),
        SAMPLE=samples,
        MERGETAG=config["merge_samples"]["tags"],
        PAIRTAG=config["paired_end"]["tags"],
        EXT=["html", "zip"]
    )
    outputs.extend(fqc_trim)
    fqc_align = expand(
        os.path.join("results", align_dir, "FastQC/{SAMPLE}{MERGETAG}_fastqc.{EXT}"),
        SAMPLE=samples,
        MERGETAG=config["merge_samples"]["tags"],
        EXT=["html", "zip"]
    )
    outputs.extend(fqc_align)

    ## BQSR summary for QC purposes
    bqsr_analyzeCovariates = expand(
        os.path.join("results", bqsr_dir, "recal/{SAMPLE}.analyzeCovariates.csv"),
        SAMPLE=samples
    )
    outputs.extend(bqsr_analyzeCovariates)

    ## Variants
    variants = [os.path.join("results", variants_dir, "6_select/all_samples.vcf.gz")]
    outputs.extend(variants)

    ## ASE read counts
    if config["ase_counts"]["activate"]:
        aseRC = expand(
            os.path.join("results", aseRC_dir, "{DIR}", "{SAMPLE}.tsv"),
            DIR=["wasp", "no_wasp"],
            SAMPLE=samples
        )
        outputs.extend(aseRC)

    return(outputs)