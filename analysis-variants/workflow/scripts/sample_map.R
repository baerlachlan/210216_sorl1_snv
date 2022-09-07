library(readr)
library(magrittr)
library(stringr)
library(tibble)

## Allow user to run interactively or as part of snakemake workflow
if (exists("snakemake")) {
    variants_dir <- snakemake@params[["variants_dir"]]
    gvcf_path <- file.path("results", variants_dir, "1_gvcf")
    sample_map <- file.path("results", variants_dir, "sample_map.tsv")
} else {
    ## Edit these variables manually if running interactively
    variants_dir <- dirname(rstudioapi::getSourceEditorContext()$path) %>%
        str_remove("/workflow/scripts") %>%
        file.path("results", "07_variants")
    gvcf_path <- file.path(variants_dir, "1_gvcf")
    sample_map <- file.path(variants_dir, "sample_map.tsv")
}

files <- list.files(gvcf_path, pattern = ".vcf.gz$", full.names = TRUE)
samples <- basename(files) %>%
    str_remove(".g.vcf.gz")

if (!dir.exists(dirname(sample_map))) {
 dir.create(dirname(sample_map), recursive = TRUE)
}

# save.image("/hpcfs/users/a1647910/210216_sorl1_snv/analysis-variants/data.R")

tibble(sample = samples, file = files) %>%
    write_tsv(file = sample_map, col_names = FALSE)
