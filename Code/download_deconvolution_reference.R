#!/usr/bin/env Rscript

zenodo_url <- "https://zenodo.org/api/records/21362490/files/Deconvolution_reference.zip/content"
expected_md5 <- "82ac3788f757c720125161f42e08fca1"

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[[1]]) else "Code/download_deconvolution_reference.R"
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
code_dir <- file.path(repo_root, "Code")
zip_file <- file.path(code_dir, "Deconvolution_reference.zip")
target_dir <- file.path(code_dir, "Deconvolution_reference")

message("Downloading deconvolution reference data from Zenodo...")
dir.create(code_dir, recursive = TRUE, showWarnings = FALSE)
download.file(zenodo_url, zip_file, mode = "wb", quiet = FALSE)

observed_md5 <- tools::md5sum(zip_file)[[1]]
if (!identical(unname(observed_md5), expected_md5)) {
  stop(
    "MD5 checksum mismatch for downloaded Zenodo archive.\n",
    "Expected: ", expected_md5, "\n",
    "Observed: ", observed_md5
  )
}

if (dir.exists(target_dir)) {
  message("Reference directory already exists; files will be overwritten if present.")
}

utils::unzip(zip_file, exdir = code_dir)
message("Reference data ready at: ", target_dir)
