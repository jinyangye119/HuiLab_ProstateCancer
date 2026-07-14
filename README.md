# HuiLab Prostate Cancer

Analysis code for the Hui Lab prostate cancer project, including final figure generation and cell deconvolution workflows.

## Repository Contents

```text
Code/
  Edec_deconvolution.R              # EDec methylation deconvolution parameter sweep
  Final_Figures.Rmd                 # Final figure and downstream analysis notebook
  xCell_deconvolution.Rmd           # Build a custom prostate xCell2 reference
  xCellv_deconvolution_JHU.Rmd      # JHU xCell2 reference-building workflow
  download_deconvolution_reference.R
```

Large reference files are not stored in GitHub. They are archived on Zenodo:

- Dataset: HuiLab_CellDeconvolution
- DOI: https://doi.org/10.5281/zenodo.21362490
- Record: https://zenodo.org/records/21362490
- Archive: `Deconvolution_reference.zip`
- MD5: `82ac3788f757c720125161f42e08fca1`

## Reference Data

Download and unpack the reference files into `Code/Deconvolution_reference/`:

```bash
Rscript Code/download_deconvolution_reference.R
```

After extraction, the expected layout is:

```text
Code/Deconvolution_reference/
  eDec_reference/
    reference_dataset_list.csv
    reference_dataset_list_strict_luminal.csv
    reference_dataset_manifest_wide.csv
    reference_dataset_manifest_wide_strict_luminal.csv
    prostate_*_beta_matrix.tsv
    prostate_*_metadata.csv
    GSE*_series_matrix.txt.gz
  xCell_reference/
    dge_combined_annotated.rds
    prostate_custom_xCell2_reference.rds
```

## R Environment

The workflows use R and Bioconductor packages, including:

```r
data.table
matrixStats
EDec
ggpubr
edgeR
tidyverse
ggplot2
DESeq2
openxlsx
org.Hs.eg.db
ggrepel
reshape2
RColorBrewer
DEP
readxl
pheatmap
VennDiagram
limma
AnnotationDbi
Matrix
Seurat
xCell2
biomaRt
BiocParallel
```

Install missing Bioconductor packages with `BiocManager::install()` and CRAN packages with `install.packages()`.

## Running The Workflows

Run the xCell2 reference-building notebooks from RStudio or render them with `rmarkdown::render()`:

```r
rmarkdown::render("Code/xCell_deconvolution.Rmd")
rmarkdown::render("Code/xCellv_deconvolution_JHU.Rmd")
```

Run the EDec methylation deconvolution workflow from the repository root:

```bash
Rscript Code/Edec_deconvolution.R
```

Optional environment variables:

```bash
EDEC_METH_FILE=/path/to/methylation_matrix.csv \
EDEC_OUTPUT_TAG=run_name \
Rscript Code/Edec_deconvolution.R
```

Generate final figures by rendering:

```r
rmarkdown::render("Code/Final_Figures.Rmd")
```

## Notes

Some scripts contain absolute local paths from the original analysis environment. Before running on a different machine, update the path settings near the top of each script/notebook so they point to your local project directory, methylation matrix, metadata files, output directory, and the downloaded `Code/Deconvolution_reference/` folder.

Please cite the Zenodo dataset when using the deconvolution reference files:

Jin, Yang. HuiLab_CellDeconvolution. Zenodo. https://doi.org/10.5281/zenodo.21362490
