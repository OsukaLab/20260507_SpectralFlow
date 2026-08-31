# Spectral Flow Cytometry UMAP and FlowSOM Analysis

This repository contains the R Markdown workflow used for spectral flow cytometry analysis across three experimental comparisons:

- **GS vs GSRR**
- **Collagen knockout (COL KO) vs control**
- **BKM120 vs control**

The workflow analyzes pre-gated spectral flow cytometry populations using **UMAP** for visualization and **FlowSOM** for unsupervised clustering. Depending on the experiment, separate analysis sections are provided for T-cell/T-NK-cell and myeloid populations.

The general workflow is:

1. Import compensated, pre-gated FCS files.
2. Parse treatment group and mouse/sample identity from filenames.
3. Rename fluorescence channels to biological marker names.
4. Apply an arcsinh transformation.
5. Downsample equally across biological samples.
6. Generate a shared UMAP containing all experimental groups.
7. Cluster cells with FlowSOM.
8. Visualize marker expression and cluster phenotypes.
9. Merge phenotypically similar FlowSOM clusters when appropriate.
10. Calculate cluster proportions by individual sample and experimental group.
11. Save publication-ready UMAPs, heatmaps, and CSV tables.

---

## Important: what this code expects

This pipeline is **not intended to start from completely ungated raw cytometry files**.

The analysis assumes that:

- compensation has already been performed upstream;
- the desired parent population has already been gated and exported from FlowJo;
- each biological sample is exported as an individual `.fcs` file;
- files from the same analysis use the same antibody panel and channel assignments;
- filenames contain enough information to identify the experimental group and mouse/sample number.

The scripts import FCS files with:

```r
transformation = FALSE
truncate_max_range = FALSE
```

and use compensated fluorescence parameters such as `FJComp-B510-A`. Therefore, **do not compensate the data a second time in R**.

The exact upstream gate depends on the analysis section. For example, T-cell sections use gated lymphocyte/T-cell populations, whereas myeloid sections use the corresponding gated myeloid population.

---

## Repository file

The main analysis is contained in the R Markdown file:

```text
20260811_Manuscript.Rmd
```

The file contains multiple experiment-specific analysis sections. These sections reuse some object names, so the safest way to use the repository is to:

> **Run one experiment and cell-compartment section at a time, from top to bottom.**

Do not assume that running the entire R Markdown document at once will reproduce every analysis cleanly in the same R session.

For a completely independent analysis, restarting the R session before beginning a new experiment is recommended.

---

## Software requirements

R and RStudio are recommended.

Core packages used by the workflow include:

- `tidyverse`
- `dplyr`
- `ggplot2`
- `uwot`
- `FlowSOM`
- `flowCore`
- `pheatmap`

Install CRAN packages with:

```r
install.packages(c(
  "tidyverse",
  "dplyr",
  "ggplot2",
  "uwot",
  "pheatmap"
))
```

Install the Bioconductor packages with:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "FlowSOM",
  "flowCore"
))
```

Then load the packages used by the spectral flow sections:

```r
library(tidyverse)
library(uwot)
library(ggplot2)
library(dplyr)
library(FlowSOM)
library(pheatmap)
library(flowCore)
```

---

## Recommended input directory structure

The original script uses absolute paths beginning with `/Data/`. You do **not** need to reproduce that exact directory structure, but you must change the input and output paths in the R Markdown file to match your computer.

A convenient organization is:

```text
project/
├── README.md
├── 20260811_Manuscript.Rmd
├── data/
│   ├── GSvsGSRR/
│   │   ├── Tcell/
│   │   └── Myeloid/
│   ├── COLvsCOLKO/
│   │   ├── Tcell/
│   │   └── Myeloid/
│   └── BKM120/
│       ├── Tcell/
│       └── Myeloid/
└── output/
    ├── GSvsGSRR/
    ├── COLvsCOLKO/
    └── BKM120/
```

Then replace paths such as:

```r
fcs_dir <- "/Data/Input/SpectralFlow_GSvsGSRR/Tcell"
```

with your local path, for example:

```r
fcs_dir <- "data/GSvsGSRR/Tcell"
```

Do the same for each `out_dir`, `fcs_dir`, or `input_dir` definition in the section you are running.

---

# How to use the code

## 1. Export compensated and gated FCS files

Perform the required compensation and gating in FlowJo before exporting.

Export **one FCS file per biological sample** for the population being analyzed.

Do not combine all mice into a single FCS file. The R workflow uses individual files so that it can preserve biological replicate identity, downsample equally per mouse, and calculate cluster proportions per mouse.

---

## 2. Check the filenames

Sample metadata are generated directly from filenames.

The original analyses use filename patterns that distinguish the experimental groups.

Examples include:

### GS vs GSRR

The filename parser searches for:

```text
_GS_
_GSRR_
```

and extracts the mouse number from the filename.

### COL KO vs control

The T-cell analysis searches for:

```text
COL_KO
_Ctrl_
```

and the original mouse-number parser expects a three-digit identifier in the relevant portion of the filename.

### BKM120 vs control

The parser searches for:

```text
_BKM_
_Ctrl_
```

and extracts a three-digit mouse/sample identifier.

After the `sample_meta` object is created, **always inspect it before continuing**:

```r
sample_meta
table(sample_meta$group)
```

or the corresponding experiment-specific metadata object.

If the script contains checks such as:

```r
stopifnot(!any(is.na(sample_meta$group)))
stopifnot(!any(is.na(sample_meta$mouse)))
```

they should pass before continuing.

### Using different filenames

If your filenames do not follow the original naming convention, edit the `case_when()`, `str_detect()`, `str_extract()`, or `str_match()` commands used to build the metadata.

For example:

```r
group = case_when(
  str_detect(filename, "_Treatment_") ~ "Treatment",
  str_detect(filename, "_Ctrl_") ~ "Ctrl",
  TRUE ~ NA_character_
)
```

Do not continue until every file has a valid group and sample ID.

---

## 3. Inspect the FCS channel names

Before importing the full experiment, the workflow can read a test FCS file and display the parameter names:

```r
test_fcs <- read.FCS(
  sample_meta$file[1],
  transformation = FALSE,
  truncate_max_range = FALSE
)

pData(parameters(test_fcs))[, c("name", "desc")]
```

Use this output to confirm that your FCS channel names match the `rename()` block in the script.

If the names are different, update the channel mapping before continuing.

---

# Marker/channel mapping used in the original analyses

## T-cell/T-NK panel

The original compensated channel mapping is:

| Marker | FCS parameter |
|---|---|
| PD1 | `FJComp-B510-A` |
| CD44 | `FJComp-B675-A` |
| CTLA4 | `FJComp-B710-A` |
| CD3 | `FJComp-R675-A` |
| CD45 | `FJComp-R780-A` |
| Live/dead | `FJComp-UV379-A` |
| CD4 | `FJComp-UV660-A` |
| CD62L | `FJComp-V427-A` |
| KI67 | `FJComp-V510-A` |
| CD8 | `FJComp-V595-A` |
| NK1.1 | `FJComp-V660-A` |
| CD19 | `FJComp-V710-A` |
| CD69 | `FJComp-V785-A` |
| mCherry | `FJComp-YG585-A` |
| FOXP3 | `FJComp-YG780-A` |

These are renamed in R to simplified variables such as `PD1`, `CD4`, `CD8`, `NK11`, and `FOXP3`.

If you use a different panel or cytometer configuration, change the `rename()` block to match your FCS parameters.

---

## Myeloid panel

The original myeloid channel mapping is:

| Marker | FCS parameter |
|---|---|
| LY6C | `FJComp-B510-A` |
| CD172 | `FJComp-B710-A` |
| XCR1 | `FJComp-R675-A` |
| CD45 | `FJComp-R780-A` |
| Live/dead | `FJComp-UV379-A` |
| CD317 | `FJComp-UV585-A` |
| CD11B | `FJComp-UV660-A` |
| CD206 | `FJComp-V427-A` |
| CCR2 | `FJComp-V510-A` |
| LY6G | `FJComp-V595-A` |
| F480 | `FJComp-V660-A` |
| MHCII | `FJComp-V710-A` |
| CD86 | `FJComp-V785-A` |
| mCherry | `FJComp-YG585-A` |
| CD11C | `FJComp-YG780-A` |

Again, these mappings must be changed if the channels in your FCS files differ.

---

## 4. Import all FCS files

The script imports each FCS file with `flowCore::read.FCS()` and converts the expression matrix into a data frame.

Sample-level metadata are then added to every event:

- filename
- experimental group
- mouse/sample identifier
- unique sample ID

All samples are combined into one data frame using `bind_rows()`.

Immediately inspect the cell counts:

```r
cell_counts
```

The purpose of this step is to identify unexpectedly small samples and confirm that all biological replicates were imported.

---

## 5. Apply any analysis-specific cell filtering

Some T-cell/T-NK sections contain an additional filter that removes cells negative for CD4, CD8, and NK1.1 using an exploratory expression cutoff.

For example, the logic is:

```r
CD4 > 0 | CD8 > 0 | NK11 > 0
```

This filtering step is specific to the original dataset and gating strategy.

If you are applying the code to a new experiment, decide whether this filter is biologically and technically appropriate for your data rather than automatically retaining it.

---

## 6. Define the markers used for dimensionality reduction

Each analysis section contains a vector defining the markers that drive UMAP and/or FlowSOM.

For example:

```r
markers <- c(
  "CD4",
  "CD8",
  "FOXP3",
  "PD1",
  "CTLA4",
  "CD44",
  "CD62L",
  "KI67",
  "CD69",
  "NK11"
)
```

The exact marker list differs between T-cell and myeloid analyses.

Only use markers that are biologically appropriate for the population being analyzed and that are present in every sample.

Markers used only for visualization can be transformed and plotted without necessarily being included as UMAP-driving variables.

---

## 7. Arcsinh-transform the fluorescence values

The manuscript workflow uses an arcsinh transformation with:

```r
cofactor <- 150
```

and creates transformed columns using:

```r
asinh(expression / cofactor)
```

For example:

```r
data <- data %>%
  mutate(
    across(
      all_of(markers),
      ~ asinh(.x / cofactor),
      .names = "{.col}_asinh"
    )
  )
```

The transformed variables are then used for UMAP, FlowSOM, feature plots, and cluster phenotype summaries.

If your data were acquired or exported differently, the most appropriate cofactor may differ. To reproduce the repository analyses, retain the cofactor encoded in the relevant analysis section.

---

## 8. Downsample equally across mice/samples

The workflow downsampled cells **within each biological sample** before constructing the shared UMAP.

A typical pattern is:

```r
data_ds <- data %>%
  group_by(sample_id) %>%
  slice_sample(n = downsample_n) %>%
  ungroup()
```

Some sections calculate the usable target as the smaller of a preset target and the smallest sample:

```r
downsample_n <- min(
  5000,
  min(cell_counts$n_cells)
)
```

Equal downsampling prevents samples with very large event counts from dominating the UMAP and FlowSOM solution.

Always confirm the result:

```r
data_ds %>%
  count(group, sample_id)
```

Every sample should contain the intended number of cells after downsampling.

---

## 9. Generate the shared UMAP

UMAP is generated using `uwot::umap()` on the transformed marker matrix.

The exact UMAP settings are defined within each experiment-specific code block.

To reproduce the repository analysis:

- retain the marker set used in that section;
- retain the encoded `n_neighbors`;
- retain the encoded `min_dist`;
- retain the Euclidean metric;
- retain the random seed.

The workflow uses:

```r
set.seed(123)
```

before stochastic steps.

### Important: treatment groups share the same UMAP

The treatment and control cells are embedded **together in a single UMAP**.

Group-specific or split UMAPs are subsequently created by plotting subsets of the same UMAP coordinates.

Do **not** independently recalculate separate UMAPs for treatment and control if the goal is to compare their distributions in the same embedding.

---

## 10. Inspect feature UMAPs before interpreting clusters

The scripts generate feature UMAPs for individual markers.

These plots are important for determining the phenotype of different regions of the UMAP and for interpreting FlowSOM clusters.

Typical outputs include markers such as:

### T cells

- CD4
- CD8
- FOXP3
- PD1
- CTLA4
- CD44
- CD62L
- KI67
- CD69
- NK1.1
- CD3

### Myeloid cells

- CD11B
- CD11C
- XCR1
- LY6C
- CD172
- CD206
- CCR2
- LY6G
- F480
- MHCII
- CD86
- CD317

Feature maps should be interpreted together with the cluster heatmap rather than using UMAP location alone to assign cell identity.

---

## 11. Run FlowSOM clustering

FlowSOM is run on the transformed marker matrix used for the relevant cell population.

The general structure is:

```r
fsom <- FlowSOM::FlowSOM(
  input_matrix,
  compensate = FALSE,
  transform = FALSE,
  scale = FALSE,
  colsToUse = colnames(input_matrix),
  xdim = ...,
  ydim = ...,
  nClus = ...
)
```

The exact SOM grid size and requested number of metaclusters differ between analysis sections.

For reproduction of the repository figures, use the settings already encoded in the relevant section.

For a new dataset, the appropriate number of clusters should be evaluated rather than copied automatically.

---

## 12. Examine cluster phenotypes

After FlowSOM clustering, median transformed marker expression is calculated for each cluster.

These values are visualized with `pheatmap`.

The heatmaps generally use marker-wise scaling so that the relative expression of each marker can be compared across clusters.

Use the following together when assigning a cluster phenotype:

1. cluster location on the UMAP;
2. marker feature UMAPs;
3. median-expression heatmap;
4. known biology of the gated parent population.

Do not assign a phenotype based only on physical proximity between clusters on the UMAP.

---

## 13. Understand manual cluster merging

Several manuscript analyses intentionally merge FlowSOM metaclusters after inspection of their phenotypes.

These merges are encoded directly in the corresponding analysis section.

### If you are reproducing the original dataset

Keep the final merge definitions used by the script.

### If you are analyzing new data

Do **not** blindly reuse the same cluster numbers or merge definitions.

FlowSOM cluster numbers are arbitrary labels and can change when:

- the input dataset changes;
- the marker set changes;
- the number of clusters changes;
- the random seed changes;
- package versions change;
- filtering or downsampling changes.

For new data, first inspect the UMAP, feature maps, and cluster median-expression heatmap, then define biologically justified merges for that dataset.

---

## 14. Calculate cluster proportions

The downstream workflow calculates cluster representation at the biological-sample level.

The general logic is:

```r
cluster_proportions <- data %>%
  count(group, sample_id, cluster) %>%
  group_by(group, sample_id) %>%
  mutate(
    proportion = n / sum(n)
  ) %>%
  ungroup()
```

Group-level summaries are then calculated from the individual biological replicates.

The biological replicate is the **mouse/sample**, not the individual cell.

Cluster proportions should therefore be interpreted as the fraction of the analyzed parent population assigned to a cluster for each mouse.

They are not absolute tumor cell counts unless absolute counting information was measured separately.

---

## 15. Save final figures and tables

The final save blocks generate publication-ready PDF figures and CSV tables.

Depending on the analysis section, outputs include:

- combined gray UMAP;
- individual marker feature UMAPs;
- FlowSOM cluster UMAP;
- treatment/control split UMAP;
- treatment-only and control-only UMAPs;
- cluster phenotype heatmap;
- cluster median-expression table;
- cluster proportions for each biological sample;
- group-level cluster-proportion summaries;
- statistical comparison tables where included.

Output directories are created automatically with:

```r
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
```

Change `out_dir` to the location where you want your results saved.

---

# Reproducing the repository analysis vs adapting the workflow

| Setting | Reproduce repository analysis | Analyze a new dataset |
|---|---|---|
| Input directory | Change only to the location of the original FCS files | Change |
| Filename parser | Keep if original filenames are unchanged | Usually change |
| Channel `rename()` mapping | Keep | Change if channels differ |
| Marker list | Keep | Review/change |
| Arcsinh cofactor | Keep section-specific value | Validate for your data |
| Downsampling | Keep section-specific value | Choose based on sample size |
| Random seed | Keep | Keeping `123` is recommended |
| UMAP parameters | Keep section-specific settings | Tune only if justified |
| FlowSOM settings | Keep | Evaluate |
| Cluster numbers | Keep | Re-identify |
| Manual cluster merges | Keep | Rebuild from your own phenotypes |
| Output path | Change as needed | Change |

---

# Recommended order for adapting the code to a new experiment

If you are using this repository as a template for a new spectral flow experiment, make changes in this order:

1. Copy the analysis section most similar to your experiment.
2. Change the input directory.
3. Change filename parsing.
4. Confirm all samples receive a valid group and sample ID.
5. Inspect FCS parameter names.
6. Change the `rename()` channel mapping.
7. Define the markers appropriate for your panel.
8. Run the transformation.
9. Inspect sample cell counts.
10. Choose an equal downsampling target.
11. Generate the shared UMAP.
12. Examine feature UMAPs.
13. Run FlowSOM.
14. Inspect cluster medians and the phenotype heatmap.
15. Determine whether any clusters should be merged.
16. Calculate per-sample cluster proportions.
17. Generate split/group UMAPs.
18. Save the final PDFs and CSV files.

---

# Reproducibility notes

## Random seed

The analyses use:

```r
set.seed(123)
```

before stochastic operations such as downsampling, UMAP generation, FlowSOM, and plotting-order randomization.

Keep these seeds unchanged when attempting to reproduce the original analysis.

---

## UMAP orientation

The absolute direction of the UMAP axes has no biological meaning.

A UMAP may occasionally appear rotated, reflected, or flipped between software/package versions even when the underlying neighborhood structure is similar.

Interpret cluster relationships and marker distributions rather than the sign or orientation of `UMAP1` and `UMAP2`.

---

## Cluster numbering

FlowSOM cluster numbers are identifiers, not biological rankings.

For example, cluster 2 is not inherently more or less related to cluster 3 than to cluster 8.

Cell identity must be assigned using marker expression and biological context.

---

## Parent-population dependence

This workflow begins with pre-gated populations.

Therefore, all UMAPs, FlowSOM clusters, and cluster proportions describe the **population that was exported into the analysis**, not the entire original tumor or tissue sample.

The upstream gating strategy must be considered when interpreting downstream proportions.

---

# Common problems

## No FCS files are found

Check:

```r
length(fcs_files)
basename(fcs_files)
```

If the length is zero, the directory path or file extension is incorrect.

---

## `rename()` says a column does not exist

Your FCS parameter names do not match the original panel.

Run:

```r
pData(parameters(test_fcs))[, c("name", "desc")]
```

and update the `rename()` block.

---

## Group or mouse IDs are `NA`

The filename does not match the regular expression used by the script.

Inspect:

```r
basename(fcs_files)
sample_meta
```

and modify the metadata parsing code.

---

## A sample has fewer cells than the downsampling target

Either:

- reduce the target number of cells per sample; or
- use the smallest sample to define the common target.

For example:

```r
downsample_n <- min(
  5000,
  min(cell_counts$n_cells)
)
```

---

## `FlowSOM` is not installed

Install it through Bioconductor:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("FlowSOM")
```

---

## The UMAP does not look identical to a previous run

Confirm that all of the following are unchanged:

- input FCS files;
- sample filtering;
- marker list;
- arcsinh cofactor;
- downsampling;
- random seed;
- `uwot` version;
- UMAP parameters.

Small numerical or software-version differences can also alter the visual orientation of a UMAP.

---

## Cluster merges or downstream objects produce errors

Because the R Markdown file contains iterative clustering and final manuscript-specific merge steps, stale objects from an earlier run can cause confusing errors.

Restart R and rerun the relevant experiment section from the beginning.

---

# Analysis philosophy

This workflow separates **unsupervised clustering** from **biological annotation**.

UMAP and FlowSOM are used to identify structure in high-dimensional cytometry data, but final interpretation is based on the combined evidence from:

- transformed marker expression;
- feature UMAPs;
- FlowSOM cluster medians;
- heatmaps;
- gated parent population;
- biological context.

Manual cluster merging is therefore a downstream phenotype-based decision and should be documented explicitly when the workflow is adapted to a new dataset.

---

# Citation

If you use or adapt this workflow in a publication, please cite the associated repository and manuscript.

Repository/manuscript citation information can be added here when the final DOI or publication information is available.

---

# Questions and issues

If you encounter a reproducible problem with the code, open an issue in the repository and include:

- the analysis section being run;
- the relevant error message;
- your R version;
- package versions;
- the output of `colnames()` for the imported FCS data;
- an example filename with any identifying information removed.
