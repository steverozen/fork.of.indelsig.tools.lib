# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

`indelsig.tools.lib` is an R package for analyzing indel (insertion/deletion) mutation signatures from genomic sequencing data. The package segments indels, classifies them into different types, generates mutation catalogues (89-channel or 476-channel), and provides visualization functions.

## Installation and Development

Install the package using devtools:
```r
devtools::install()
```

Or install directly from GitHub:
```r
devtools::install_github("Nik-Zainal-Group/indelsig.tools.lib")
```

For MacOS users encountering gfortran issues during NNLM installation, install gfortran from: https://github.com/R-macos/gcc-12-branch/releases

## Core Architecture

### Classification Pipeline

The package follows a multi-stage pipeline for indel analysis:

1. **Indel Preparation** (`R/auxiliary_functions.R`): Extracts flanking sequences from reference genome and determines indel type (I/D/DI)
2. **Indel Segmentation** (`src/code.cpp`): C++ implementation for identifying repeat units and microhomology using longest tandem repeat search (LTRS)
3. **Channel Assignment** (`R/indel_classification.R` and `R/indel_classification_full.R`): Assigns each indel to specific channels based on type, length, repeat content, and flanking bases
4. **Catalogue Generation** (`R/gen_indelcatalogue.R`): Aggregates classified indels into 89 or 476-channel catalogues
5. **Visualization** (`R/plot_func.R`, `R/plot_func_full.R`): Generates publication-ready plots

### Two Classification Schemes

- **89-channel** (`indel_classifier89`): Condensed classification for standard analysis
- **476-channel** (`indel_classifier_full`): Comprehensive classification with finer granularity

Both use `assign_channels_m5` and `assign_channels_mf` respectively to categorize indels into:
- Single base indels: `[+C]`, `[+T]`, `[-C]`, `[-T]` with repeat context
- Multi-base indels: Insertions/deletions with repeat units, microhomology, spacers
- Complex indels: Simultaneous deletions and insertions (DI)

### Critical C++ Component

`src/code.cpp` contains the performance-critical segmentation logic written in Rcpp:
- `ltrs()`: Longest tandem repeat search to identify repeat units
- `segment()`: Main segmentation function that analyzes flanking sequences for repeats and microhomology

Modifications to segmentation logic require rebuilding the package with `devtools::install()`.

### Data Files

`data/` contains pre-computed templates:
- `indel_template_type_4.rda`: 89-channel template
- `indel_template_type_4_full.rda`: 476-channel template
- `*_figurelabel.rda`: Label definitions for plotting

These templates ensure consistent channel ordering in catalogues.

## Common Development Commands

### Build and Check
```r
# Build package
devtools::build()

# Run R CMD check
devtools::check()

# Load package for development
devtools::load_all()
```

### Running Tests
```r
# Run all tests (if testthat is configured)
devtools::test()
```

### Generate Documentation
```r
# Update documentation from roxygen2 comments
devtools::document()
```

### Rebuild C++ Code
After modifying `src/code.cpp`:
```r
devtools::clean_dll()  # Clean old compiled code
devtools::install()     # Rebuild and install
```

## Important Implementation Notes

### Genome Support
The package supports multiple reference genomes: hg19, hg38, mm10, canFam3. Genome selection affects:
- Reference sequences for flanking region extraction
- Chromosome naming conventions

### Verbose Parameter
Recent changes added a `verbose` parameter to `indel_classifier89()` to control output during classification. When `verbose=FALSE`, progress messages are suppressed (see commit d0bdb9b).

### Filtering High-Repeat Regions
`indel_highspecific()` filters out:
- Indels in highly repeated regions (≥10 repeats)
- Indels longer than 100 bp

This is recommended after classification to remove low-confidence mutations.

### Channel Naming Convention
Channel names follow strict patterns defined in `InDel_sigs_order`:
- Format: `[Type(Base):Context]Flanking`
- Example: `A[Ins(T):R(5,7)]G` = T insertion with 5-7 repeats flanked by A and G

When modifying channel assignments, ensure consistency with plotting templates.

## Package Dependencies

Key dependencies:
- **Bioconductor packages**: BSgenome.Hsapiens.UCSC.hg38/hg19, VariantAnnotation, GenomicFeatures
- **Rcpp**: C++ integration for performance
- **ggplot2, gridExtra**: Visualization
- **reshape2, plyr**: Data transformation
