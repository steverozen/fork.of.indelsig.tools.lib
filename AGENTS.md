# Repository Guidelines

Use this guide to keep indelsig.tools.lib segmentation, catalogue, and plotting workflows consistent.

## Project Structure & Module Organization
Core exported functions live under `R/` with roxygen2 headers feeding `man/`; rerun `devtools::document()` after any code or example update. Keep Rcpp helpers inside `src/` and rebuild the `.so` through the standard toolchain rather than editing binaries. Templates and plotting palettes belong in `data/`, with regeneration scripts inside `data-raw/` and `indelsig.tools.lib.develop.R`. Workflow utilities sit in `signature_extraction_workflow/`, demo inputs and rendered figures live in `example/`, and notebooks or replication material belong in `vignettes/` or `Replication_Strand_Bias/`.

## Build, Test, and Development Commands
- `R -q -e "devtools::document()"` – regenerate Rd files and the `NAMESPACE`.
- `R -q -e "devtools::install()"` – install locally for iterative work.
- `R CMD build .` – produce the source tarball with compiled code and data.
- `R -q -e "devtools::check()"` – run the full `R CMD check`; required before opening a PR.
- `R -q -e "source('indelsig.tools.lib.develop.R')"` – refresh internal catalogue templates prior to release tagging.

## Coding Style & Naming Conventions
Use two-space indents, 80-character lines, and whitespace around `=` in named arguments, mirroring `R/indel_classification*.R`. Prefer descriptive snake_case or verb phrases (`indel_classifier89`, `gen_catalogue89`) and keep exported names aligned with catalogue types. Maintain roxygen2 headers with runnable examples and genome-version notes, and reuse the tidyverse idioms already imported.

## Testing Guidelines
Adopt `testthat` (declared in `Suggests`). Place specs under `tests/testthat/test_<topic>.R`, mirroring function or workflow names. Use deterministic fixtures from `data/` or minimal `.vcf` snippets in `example/`. Run `R -q -e "devtools::test()"` for quick checks and `R -q -e "devtools::check(filter = 'classifier')"` when scoping to one module. Cover hg19/hg38 branches, plotting fallbacks, and add regression tests whenever catalogue templates change.

## Commit & Pull Request Guidelines
Git history (`git log --oneline`) favors concise, imperative commits such as “added verbose option to silence output during classification”; follow that tone and optionally prefix scope tags (`classifier: ...`). Each PR should summarize the biological rationale, list API or data updates, link the motivating issue, and call out regenerated `data/` or `man/` artifacts. Attach updated plot paths or screenshots when touching visualization helpers.

## Security & Data Handling
Never upload patient-identifiable catalogues or raw VCFs; stage synthetic or redacted samples in `example/`. Large intermediates (e.g., `.o`, `.so`) must remain reproducible from `src/`. Keep credentials, genome indexes, and BSgenome downloads outside the repo and reference them through documented environment variables or vignettes.
