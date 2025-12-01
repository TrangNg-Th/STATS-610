# STATS-610 Final Statistical Project

## Project Overview

This repository contains the final statistical project for STATS 610. The project replicates key results from a published paper using R for statistical analysis and LaTeX for the final report.

## Project Structure

```
STATS-610/
├── R/
│   └── analysis.R          # Main R script for statistical analysis
├── data/
│   └── .gitkeep            # Directory for paper's supplementary data
├── figures/
│   └── .gitkeep            # Generated figures for the report
├── report/
│   └── final_report.tex    # LaTeX report (3-5 pages)
├── .gitignore
└── README.md
```

## Getting Started

### Prerequisites

- **R** (version 4.0 or higher)
- Required R packages:
  - tidyverse
  - ggplot2
- **LaTeX** distribution (e.g., TeX Live, MiKTeX)

### Running the Analysis

1. Download the paper's supplementary data and place it in the `data/` directory
2. Open `R/analysis.R` and update the data loading section
3. Run the analysis:
   ```bash
   cd R
   Rscript analysis.R
   ```
4. Generated figures will be saved to `figures/`

### Compiling the Report

```bash
cd report
pdflatex final_report.tex
```

For bibliography support:
```bash
pdflatex final_report.tex
bibtex final_report
pdflatex final_report.tex
pdflatex final_report.tex
```

## Report Requirements

- 3-5 pages
- Includes replication of paper's statistical results
- Comparison with original findings

## Author

[Student Name]

## Course

STATS 610 - [Course Title]