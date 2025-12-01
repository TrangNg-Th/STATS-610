# STATS-610 Final Statistical Project

## Project Overview

This repository contains the final statistical project for STATS 610. 
The project replicates key results from a published paper using R for statistical analysis of the paper **Simulation-based model selection for dynamical systems in systems and population biology** and LaTeX for the final report.

## Project Structure

```
├── R/
│   └── codes.R          # Main R script for statistical analysis
├── supplement_data/
│   └── bioinfor-2009-1529-xxxxxxx.pdf # Tutorial on the algorithm and the paper's supplementary data
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

1. Download the paper's supplementary data and place it in the `supplement_data/` directory
2. Open `R/codes.R` and update the data loading section
3. Run the analysis:
   ```bash
   cd R
   Rscript codes.R
   ```
4. Generated figures will be saved to `figures/`


## Report Requirements

- 3-5 pages
- Includes replication of paper's statistical results
- Comparison with original findings

## Author

Trang Nguyen

## Course

STATS 610 - Statistical Computing