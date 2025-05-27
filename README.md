# UTMOST_fa

`UTMOST_fa` is an R package for joint gene expression imputation across multiple contexts, incorporating functional annotations. It is particularly tailored for single-cell transcriptome-wide association studies (sc-TWAS), extending the original [UTMOST](https://github.com/yiminghu/CTIMP/tree/master) training process with improved modeling of regulatory SNPs.

## 1. Installation

`UTMOST_fa` can be installed from this [GitHub repository](https://github.com/leaffur/UTMOST_fa).

``` r
# Install the devtools package if you haven't already
# install.packages("devtools")

devtools::install_github("leaffur/UTMOST_fa")
```

## 2. Tutorial

A detailed tutorial is under development and will be provided soon.

## 3. Association tests using external GWAS

`UTMOST_fa` focuses on building and training gene expression imputation models for TWAS analyses. To apply these trained models on external GWAS summary statistics and identify disease-associated genes, please refer to [UTMOST](https://github.com/Joker-Jerome/UTMOST) framework.

** Note for Python 3 users:
We provide a Python 3-compatible alternative for UTMOST. 

To use it:
**1. Clone the UTMOST repository**
```bash
$ git clone https://github.com/Joker-Jerome/UTMOST
```
**2. Download and extract the `UTMOST_py3.zip` file from the main `UTMOST_fa` repository folder.**
**3. Replace the `single_tissue_covariance.py` file and the `metax` folder inside the UTMOST directory with the provided Python 3 versions.**

## 4. Reference

  - *Hu, Y. et al. A statistical framework for cross-tissue transcriptome-wide association analysis. Nature Genetics 51, 568-576 (2019).* [Link](https://www.nature.com/articles/s41588-019-0345-7).
