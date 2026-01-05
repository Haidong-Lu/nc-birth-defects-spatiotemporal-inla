# Spatiotemporal Trends of Birth Defects in North Carolina (2003–2015)

This repository provides example R code for fitting a spatiotemporal Poisson model
using **Integrated Nested Laplace Approximation (INLA)**, as used in the paper:

> **Spatiotemporal Trends of Birth Defects in North Carolina, 2003–2015**  
> Haidong Lu, Andrew F. Olshan, Marc L. Serre, Kurtis M. Anthony,  
> Rebecca C. Fry, Nina E. Forestieri, Alexander P. Keil

The modeling framework includes:
- BYM spatial random effects
- RW1 temporal random effects
- IID temporal effects
- Space–time interaction terms

## Data availability
No data are included in this repository.

Due to data use restrictions, the original birth defect data and geographic
identifiers from North Carolina cannot be shared publicly. Users are expected
to supply their own analysis-ready area–time count data and corresponding
adjacency structure.

## Intended use
This code is provided to illustrate the **modeling strategy and implementation**
used in the paper. It is suitable for adaptation to other spatiotemporal
epidemiologic applications involving count data.

## Requirements
- R (≥ 4.0)
- INLA (`https://www.r-inla.org`)

## Disclaimer
This repository reproduces the statistical modeling approach only.
It does not reproduce the study results without access to the original data.

## Citation
If you use this code, please cite:

Lu H, Olshan AF, Serre ML, Anthony KM, Fry RC, Forestieri NE, Keil AP.  
*Spatiotemporal Trends of Birth Defects in North Carolina, 2003–2015.*
