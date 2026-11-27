# Code for McCullough et al. 2026
## Drift drives phenotypic evolution in a rapid island radiation of Kingfishers 

This repository holds the code necessary to perform all analyses for this study. 

Study from the publisher's website: https://www.cell.com/current-biology/abstract/S0960-9822(26)01443-0 

If you need a copy, PLEASE email me! If I don't respond in a couple days, please email me again! 

Website for the Dryad: 
https://datadryad.org/dataset/doi:10.5061/dryad.dv41ns2f4 

Website for the additional results (discussing time calibration and sensitivity analyses): 
https://zenodo.org/records/24322666

What is not included in this code are details for time calibration. Please see the Dryad for these files. 

There are two files that are on the Dryad that must be downloaded to follow along with the code, specifically the .rda file "specs4_metaclean4-final-extra-cleaning.rda" in the 02_color-data zipped file and Todiramphus_specific.zip. 

Once these are downloaded an in your working directory, you can follow along with this code. 

 There are six main parts: 
 
## PART 1. Calculating diversification with CLaDS and ES metrics. 
This section begins with a time calibrated tree and will end with comparisons on these metrics. You don't have to do the actual analysis (which took me several days when I originally did it myself on my laptop) and instead can load the results from diversification-output/mcmc_phy-200K-thin100.rda 
 
## PART 2. Creation of excel files from spectra that have averaged X,Y,Z coordinates.  
Here we will be projecting the spectra into tetrahedral colorspace and averaging within individuals and sexes to create the XYZ coordinates and luminance values that is used for the rest of the study. This is calculated from the specs4_metaclean4-final-extra-cleaning.rda file that is available on Dryad. This will recreate the excel files that are on Dryad and this GitHub (i.e., color-data/tax95-f-chrom-avg-specs3.xyz.nolores.csv etc.) 

## PART 3. Estimating color evolutionary rates with RRphylo. 
Using the sex-specific files for chromatic and achromatic color, we will estimate evolutionary rates of color across the phylogeny. 

## PART 4. Preparing the adaptive explanatory variables. 
We will A) calculating dichromatism, as well as testing whether colordist vs. Euclidian distances give us the same general answers; B) calculating island size metrics, and C) log transform island sizes and sympatry values.   
 
## PART 5. Testing drivers of rapid color evolution with linear and phylogenetic linear models. 
We will also plot many of the supplemental figures, including figures S7 - S14. 

## PART 6. Sensitivity tests for heterozygosity. 
We specifically will look at different ways to account for the potential bias that short branches may have when trying to understand the rate of trait changes. We address three questions: 1. Are heterozygosity estimates are uncorrelated with taxon age (i.e. terminal branch length in the tree)? 2. Do relationships persist when using colour rates derived from lambda-scaled trees? 3. Do relationships  persist when including taxon age (i.e. terminal branch length in the tree) as an additional predictor in the multipredictor models to account for the fact that short branches are invariably inferred to have faster rates? 

