# pMCMC-SIR

This is an R version of the particle MCMC code used in Rasmussen, D. A., Ratmann, O., & Koelle, K. (2011). [Inference for Nonlinear Epidemiological Models Using Genealogies and Time Series](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002136). PLoS computational biology, 7(8), e1002136. The raw code was written in MATLAB and is available from website of [Koelle research group](http://www.biology.duke.edu/koellelab). I translated the raw MATLAB code into R and the file structures of raw code and mock data are preserved.

To run the particle MCMC algorithm for SIR epidemiological model with mock data, type following command in R terminal

```r
source('main_Inference.R')
```

