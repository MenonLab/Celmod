# Celmod

Celmod is an R package designed to estimate cell type proportions in a bulk data set, based on training using matched bulk and single-cell data. It averages feature-by-feature regression estimates to generate final prediction values that reflect the strongest associations in the bulk data. 
The requirement for matched bulk and single-cell data is stringent, so it is not as broadly applicable as general deconvolution methods. In the absence of such matched data, it can be trained on "pseudo-bulk" data created from aggregated single-cell data, but these pseudo-bulk distributions may differ from those in real bulk samples.
The main advantage of Celmod is that it can "learn" cell type proportion relationships with genes from a relatively small training set, and does not require the bulk data to be in any specific format (e.g. it works fine with negative values or residuals).

To install Celmod, use the devtools package:
```
require(devtools)
devtools::install_github("MenonLab/Celmod")
```

To browse a basic vignette outlining the major commands (train_model and predict_estimates):
```
browseVignettes(Celmod)
```
