# Didactic code example of "Three perspectives on model building and calibration for ecological risk assessment"



## Update html

To compile the html file run in R:
```r
library(rmarkdown)
render("toxicokinetic_model_example.rmd")
```

The packages `ggplot2`, `maxLik`, and `adaptMCMC` must be installed.

## Prior distributions

The file `prior_calibration.r` contains derivations of the prior distribution
together with some comments and questions.
