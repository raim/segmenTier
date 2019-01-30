
pkg <- "~/programs/segmenTier/pkg"

## generate vignette

rmarkdown::render(file.path(pkg,"vignettes","segmenTier.Rmd"))


## generate documentation build/check package
library(devtools)

Rcpp::compileAttributes(pkg)
document(pkg) # 
check(pkg,cran=TRUE) # check if all is OK
build(pkg) # generate .tar.gz

library(goodpractice)
g <- gp(pkg)


