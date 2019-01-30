
pkg <- "~/programs/segmenTier/pkg"



## generate documentation build/check package
library(devtools)

Rcpp::compileAttributes(pkg)
document(pkg) # 
check(pkg,cran=TRUE) # check if all is OK
build(pkg) # generate .tar.gz

library(goodpractice)
g <- gp(pkg)


## generate vignette
rmarkdown::render(file.path(pkg,"vignettes","segmenTier.Rmd"))
# generate doc pdf in folder 
## DOESNT WORK: check(pkg=pkg, path=pkg, args=c('--no-examples'), manual=TRUE)
system(paste("rm pkg.pdf -f; R CMD Rd2pdf", pkg))
