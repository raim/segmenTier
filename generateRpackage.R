
pkg <- "~/programs/segmenTier/pkg"

## generate documentation build/check package
Rcpp::compileAttributes(pkg)
devtools::document(pkg) # 
devtools::check(pkg,cran=TRUE) # check if all is OK
devtools::build(pkg) # generate .tar.gz

library(goodpractice)
g <- gp(pkg)


## generate vignette
pkg <- "~/programs/segmenTier/pkg"
rmarkdown::render(file.path(pkg,"vignettes","segmenTier.Rmd"))
# generate doc pdf in folder 
## DOESNT WORK: check(pkg=pkg, path=pkg, args=c('--no-examples'), manual=TRUE)

pkg <- "~/programs/segmenTier/pkg"
system(paste("rm pkg.pdf -f; R CMD Rd2pdf", pkg))
