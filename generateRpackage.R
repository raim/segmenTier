


## generate documentation build/check package
pkg <- "~/programs/segmenTier/pkg"
library(devtools)

document(pkg) # 
check(pkg) # check if all is OK
build(pkg) # generate .tar.gz

library(goodpractice)
g <- gp(pkg)
