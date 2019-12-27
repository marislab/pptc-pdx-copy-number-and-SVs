
if (!require("dplyr")){
  install.packages("dplyr", repos='http://cran.us.r-project.org', dependencies = TRUE)
  library(dplyr)
}

if (!require("BSgenome.Hsapiens.UCSC.hg19")){
  install.packages("https://bioconductor.org/packages/release/data/annotation/src/contrib/BSgenome.Hsapiens.UCSC.hg19_1.4.0.tar.gz", repo=NULL, type="source", dependencies = TRUE)
  library(BSgenome.Hsapiens.UCSC.hg19)
}

if (!require("rmatio")){
  install.packages("rmatio", repos='http://cran.us.r-project.org', dependencies = TRUE)
  library(rmatio)
}