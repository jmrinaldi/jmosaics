\name{jmosaics-package}
\alias{jmosaics-package}
\alias{jmosaics}
\docType{package}
\title{Joint analysis of multiple ChIP-seq data sets}
\description{
Joint analysis of multiple ChIP-seq data sets
}
\details{
\tabular{ll}{
Package: \tab jmosaics\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-06-24\cr
License: \tab GPL(>=2)\cr
LazyLoad: \tab yes\cr

}
\code{readBinsMultiple} match coordinates for multiple datasets,\code{jmosaicsPattern} call E_LAYER and B_LAYER peaks
}
\author{
Xin Zeng, Sunduz Keles

Maintainer: Xin Zeng <xinzeng@stat.wisc.edu>

}
\references{
jMOSAiCS: Joint Analysis of Multiple ChIP-seq Datasets
}
\examples{

## step by step not run
\dontrun{
bin1 <- readBins(type = c("chip", "M", "GC", "N","input"), 
 fileName = c("h3k27me3_chip_chr10.txt",
 "./M_chr10.txt", "./GC_chr10.txt", "./N_chr10.txt",
"h3k27me3_input_chr10.txt"))
bin2 <- readBins(type = c("chip", "M", "GC", "N","input"), 
 fileName = c("h3k4me1_chip_chr10.txt",
 "./M_chr1.txt", "./GC_chr10.txt", "./N_chr10.txt",
 "h3k4me1_input_chr10.txt"))
origin_bin <- list(bin1, bin2)
}

data("jmosaics_example_data")
bin <- readBinsMultiple(origin_bin)
fit1 <- mosaicsFit(bin[[1]], analysisType = "IO")
fit2 <- mosaicsFit(bin[[2]], analysisType = "IO")
fit <- list(fit1,fit2)
result <- jmosaicsPattern(fit, region_length=1, FDR=0.01, thres=c(10,10), type=c('B','E','Pattern'), patternInfo='FALSE')


}
