\name{jmosaicsPattern}
\alias{jmosaicsPattern}

\title{
Call peaks and obtain combinatorial enrichment patterns
}
\description{
Call peaks and obtain combinatorial enrichment patterns
}
\usage{
jmosaicsPattern(fit_all,region_length,FDR,thres=NULL,type=c('B','E','Pattern'),patternInfo='FALSE')
}
\arguments{
  \item{fit_all}{ A list of fitted MOSAiCS models obtained using function mosaicsFit. }
\item{region_length}{ Region length, the number of bins covered by each region. }
 \item{FDR}{ False discovery rate. } 
 \item{thres}{ A peak within initial peak is removed if its ChIP tag counts are less than thres. thres could be NULL or a numeric vector, corresponding to each dataset.
   If NULL, no threshold for average ChIP tag counts for all data sets would be used. Default is NULL. } 
 \item{type}{ a vector of characters to decide which objects would be returned. Possible values are 'E': return E_LAYER result; 'B': return B_LAYER result and 'Pattern': return Pattern (enrichment pattern).}
 \item{patternInfo}{Parameter for returning information on regions. Possible values are  TRUE (return information on regions) or FALSE (do not return information on regions). Default is FALSE. }


}
\details{
The function returns objects based on `type'. 'B_LAYER' object is a list of regions which are enriched
in at least one dataset. Peak information can be accessed by `chrID', `PeakStart', `PeakStop',
`Postprob'(P(B=0jdata information)), `aveChipCount E_*', `aveInputCount E_*'. * indicates the
index of the datasets, for example: aveChipCount E_1 is the average tagCount for the first dataset. 
Each list of 'E_LAYER' object reports enriched regions for each dataset which
can be accessed by `chrID', `PeakStart', `PeakStop', `Postprob'(P(E=0j data sets)), `aveChipCount',
`maxChipCount', `aveInputCount', '`aveInputCountScaled', `aveLog2Ratio'. If region length is 1, it
can be accessed by `chrID', `PeakStart', `PeakStop', `Postprob', `ChipCount', `InputCount', `Input-
CountScaled', `Log2Ratio'. In the object of Pattern, it reports the enrichment patterns to the regions
which cover the whole genome. When the region covering more than one bin, to get the average ChIP
and input tagCount would be time consuming, the argument of 'patternInfo' let users decide whether to
report the average ChIP and input tagCount.
}
\value{
A list with following components:
\item{E_LAYER}{lists of enriched E regions, each list includes enriched regions for each data set.}
\item{B_LAYER}{list of enriched B regions.}
\item{Pattern}{list of regions annotated with patterns.}
}

\author{
Xin Zeng
}
\examples{
data("jmosaics_example_data")
bin <- readBinsMultiple(origin_bin)
fit1 <- mosaicsFit(bin[[1]], analysisType = "IO")
fit2 <- mosaicsFit(bin[[2]], analysisType = "IO")
fit <- list(fit1,fit2)
result <- jmosaicsPattern(fit, region_length=1, FDR=0.01, thres=c(10,10), type=c('B','E','Pattern'), patternInfo='FALSE')

}

\keyword{ models }
\keyword{ methods }
