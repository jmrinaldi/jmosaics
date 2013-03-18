.call_E=function(fit,betapH,result,thres,region_length,FDR,binsize){
  betapH=1-betapH
   betapH_s <- sort(betapH)
    sbetapH <- cumsum(betapH_s) / c(1:length(betapH))       # expected rate of false discoveries
    id <- which( sbetapH <= FDR )
    cutoff <- betapH_s[max(id)]    
    bd_bin_E <- rep( 0, length(betapH) )
    bd_bin_E[ betapH<=cutoff  ] <- 1
    if(region_length==1){
     Y=fit@tagCount
     temp=which(Y>=thres)
     id_E=intersect(which(bd_bin_E==1),temp)
     nRatio <- sum(fit@tagCount)/sum(fit@input)
     region_E=cbind(as.data.frame(fit@chrID[id_E]),as.data.frame(fit@coord[id_E]),as.data.frame(fit@coord[id_E]+fit@coord[2]-fit@coord[1]-1),as.data.frame (betapH[id_E]),as.data.frame(fit@tagCount[id_E]),as.data.frame(fit@input[id_E]),as.data.frame(fit@input[id_E]*nRatio),as.data.frame(log2((fit@tagCount [id_E]+1)/(fit@input[id_E]*nRatio+1))))
colnames(region_E)=c('chrID','PeakStart','PeakStop','Postprob','ChipCount',  'InputCount', 'InputCountScaled', 'Log2Ratio' )
   
     }else{
            chrID <- result[,1]
            Y=fit@tagCount
            X=fit@input
            chrList <- sort(unique(chrID))
             
            region_E=NULL
            for ( chr in 1:length(chrList) ) {
                
                # extract data for given chromosome
                betapH_chr<-betapH[chrID==chrList[chr]]
                coord_region_chr<-result[chrID==chrList[chr],2]
                bd_bin_chr=bd_bin_E[chrID == chrList[chr] ]
                coord_chr <- fit@coord[ fit@chrID  == chrList[chr] ]
                 
                Y_chr <- fit@tagCount[ fit@chrID == chrList[chr] ]
               X_chr<- fit@input[ fit@chrID == chrList[chr] ]
                
                
                #initial peaks
                peak_start=as.numeric(coord_region_chr[bd_bin_chr==1])
                
                            # extract info
                peak_stop=peak_start+region_length*binsize-1                  
                            peak_info<- t(apply( cbind(peak_start,peak_stop), 1, 
                                function(x) {
                                    start_i <- x[1]
                                    end_i <- x[2]
                                    ind_i <- which( coord_chr>=start_i & coord_chr<=(end_i-binsize+1) )
                                    
                                    aveChipCount <- mean(Y_chr[ind_i])
                                    maxChipCount <- max(Y_chr[ind_i])
                                    aveInputCount <- mean(X_chr[ind_i])
                                    nRatio <- sum(Y) / sum(X)
                                        # genome-wide sequencing depth adjustment
                                    aveInputCountScaled <- aveInputCount * nRatio
                                    aveLog2Ratio <- mean( log2( (Y_chr[ind_i]+1) / (X_chr[ind_i]*nRatio+1) ) )
                                    return( c(aveChipCount, maxChipCount, 
                                        aveInputCount, aveInputCountScaled, aveLog2Ratio ))
                                }
                            ))
               colnames(peak_info)=c('aveChipCount', 'maxChipCount', 
                                        'aveInputCount', 'aveInputCountScaled', 'aveLog2Ratio' ) 
                peak_info=as.data.frame(peak_info)
                 
                 id=which(peak_info$aveChipCount>=thres)
                 region_E=rbind( region_E,cbind(as.data.frame(chrList[chr]),as.data.frame(peak_start[id]),as.data.frame(peak_stop[id]),as.data.frame (betapH_chr[bd_bin_chr==1][id]),peak_info[id,]))
                   }
                 colnames(region_E)=c('chrID','PeakStart','PeakStop','Postprob','aveChipCount', 'maxChipCount', 'aveInputCount', 'aveInputCountScaled',  'aveLog2Ratio' )}
region_E=as.data.frame(region_E)
return(region_E)
}

     