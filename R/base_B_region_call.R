.call_B=function(fit,betapH,result,thres,region_length,FDR,binsize){
  
   betapH=1-betapH
   
   betapH_s <- sort(betapH)
    sbetapH <- cumsum(betapH_s) / c(1:length(betapH))       # expected rate of false discoveries
    id <- which( sbetapH <= FDR )
     dataS=length(fit)
    cutoff <- betapH_s[max(id)]    
    bd_bin_B <- rep( 0, length(betapH) )
    bd_bin_B[ betapH<=cutoff  ] <- 1
    if(region_length==1){
     temp=NULL
     for( d in 1:dataS){
     Y=fit[[d]]@tagCount
     temp=c(temp,which(Y>=thres[d]))
     }
     temp=unique(temp)
     id_B=intersect(which(bd_bin_B==1),temp)
     region_B=NULL
     for(d in 1:dataS)
     region_B=cbind(region_B,fit[[d]]@tagCount[id_B],fit[[d]]@input[id_B])
   
     region_B=cbind(as.data.frame(fit[[1]]@chrID[id_B]),as.data.frame(fit[[1]]@coord[id_B]),as.data.frame(fit[[1]]@coord[id_B]+binsize-1),as.data.frame(betapH[id_B]),as.data.frame(region_B))
     name=NULL
     for(i in 1:dataS)
       name=c(name,paste(c('ChipCount_E_','InputCount_E_'),i,sep=''))
     colnames(region_B)=c('chrID','PeaksStart','PeakStop','Postprob',name)
     }else{
            region_B=list()
            id_B=list()
            chrID <- result[,1]
            for(d in 1:dataS){
            chrID <- result[,1]
            Y=fit[[d]]@tagCount
            X=fit[[d]]@input
            chrList <- sort(unique(chrID))
         
   
             region_B[[d]]=id_B[[d]]=vector('list',length(chrList))
            for ( chr in 1:length(chrList) ) {
                 region_B[[d]][[chr]]=NULL
                # extract data for given chromosome
                betapH_chr<-betapH[chrID==chrList[chr]]
                coord_region_chr<-result[chrID==chrList[chr],2]
                bd_bin_chr=bd_bin_B[chrID == chrList[chr] ]
                coord_chr <- fit[[d]]@coord[ fit[[d]]@chrID  == chrList[chr] ]
                Y_chr <- fit[[d]]@tagCount[ fit[[d]]@chrID == chrList[chr] ]
                X_chr<- fit[[d]]@input[ fit[[d]]@chrID == chrList[chr] ]
                
                
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
                                    aveInputCount <- mean(X_chr[ind_i])
                                    
                                    return( c(aveChipCount, aveInputCount ) )
                                }
                            ))
               colnames(peak_info)=c('aveChipCount',  'aveInputCount' )
              peak_info=as.data.frame(peak_info)                  
                 id=which(peak_info$aveChipCount>=thres[d])
                 region_B[[d]][[chr]]=cbind(as.data.frame(chrList[chr]),as.data.frame(peak_start),as.data.frame(peak_stop),as.data.frame(betapH_chr[bd_bin_chr==1]),as.data.frame(peak_info))
                 id_B[[d]][[chr]]=id
                 }
                                                         }
  id_B_final=list()
  for(chr in 1:length(chrList) ){
   d=1
   id_B_final[[chr]]=id_B[[d]][[chr]]
   for(d in 2:dataS){
    id_B_final[[chr]]=c(id_B_final[[chr]],id_B[[d]][[chr]])
                    }
      id_B_final[[chr]]=unique(id_B_final[[chr]])
                               }
out_B=NULL
for(chr in 1:length(chrList) ){
result_B=NULL
for(d in 1:dataS){
m=region_B[[d]][[chr]][id_B_final[[chr]],5:dim(region_B[[d]][[chr]])[2]]
m=as.matrix(m)
result_B=cbind(result_B,m)
}
result_B=cbind(region_B[[d]][[chr]][id_B_final[[chr]],1:4],result_B)
out_B=rbind(out_B,result_B)
}
name=NULL
     for(i in 1:dataS)
     name=c(name,paste(c('aveChipCount_E_','aveInputCount_E_'),i,sep=''))

colnames(out_B)=c('chrID','PeakStart','PeakStop','Postprob',name)
region_B=out_B
}
region_B=as.data.frame(region_B)
return(region_B)
}

     
    



