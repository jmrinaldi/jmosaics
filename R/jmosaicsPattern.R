# fit bin_all 
# fit_all list of fitted for all datsets
#current version, if readBins read N score, region length has to be 1. without N, can be 1 or larger than 1
#thres vector of cut-off at tagCount for each region for each dataset

jmosaicsPattern=function(fit_all,region_length,FDR,thres=NULL,type=c('B','E','Pattern'),patternInfo='FALSE'){

binsize=min(diff(fit_all[[1]]@coord[fit_all[[1]]@chrID==fit_all[[1]]@chrID[1]]))
 dataS=length(fit_all)
# no thres for tagCount
if(length(thres)!=dataS){thres=rep(0,dataS)}
 f_background=f_bound=vector("list",dataS)
for(d in 1:dataS){
postprob=.mosaics_post(fit_all[[d]])
f_background[[d]]=postprob$MDZ0
f_bound[[d]]=postprob$MDZ1+postprob$MDZ2
}
N=length(f_background[[d]])
if(region_length==1){
 regionL=rep(1,N)
 p0=p1=NULL
   for(d in 1:dataS){
    p0=cbind(p0,f_background[[d]])
    p1=cbind(p1,f_bound[[d]])
                    }
 postprob=.EM(N, dataS,p0,p1)
 result=cbind(as.data.frame(fit_all[[1]]@chrID),as.data.frame(fit_all[[1]]@coord))
  colnames(result)=c('chrID','coord')
 }else{
         chrID <- fit_all[[1]]@chrID
         chrList <- sort(unique(chrID))
         result <-NULL   
         coord <- fit_all[[1]]@coord
         p0_list=p1_list=list()
         p0_list_bychr=p1_list_bychr=list()
         for(d in 1:dataS){
         p0_list[[d]] <- split( f_background[[d]], chrID )
         p1_list[[d]] <- split( f_bound[[d]], chrID )}
         coord_list <- split( coord, chrID )
         result=p0=p1=NULL
         postprob=list()
         for ( chr in 1:length(chrList) ) {
                # extract data for given chromosome
                  for( d in 1:dataS){
                  p0_list_bychr[[d]] <-p0_list[[d]][[ chrList[chr] ]]
                  p1_list_bychr[[d]] <-p1_list[[d]][[ chrList[chr] ]]}
                                 
                I=floor(length(which(chrID==chrList[chr]))/region_length) 
                regionL=rep(region_length,I)                       
                p=.region_level_postprob(I, dataS,regionL,p0_list_bychr,p1_list_bychr)
                result=rbind(result,cbind(as.data.frame(rep(chrList[chr],I)),as.data.frame(seq(0,length.out=I,by=region_length*binsize))))
                p0=rbind(p0,p$p0)
                p1=rbind(p1,p$p1)}
         
      
     
 postprob=.EM(dim(result)[1], dataS,p0,p1)
colnames(result)=c('chrID','coord')
 }

# call peaks
    betapH_B <- postprob$p_B
    bd_bin_B=bd_bin_E<-list()
out=list()
if(length(which(type=='E'))>0){
    for( d in 1:dataS){
   
    bd_bin_E[[d]]=.call_E(fit_all[[d]],postprob$p_E[,d],result,thres[d],region_length,FDR,binsize)
    }
out[[length(out)+1]]=bd_bin_E
names(out)[length(out)]='E_LAYER'}
if(length(which(type=='B'))>0)
    {bd_bin_B=.call_B(fit_all,betapH_B,result,thres,region_length,FDR,binsize)
     out[[length(out)+1]]=bd_bin_B
names(out)[length(out)]='B_LAYER'}
if(length(which(type=='Pattern'))>0){
 m=cbind(postprob$p_B,postprob$p_E)
pattern_temp=t(apply(m,1,.category_fn))
pattern=pattern_temp[,1]
for(d in 2:dataS){
pattern=paste(pattern,pattern_temp[,d],sep='')

}
pattern=as.data.frame(pattern)
colnames(pattern)='Enrichment Pattern'

name=NULL
 for(i in 1:dataS){
   name=c(name,paste(c('aveChipCount_E_','aveInputCount_E_'),i,sep=''))}
if(region_length==1){
   

   for(i in 1:dataS){
   pattern=cbind(pattern,fit_all[[i]]@tagCount)
   if(length(fit_all[[i]]@input)!=0)
   pattern=cbind(pattern,fit_all[[i]]@input)
  
   }
    colnames(pattern)=c('Enrichment Pattern',name)
   pattern=cbind(result,result[,2]+binsize*region_length,as.data.frame(pattern))
   colnames(pattern)[1:3]=c('chrID','RegionStart','RegionStop')
   out[[length(out)+1]]=pattern
   names(out)[length(out)]='Pattern'
   }else{
             if(patternInfo=='TRUE'){
                   peak_infor=.peak_infor_f(result,fit_all,region_length,binsize)
                   pattern=cbind(pattern,peak_infor)
                   colnames(pattern)=c('Enrichment Pattern',name)}
   pattern=cbind(result,result[,2]+binsize*region_length,data.frame(pattern))
   colnames(pattern)[1:3]=c('chrID','RegionStart','RegionStop')
   out[[length(out)+1]]=pattern
   names(out)[length(out)]='Pattern'
                }


}
     
return(out)
}
    

.f_temp=function(p_B,x){
 which.max(c(x,p_B-x))
 }
 .category_fn=function(post){
 p_B=post[1]
 p_E=post[2:length(post)]
 
 s1=prod(apply(as.matrix(p_E),1,function(x){if(.f_temp(p_B,x)==1){return(x)}else{return(p_B-x)}}))/p_B^(length(post)-2)
 s0=prod(p_B-p_E)/p_B^(length(post)-2)
 s2=1-p_B+s0
 if(s1>s2){return(apply(as.matrix(p_E),1,function(x){if(.f_temp(p_B,x)==1){return(1)}else{return(0)}}))}
 else{return(rep(0,length(p_E)))}
 }
 

.peak_infor_f=function(result,fit,region_length,binsize){
           
            dataS=length(fit)
            chrID <- result[,1]
            result=cbind(as.data.frame(result),1:dim(result)[1])
            region=vector('list',dataS)
            for(d in 1:dataS){
            chrID <- result[,1]
            Y=fit[[d]]@tagCount
            X=fit[[d]]@input
            chrList <- sort(unique(chrID))
            for ( chr in 1:length(chrList) ) {
                
                # extract data for given chromosome
                id_chr<-result[chrID==chrList[chr],3]
                coord_region_chr<-result[chrID==chrList[chr],2]
                coord_chr <- fit[[d]]@coord[ fit[[d]]@chrID  == chrList[chr] ]
                Y_chr <- fit[[d]]@tagCount[ fit[[d]]@chrID == chrList[chr] ]
                X_chr<- fit[[d]]@input[ fit[[d]]@chrID == chrList[chr] ]
                
                
                #initial peaks
                peak_start=as.numeric(coord_region_chr)
                
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
                 region[[d]]=rbind(region[[d]],cbind(id_chr,as.data.frame(peak_info)))
                 
                 }
                                                         }
final_result=cbind(as.data.frame(result[,1:2]),result[,2]+binsize*region_length-1)
for(i in 1:dataS){
temp=region[[i]][order(region[[i]][,1],decreasing=FALSE),2:3]
final_result=cbind(final_result,temp)
}

name=NULL
     for(i in 1:dataS)
     name=c(name,paste(c('aveChipCount_E_','aveInputCount_E_'),i,sep=''))

colnames(final_result)=c('chrID','PeakStart','PeakStop',name)
return(final_result)
}







