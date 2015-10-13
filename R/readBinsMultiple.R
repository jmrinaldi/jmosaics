# read into bins

#bin_all: list of bin level data (chip, input)for all datasets
# # read bin-level data & process it
#bin: list of all bins readed by mosaics function of readBins, if one with mappability and gc, all should include them in bin data
readBinsMultiple<- function( dataset)
{
#dataset=vector(list,length(which(type=='chip')))
#for( d in 1:length(which(type=='chip')))
#dataset[[i]]=readBins(type=unique(type), fileName=filName, excludeChr=excludeChr, dataType='unique', rounding=100, parallel=FALSE, nCore=8)

data_n=length(dataset)
existInput= existM=existGC=bin_all=vector('list',data_n)
for(i in 1:data_n){
existInput[[i]] <- ( length(dataset[[i]]@input) > 0 )
    existM[[i]] <- ( length(dataset[[i]]@mappability) > 0 )
    existGC[[i]] <- ( length(dataset[[i]]@gcContent) > 0 )
    
    }
existM=unlist(existM)
existGC=unlist(existGC)
existInput=unlist(existInput)
if(length(which(existM=='TRUE'))>0){M_d=which(existM=='TRUE')[1]}
if(length(which(existGC=='TRUE'))>0){GC_d=which(existGC=='TRUE')[1]}
chrChIP=list()
for(i in 1:data_n){   
   chrChIP[[i]] <- unique(dataset[[i]]@chrID)}

    
 if ( length(unique(unlist(chrChIP))) > 1){
 chipByChr=inputByChr=list()
 #multiple chromosome
    for(d in 1:data_n){ 
      chipByChr[[d]] <- split( as.data.frame(cbind(dataset[[d]]@coord,dataset[[d]]@tagCount)), dataset[[d]]@chrID)
      if(!is.null(dataset[[d]]@input))
      inputByChr[[d]] <- split( as.data.frame(cbind(dataset[[d]]@coord,dataset[[d]]@input)), dataset[[d]]@chrID)
                       }
      if(length(which(existM=='TRUE'))>0) 
        mapByChr <- split( as.data.frame(cbind(dataset[[M_d]]@coord,dataset[[M_d]]@mappability)), dataset[[M_d]]@chrID)
      if(length(which(existGC=='TRUE'))>0) 
        gcByChr <- split( as.data.frame(cbind(dataset[[GC_d]]@coord,dataset[[GC_d]]@gcContent)), dataset[[GC_d]]@chrID)
        chrCommon <- Reduce( intersect,chrChIP)
                    chrCommon <- sort(chrCommon)
                    
                    # provide error & stop if there is no common chromosome
                    
                    if ( length(chrCommon) == 0 ) {
                        stop( "No chromosome is common among datasets." )
                    }

                    out<-list()
                    dataList <- list()
                    for ( i in 1:length(chrCommon) ) {
                        dataList[[i]] <- list()
                        dataList[[i]]$chip=list()
                        dataList[[i]]$input=list()
                        for( d in 1:data_n){
                           dataList[[i]]$chip[[d]] <- chipByChr[[d]][[ chrCommon[i] ]]
                           if(!is.null(inputByChr[[d]][[ chrCommon[i] ]])){
                              dataList[[i]]$input[[d]] <- inputByChr[[d]][[ chrCommon[i] ]]
                                                }
                        if(length(which(existM=='TRUE'))>0)
                          if(!is.null(mapByChr[[ chrCommon[i] ]]))
                            dataList[[i]]$mapScore <- mapByChr[[ chrCommon[i] ]]
                        if(length(which(existGC=='TRUE'))>0)
                          if(!is.null(gcByChr[[ chrCommon[i] ]]))
                            dataList[[i]]$gcScore <- gcByChr[[ chrCommon[i] ]]}}

                        out <- lapply( dataList, function(x) { .match_coord(x)})

 
chrID <-c()
                    
                    for ( i in 1:length(chrCommon) ) {
                        chrID <- c( chrID, rep( chrCommon[i], length(out[[i]]$data$chip[[1]]) ) )}

coord <-c()
                    
                    for ( i in 1:length(chrCommon) ) {
                        coord <- c( coord,out[[i]]$coord)}

bin_all=dataset
i=1
for(j in 1:data_n){
bin_all[[j]]@tagCount=out[[i]]$data$chip[[j]]
if(!is.null(out[[i]]$data$input[[j]]))
bin_all[[j]]@input=out[[i]]$data$input[[j]]
if(!is.null(out[[i]]$data$mapScore))
bin_all[[j]]@mappability=out[[i]]$data$mapScore
if(!is.null(out[[i]]$data$gcScore))
bin_all[[j]]@gcContent=out[[i]]$data$gcScore
}

for ( i in 2:length(chrCommon) )
for(j in 1:data_n){
bin_all[[j]]@tagCount=c(bin_all[[j]]@tagCount,out[[i]]$data$chip[[j]])
if(!is.null(out[[i]]$data$input[[j]]))
bin_all[[j]]@input=c(bin_all[[j]]@input,out[[i]]$data$input[[j]])
if(!is.null(out[[i]]$data$mapScore))
bin_all[[j]]@mappability=c(bin_all[[j]]@mappability,out[[i]]$data$mapScore)
if(!is.null(out[[i]]$data$gcScore))
bin_all[[j]]@gcContent=c(bin_all[[j]]@gcContent,out[[i]]$data$gcScore)
}
for(j in 1:data_n){
bin_all[[j]]@chrID=chrID
bin_all[[j]]@coord=coord}
}else{
 chip=input=list()

 for(d in 1:data_n){ 
      chip[[d]] <- as.data.frame(cbind(dataset[[d]]@coord,dataset[[d]]@tagCount))
      if(!is.null(dataset[[d]]@input))
      input[[d]] <-as.data.frame(cbind(dataset[[d]]@coord,dataset[[d]]@input))}
      if(length(which(existM=='TRUE'))>0) 
       mapScore <-as.data.frame(cbind(dataset[[M_d]]@coord,dataset[[M_d]]@mappability))
      if(length(which(existGC=='TRUE'))>0) 
        gcScore<-as.data.frame(cbind(dataset[[GC_d]]@coord,dataset[[GC_d]]@gcContent))
       dataList <- list()
       dataList<- list()
       dataList$chip=list()
       dataList$input=list()
       for( d in 1:data_n){
                 dataList$chip[[d]] <- chip[[d]]
                 if(!is.null(dataset[[d]]@input))
                      dataList$input[[d]] <- input[[d]]}
       if(length(which(existM=='TRUE'))>0)
                     dataList$mapScore <-  mapScore
       if(length(which(existGC=='TRUE'))>0)
                      dataList$gcScore <- gcScore

      out <-.match_coord( dataList)
      chrID <- rep( dataset[[1]]@chrID[1], length(out$coord) )
      bin_all=dataset
      for(j in 1:data_n){
        bin_all[[j]]@tagCount=out$data$chip[[j]]
        if(!is.null(dataset[[j]]@input))
           bin_all[[j]]@input=out$data$input[[j]]
        if(!is.null(out$data$mapScore))
           bin_all[[j]]@mappability=out$data$mapScore
        if(!is.null(out$data$gcScore))
           bin_all[[j]]@gcContent=out$data$gcScore
      }
      for(j in 1:data_n){
               bin_all[[j]]@chrID=chrID
               bin_all[[j]]@coord=out$coord}

}
return(bin_all)
}

.match_coord=function(data){
c=NULL
dataS=length(data$chip)
for(i in 1:dataS){c=c(c,nrow(data$chip[[i]]))}

minID <- which.min( c)
    dataList <- data
    dataMinID <- dataList$chip[[ minID ]]
    rm( dataList, minID )
    
    # match coordinates
    coord <- data$chip[[1]][!is.na(match(data$chip[[1]][,1],dataMinID[,1])),1]
    for(i in 1:dataS){
    data$chip[[i]]<-data$chip[[i]][!is.na(match(data$chip[[i]][,1],dataMinID[,1])),2 ]
    if(!is.null(data$input[[i]]))
      data$input[[i]]<-data$input[[i]][!is.na(match(data$input[[i]][,1],dataMinID[,1])), 2]
     }

if(!is.null(data$mapScore))
      data$mapScore<-data$mapScore[!is.na(match(data$mapScore[,1],dataMinID[,1])), 2]
if(!is.null(data$gcScore))
      data$gcScore<-data$gcScore[!is.na(match(data$gcScore[,1],dataMinID[,1])),2]
result=list(coord=coord,data=data)
}