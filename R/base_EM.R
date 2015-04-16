
.product=function(start,end,fun){

c=1
if(start<=end){
  c=prod(fun[start:end])
# for(i in start:end){
#   c=c*fun[i]
# }
}
return(c)
}

.f0=function(data,start,end,f_background,f_bound){
c=.product(start,end,f_background[[data]])
return(c)
}
.f1=function(data,start,end,f_background,f_bound){
c=.product(start,end,f_bound[[data]])
return(c)
}


.mosaics_post=function(object){


 mosaicsEst <- object@mosaicsEst
        tagCount <- object@tagCount
        analysisType <- object@mosaicsEst@analysisType
        
        # error treatment: invalid signal model specification
        
        
                
                pi0 <- mosaicsEst@pi0
                p1 <- mosaicsEst@p1
 fitMD <- mosaics:::.margDist_2S( mosaicsEst=mosaicsEst, tagCount=tagCount,pNfit=mosaicsEst@pNfit,k=3 )
return(fitMD)
}
.region_level_postprob=function(I, dataS,regionL,f_background,f_bound){

p1=p0=matrix(0,I,dataS)

for(i in 1:I)
for(d in 1:dataS){
  p0[i,d]=.f0(d,regionL[i]*(i-1)+1,regionL[i]*i,f_background,f_bound)
  for(v in 1:regionL[i])
  for(s in 1:(regionL[i]-v+1)){
    p1[i,d]=p1[i,d]+1/regionL[i]*1/(regionL[i]-v+1)*.f0(d,regionL[i]*(i-1)+1,regionL[i]*(i-1)+s-1,f_background,f_bound)*.f0(d,regionL[i]*(i-1)+s+v,regionL[i]*i,f_background,f_bound)*.f1(d,regionL[i]*(i-1)+s,regionL[i]*(i-1)+s+v-1,f_background,f_bound)
  }
}
p=list(p0=p0,p1=p1)
return(p)
}
.EM=function(I, dataS,p0,p1){

iter=1
q=rep(0.5,dataS)
pi=0.02
eps <- 1e-6
par_q=rep(0,dataS)
par_pi=0

L=c(q-par_q,pi-par_pi)
while( max(abs(L))>eps & iter < 100 ){
par_q=q
par_pi=pi
iter <- iter + 1
a=rep(0,I)
s1=s2=rep(1,I)

for(d in 1:dataS) {
s1=s1*(q[d]*p1[1:I,d]+(1-q[d])*p0[1:I,d])
s2=s2*p0[1:I,d] }
a=pi*s1/(pi*s1+(1-pi)*s2)

b=matrix(0,I,dataS)

for(d in 1:dataS){
b[,d]=q[d]*p1[,d]*a/(q[d]*p1[,d]+(1-q[d])*p0[,d])
}

pi=sum(a)/I
for(d in 1:dataS)
q[d]=sum(b[,d])/sum(a)
L=c(q-par_q,pi-par_pi)
#print(iter)
}


p_B=rep(0,I)
s1=s2=rep(1,I)
for(d in 1:dataS) {
s1=s1*(q[d]*p1[1:I,d]+(1-q[d])*p0[1:I,d])
s2=s2*p0[1:I,d] }
p_B=pi*s1/(pi*s1+(1-pi)*s2)

p_E=matrix(0,I,dataS)
for(d in 1:dataS)
 p_E[,d]=p_B*q[d]*p1[,d]*a/(q[d]*p1[,d]+(1-q[d])*p0[,d])


postprob_result=list(p_B=p_B,p_E=p_E,pi=pi,q=q,iter =iter )

return(postprob_result)
}