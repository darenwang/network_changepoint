rm(list = ls())
#setwd()
source("functions.R")
source("BinarySegmentation.R")
source("localrefine.R")
# n is the length of between each two change points 
# p is the dimension
#B is the number of blocks
#K is the number of change point+1
#rho is the sparsity
#dist.mat.com is a matrix with 4 rows and RR number of columns




set.seed(0)
#p dimension of network
p=60
#number of blocks
B=3
#number of change point
K=3
# sparsity of the network; used in generating data
rho=0.05
#number of repetition
RR=100


#error matrices
dist.mat.com=list()
count.mat=matrix(0,nrow=4,ncol=RR)

#candidate for length of time series
N.candidate= c(80,100,120,140,160)*K

#start simulation
for ( cc in 1:length(N.candidate) ){

    dist.mat.com[[cc]]= matrix(0,nrow=4,ncol=RR)

#N is the length of time series
  N =N.candidate[cc]

# spacing parameter, only used in computing errors
n=N/K



 for(rr in 1:RR ){
   #candidata vector for SBM
   print(rr)
can.vec=sample(1:p, replace = F)

# connetivity matrix are stored in the list
con.mat.list=list()

data=matrix(nrow=p^2, ncol=0)
# generate data at each of the intervals
# The data generated are lower diagonal as the matrices are symmetry
# need some control of sbm.mean i.e. set it to be deterministic

A1 = matrix(c(0.2,1,0.2,1,0.2,0.2,0.2,0.2,0.2),nrow=3)
A2 = matrix(c(0.2,0.2,1,0.2,0.2,0.2,1,0.2,0.2),nrow=3)

for( i in 1:K){
  tem.mat=A1
  if ( i %%2 ==0){
  tem.mat = A2}
    
 
  con.mat.list[[i]]=rho* tem.mat 
  sbm.mean=generate.SBM.mean(con.mat.list[[i]] , can.vec,B,p)
  sbm.mean[upper.tri(sbm.mean,diag = T)]=0
  data=cbind(data,generate.sec.net( as.vector(sbm.mean), n,p ))
  
}
#END OF GENERATEING data
#vectorize data: only take the lower triangular part of the matrix to reduce computation complexity
data.vec=data[gen.lower.coordinate(p),]
 #data splitting
data1.vec=data.vec[,seq(2,N,2)]
data2.vec=data.vec[,seq(1,N-1,2)]
 
#rho.hat is an estimation of sparsity
rho.hat=4*mean(rowMeans(data.vec))

 

population.change=c(1, (1:K)*n)
#start with NBS
nbs=Binary.Segmentation(data1.vec,data2.vec,tau=p*rho.hat)
nbs =c(1,nbs ,N)
print("nbs = ")  ;print(nbs)
#end of NBS

#compute errors for NBS
dist.mat.com[[cc]][1,rr]=ifelse(length(nbs)==0, n,
hausdorff.distance(nbs,population.change))
 
#start with LR
#
 lr= local.refine (data.vec ,tau=sqrt(B*p*rho.hat) ,nbs,p,rho.hat,n)
 lr=c(1,lr,K*n)
 #end of LR
 #compute errors for LR
 dist.mat.com[[cc]][2,rr]=ifelse(length(lr)==0, n,
                           hausdorff.distance(lr ,population.change))
 print("lr = ")  ;print(lr)
 

 print(rowMeans(dist.mat.com[[cc]]  )/(K*n) *RR/rr)
 
 }
}
 

