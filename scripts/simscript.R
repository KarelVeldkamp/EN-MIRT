library(mirt)
args = commandArgs(trailingOnly=TRUE)
ndim = as.numeric(args[1])
iteration = as.numeric(args[2])
sparsity = as.numeric(args[3])

theta1 <- matrix(rnorm(10000*ndim), ncol=ndim)
d1 <- matrix(runif(400, -2, 2))
a1 <- matrix(runif(400*ndim,.5,2), ncol=ndim)

a1[1] <- 0
if(ndim==1){
  a1[1] = 0
} else if (ndim==2){
  a1[1:5, 1] = 0
  a1[6:10, 2] = 0
}else if (ndim==3){
  a1[1:10, 1] = 0
  a1[11:20, 2] = 0
  a1[21:30, 3] = 0
}

save(theta1, a1, d1, file=paste('~/data_parameters/parameters/true/true',  ndim, iteration, sparsity, '.RData', sep='_'))
      

print('simulating data...')
data1 = simdata(a1, d1, itemtype = '2PL', Theta = theta1)
data1[sample(length(data1), length(data1)*sparsity, replace = FALSE)] <- NA

write.csv(data1, file=paste('~/data_parameters/data/', ndim, iteration, sparsity, '.csv', sep='_'))
