args = commandArgs(trailingOnly=TRUE)
lambda1 = as.numeric(args[1])
lambda2 = as.numeric(args[2])
ndim = as.numeric(args[3])
iteration = as.numeric(args[4])
sparsity = as.numeric(args[5])


library(mirt)

predict.mirt <- function(a, d, theta)
{
  prob <- t(exp(a %*% t(theta) +d) / (1 + exp(a %*% t(theta)+d)))
  
  return(prob)
}

# function that inserts NA rows at each index specified in m
insertEmptyRows <- function(m, rs){
  new = matrix(NA, nrow=nrow(m)+length(rs), ncol=ncol(m))
  new[-rs,] <- m
  return(new)
}
# function that inserts NA rows at each index specified in m
insertEmptyCols <- function(m, cs){
  new = matrix(NA, nrow=nrow(m), ncol=ncol(m)+length(cs))
  new[-cs,] <- m
  return(new)
}


RMSE <- function(par, est)
{
  sqrt(sum((est-par)^2) / length(par))
}

bias <- function(par, est)
{
  mean(est-par)
}

cv.mirt <- function(data,          # matrix of responses
                    model,         # model to be passed to mirt, can be number of dimensions
                    method,        # method of esitmation, e.g. EM
                    randompars,    # mirt argument for whether to use random paramters values
                    k = c(3,2),    # number of folds for 2d cross validation
                    lambda,        # l1 and l2 penalty
                    a_true,        # true parameter values
                    d_true,        # true parameter values
                    theta_true,    # true parameter values
                    pars)          # parameters object from mirt
{
  
  # calculate to whcih fold each person and item belongs
  folds.persons <-  cut(seq(1,nrow(data)),breaks=k[1],labels=FALSE)
  folds.items <- cut(seq(1,ncol(data)),breaks=k[2],labels=FALSE)
  
  rmses = c()
  biases = c()
  for (i in 1:k[1])
  {
    for (j in 1:k[2])
    {
      print(c(i,j))
      
      # save the persons and items that belong to the hold out test set
      per.ind <- which(folds.persons==i)
      item.ind <- which(folds.items==j)
      
      test <- data[per.ind, item.ind]
      train <- data
      train[per.ind,item.ind] <- NA
      
      # if some rows contain no responses after the leave out set, remove them from the training set, and adjust person indices accordingly
      emptyrows = rowSums(is.na(train)) == ncol(train)
      emptycols = colSums(is.na(train)) == nrow(train)
      if(sum(emptyrows)>0){
        train = train[!emptyrows, ]
        per.ind = per.ind[!per.ind %in% which(emptyrows)]
      }
      if(sum(emptycols)>0){
        train = train[,!emptycols]
        item.ind = item.ind[!item.ind %in% which(emptycols)]
      }
      
      # fit model and save parameters: 
      fit <- mirt(train, model, '2PL', method=method, GenRandomPars=randompars, lambda=lambda, pars=pars, technical = list())
      itempars <- coef(fit, simplify = TRUE)$items
      a <- itempars[,1:(ncol(itempars)-3)]
      d <- itempars[, ncol(itempars)-2]
      theta <- fscores(fit)
      
      save(theta, a, d, file=paste('~/data_parameters/parameters/est/est', lambda1, lambda2, ndim, iteration, sparsity, '.RData', sep='_'))
      
      # calculate true and predicted probablities
      p.pred <- predict.mirt(a, d, theta)
      p.true <- predict.mirt(a_true, as.vector(d_true), theta_true)
      
      # if there are rows without responses, they have no theta estimates, so we have to fill in some empty rows to make the dimensions correct
      if(sum(emptyrows)>0){
        p.pred = insertEmptyRows(p.pred, which(emptyrows))
      }
      if(sum(emptycols)>0){
        p.pred = insertEmptyCols(p.pred, which(emptycols))
      }
      
      rmses <- c(rmses, RMSE(est=p.pred[per.ind, item.ind], par=p.true[per.ind, item.ind]))
      biases <- c(biases, bias(est=p.pred[per.ind, item.ind], par=p.true[per.ind, item.ind]))
      
    }
  }
  return(c(mean(rmses), mean(biases)))
}

data1 = read.csv(paste('~/data_parameters/data/', ndim, iteration, sparsity, '.csv', sep='_'), row.names=1)

load(file=paste('~/data_parameters/parameters/true/true', ndim, iteration, sparsity, '.RData', sep='_'))
      

print('initializing pars...')
pars <- mirt(data1,ndim, pars = 'values', technical = list())

# set some loadings to zero
if(ndim==1){
  pars[pars$name == 'a1' & pars$item %in% paste0('Item_', 1), ]$value = 0
  pars[pars$name == 'a1' & pars$item %in% paste0('Item_', 1), ]$est = FALSE
} else if (ndim==2){
  pars[pars$name == 'a1' & pars$item %in% paste0('Item_', 1:5), ]$value = 0
  pars[pars$name == 'a1' & pars$item %in% paste0('Item_', 1:5), ]$est = FALSE
  pars[pars$name == 'a2' & pars$item %in% paste0('Item_', 6:10), ]$value = 0
  pars[pars$name == 'a2' & pars$item %in% paste0('Item_', 6:10), ]$est = FALSE
}else if (ndim==3){
  pars[pars$name == 'a1' & pars$item %in% paste0('Item_', 1:10), ]$value = 0
  pars[pars$name == 'a1' & pars$item %in% paste0('Item_', 1:10), ]$est = FALSE
  pars[pars$name == 'a2' & pars$item %in% paste0('Item_', 11:20), ]$value = 0
  pars[pars$name == 'a2' & pars$item %in% paste0('Item_', 11:20), ]$est = FALSE
  pars[pars$name == 'a3' & pars$item %in% paste0('Item_', 21:30), ]$value = 0
  pars[pars$name == 'a3' & pars$item %in% paste0('Item_', 21:30), ]$est = FALSE
}



print('fitting model...')
start = Sys.time()
metrics <- cv.mirt(data1, ndim, 
                method = 'EM',
                randompars = F, k = c(3,2), 
                lambda = c(lambda1,lambda2), 
                a_true = a1, 
                d_true = d1, 
                theta_true = theta1,
                pars = pars)
end = Sys.time()
time = end-start
print(time)

fileConn<-file(paste0('~/results/', lambda1,'_',lambda2,'_',ndim, '_',iteration,'_', sparsity, '.txt'))
writeLines(c(as.character(metrics[1]), as.character(metrics[2])), fileConn)
close(fileConn)


