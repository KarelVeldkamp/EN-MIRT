args = commandArgs(trailingOnly=TRUE)
lambda1 = as.numeric(args[1])
lambda2 = as.numeric(args[2])
ndim = as.numeric(args[5])
fold1 = as.numeric(args[3])
fold2 = as.numeric(args[4])

library(mirt)
predict.mirt <- function(a, d, theta)
{
  prob <- t(exp(a %*% t(theta) +d) / (1 + exp(a %*% t(theta)+d)))
  
  return(prob)
}

RMSE <- function(true, est)
{
  dif = (est-true)[!is.na(est-true)]
  rmse = sqrt(sum(dif^2) / length(dif))
  
  return(rmse)
}
bias <- function(true, est)
{
  dif = (est-true)[!is.na(est-true)]
  bias = mean(dif)
  
  return(bias)
}

cv.mirt <- function(data,          # matrix of responses
                    model,         # model to be passed to mirt, can be number of dimensions
                    method,        # method of esitmation, e.g. EM
                    randompars,    # mirt argument for whether to use random paramters values
                    k = c(3,2),    # number of folds for 2d cross validation
                    lambda,        # l1 and l2 penalty
                    pars, 
                    fold1,
                    fold2)          # parameters object from mirt
{
  
  # calculate to whcih fold each person and item belongs
  folds.persons <-  cut(seq(1,nrow(data)),breaks=k[1],labels=FALSE)
  folds.items <- cut(seq(1,ncol(data)),breaks=k[2],labels=FALSE)
  
  rmses = c()
  biases = c()
  i = fold1
  j = fold2
      
  # save the persons and items that belong to the hold out test set
  per.ind <- which(folds.persons==i)
  item.ind <- which(folds.items==j)
  
  test <- data[per.ind, item.ind]
  train <- data
  train[per.ind,item.ind] <- NA
  
  # if some rows contain no responses after the leave out set, remove them from the training set, and adjust person indices accordingly
  emptyrows = rowSums(is.na(train)) == ncol(train)
  emptycols = apply(train, 2, function(x){length(unique(na.omit(x)))<2})
  if(sum(emptyrows)>0){
    train = train[!emptyrows, ]
    per.ind = per.ind[!per.ind %in% which(emptyrows)]
  }
  if(sum(emptycols)>0){
    train = train[,!emptycols]
    item.ind = item.ind[!item.ind %in% which(emptycols)]
  }
  
  # fit model and save parameters: 
  removed_items = names(which(emptycols))
  newpars = pars[!pars$item %in% removed_items, ]
  newpars$parnum = 1:nrow(newpars)
  fit <- mirt(train, model, '2PL', method=method, GenRandomPars=randompars, lambda=lambda, pars=newpars, technical = list())
  itempars <- coef(fit, simplify = TRUE)$items
  a <- itempars[,1:(ncol(itempars)-3)]
  d <- itempars[, ncol(itempars)-2]
  theta <- fscores(fit)
  
  save(theta, a, d, file=paste('~/data_parameters/parameters/ml/', lambda1, lambda2, ndim, i, j, '.RData', sep='_'))
  
  pred = predict.mirt(a, d, theta)[per.ind, item.ind]
  
  rmse = RMSE(pred, test)
  bias = bias(pred, test)
  
  return(c(rmse, bias))
}

data1 = read.csv('~/rgmirt/movielens.csv', row.names=1)



#print('initializing pars...')
#pars1 <- mirt(data1,ndim, pars = 'values', technical = list(), GenRandomPars = F)

pars = read.csv(paste0('~/rgmirt/', ndim, 'dstart.csv'), row.names=1)
pars$prior_1 = rep(NaN, nrow(pars))
pars$prior_2 = rep(NaN, nrow(pars))

# set some loadings to zero
if(ndim==1){
  pars[pars$name == 'a1' & pars$item %in% paste0('X',"1649"), ]$value = 0
  pars[pars$name == 'a1' & pars$item %in% paste0('X', "1649"), ]$est = FALSE
} else if (ndim==2){
  pars[pars$name == 'a1' & pars$item %in%  paste0('X',c("1649", "1191", "2693", "2677", "2859")), ]$value = 0
  pars[pars$name == 'a1' & pars$item %in%  paste0('X',c("1649", "1191", "2693", "2677", "2859")), ]$est = FALSE
  pars[pars$name == 'a2' & pars$item %in%  paste0('X',c("3016", "1347", "1982", "1321", "1345")), ]$value = 0
  pars[pars$name == 'a2' & pars$item %in%  paste0('X',c("3016", "1347", "1982", "1321", "1345")), ]$est = FALSE
}else if (ndim==3){
  pars[pars$name == 'a1' & pars$item %in%  paste0('X',c("1649", "1191", "2693", "2677", "2859", "1189", "1147", "3007", "162",  "246")), ]$value = 0
  pars[pars$name == 'a1' & pars$item %in%  paste0('X',c("1649", "1191", "2693", "2677", "2859", "1189", "1147", "3007", "162",  "246")), ]$est = FALSE
  pars[pars$name == 'a2' & pars$item %in%  paste0('X',c("3016", "1347", "1982", "1321", "1345", "1333", "3499", "1997", "1258", "2710")), ]$value = 0
  pars[pars$name == 'a2' & pars$item %in%  paste0('X',c("3016", "1347", "1982", "1321", "1345", "1333", "3499", "1997", "1258", "2710")), ]$est = FALSE
  pars[pars$name == 'a3' & pars$item %in%  paste0('X',c("3481", "2706", "223",  "1394", "2683", "2918", "2599", "1136", "2791", "2997")), ]$value = 0
  pars[pars$name == 'a3' & pars$item %in%  paste0('X',c("3481", "2706", "223",  "1394", "2683", "2918", "2599", "1136", "2791", "2997")), ]$est = FALSE
}



print('fitting model...')
start = Sys.time()
metrics <- cv.mirt(data1, ndim, 
                method = 'EM',
                randompars = F, k = c(3,2), 
                lambda = c(0,0),
                pars = pars,
                fold1=fold1,
                fold2=fold2)
end = Sys.time()
time = end-start
print(time)

fileConn<-file(paste0('~/mlresults/movielens', lambda1,'_',lambda2,'_',ndim, '_', fold1, '_', fold2, '.txt'))
writeLines(c(as.character(metrics[0]), as.character(metrics[1])), fileConn)
close(fileConn)


