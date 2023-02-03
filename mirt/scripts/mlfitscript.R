# Read command line arguments
args = commandArgs(trailingOnly=TRUE)
lambda1 = as.numeric(args[1]) # l1 penalty
lambda2 = as.numeric(args[2]) # l2 penalty
ndim = as.numeric(args[5])    # number of latent dimensions
fold1 = as.numeric(args[3])   # cross falidation fold over users
fold2 = as.numeric(args[4])   # cross valdiation fold over items

# load adapted version of mirt library (rgmirt)
library(mirt)

# predict probabilities based on mirt parameters
predict.mirt <- function(a, d, theta)
{
  prob <- t(exp(a %*% t(theta) +d) / (1 + exp(a %*% t(theta)+d)))
  
  return(prob)
}

# function that inserts NA rows at each index specified in m
RMSE <- function(true, est)
{
  dif = (est-true)[!is.na(est-true)]
  rmse = sqrt(sum(dif^2) / length(dif))
  
  return(rmse)
}

# function that inserts NA rows at each index specified in m
bias <- function(true, est)
{
  dif = (est-true)[!is.na(est-true)]
  bias = mean(dif)
  
  return(bias)
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
  print(ncol(new[,-cs]))
  new[,-cs] <- m
  return(new)
}

# Function that fits rgmirt model using cross validation and returns the average RMSE and average bias
cv.mirt <- function(data,          # matrix of responses
                    model,         # model to be passed to mirt, can be number of dimensions
                    method,        # method of esitmation, e.g. EM
                    randompars,    # mirt argument for whether to use random paramters values
                    k = c(3,2),    # number of folds for 2d cross validation
                    lambda,        # l1 and l2 penalty
                    pars,          # parameters object from mirt
                    fold1,         # cross validation folds
                    fold2)         
{
  
  # calculate to which fold each person and item belongs
  folds.persons <-  cut(seq(1,nrow(data)),breaks=k[1],labels=FALSE)
  folds.items <- cut(seq(1,ncol(data)),breaks=k[2],labels=FALSE)
  
  rmses = c()
  biases = c()
  i = fold1
  j = fold2
      
  # save the persons and items that belong to the hold out test set
  per.ind <- which(folds.persons==i)
  item.ind <- which(folds.items==j)
  
  # create train and test folds
  test <- data[per.ind, item.ind]
  train <- data
  train[per.ind,item.ind] <- NA
  
  # if some rows contain no responses after the leave out set, remove them from the training set, and adjust person indices accordingly
  emptyrows = rowSums(is.na(train)) == ncol(train)
  emptycols = apply(train, 2, function(x){length(unique(na.omit(x)))<2})
  if(sum(emptyrows)>0){
    train = train[!emptyrows, ]
  }
  # do the same for columns
  if(sum(emptycols)>0){
    train = train[,!emptycols]
  }
  
  # remove empty items from pars object
  removed_items = names(which(emptycols))
  newpars = pars[!pars$item %in% removed_items, ]
  newpars$parnum = 1:nrow(newpars)
  
  # fit model
  fit <- mirt(train, model, '2PL', method=method, GenRandomPars=randompars, lambda=lambda, pars=newpars, technical = list())
  
  # save parameters
  itempars <- coef(fit, simplify = TRUE)$items
  a <- itempars[,1:(ncol(itempars)-3)]
  d <- itempars[, ncol(itempars)-2]
  theta <- fscores(fit)
  
  # calculate the predicted probabilities
  pred = predict.mirt(a, d, theta)
  # insert empty rows in prediction for items or users that did not have any responses in the training set
  if(sum(emptyrows)>0){
    pred = insertEmptyRows(pred, which(emptyrows))
  }
  if(sum(emptycols)>0){
    pred = insertEmptyCols(pred, which(emptycols))
  }
  # take only the predictions on the test set
  test_pred = pred[per.ind, item.ind]
  
  # calculate the rmse and the bias
  rmse = RMSE(test_pred, test)
  bias = bias(test_pred, test)
  
  return(c(rmse, bias))
}

# read datate
data1 = read.csv('~/rgmirt/movielens.csv', row.names=1)

# initialize pars
pars = mirt(data1, ndim, pars='values')

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

# fit model
print('fitting model...')
start = Sys.time()
metrics <- cv.mirt(data1, ndim, 
                method = 'EM',
                randompars = F, k = c(3,2), 
                lambda = c(lambda1,lambda2),
                pars = pars,
                fold1=fold1,
                fold2=fold2)
end = Sys.time()
time = end-start
print(time)

# save results
fileConn<-file(paste0('~/mlresults/movielens', lambda1,'_',lambda2,'_',ndim, '_', fold1, '_', fold2, '.txt'))
writeLines(c(as.character(metrics[1]), as.character(metrics[2])), fileConn)
close(fileConn)


