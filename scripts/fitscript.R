# read command line arguments
args = commandArgs(trailingOnly=TRUE)
lambda1 = as.numeric(args[1])   # l1 penalty 
lambda2 = as.numeric(args[2])   # l2 penalty
ndim = as.numeric(args[3])      # number of latent dimensions
iteration = as.numeric(args[4]) # number of the current iteration
sparsity = as.numeric(args[5])  # sparsity of the dataset 

# load adapted version of mirt library (rgmirt)
library(mirt)

# predict probabilities based on mirt parameters
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

# Root mean square error
RMSE <- function(par, est)
{
  sqrt(sum((est-par)^2) / length(par))
}

# bias 
bias <- function(par, est)
{
  mean(est-par)
}

# calcualte DCG based on relevance grades (ordered)
DCG = function(r){
  k = 1:length(r)
  
  D = log(k+1, base = 2)
  DCG = sum(r/D)
  return(DCG)
}

# NDCG
NDCG = function(r){
  NDCG = DCG(r)/DCG(sort(r, decreasing=T))
  return(NDCG)
}

acc_f1_at_n <- function(p_true, p_pred, n) {
  # binarise probabilities to 0/1
  true_resp_matrix <- round(p_true, 0)
  pred_resp_matrix <- round(p_pred, 0)
  
  num_users <- dim(p_true)[1]
  num_movies <- dim(p_true)[2]
  
  topn_accuracies <- vector("numeric", num_users)
  topn_f1_scores <- vector("numeric", num_users)
  ndcgs <- vector("numeric", num_users)
  
  # for each user
  for (i in 1:num_users) {
    true_response <- true_resp_matrix[i, ]
    pred_response <- pred_resp_matrix[i, ]
    
    pred_relevance <- p_pred[i, ]
    
    # Get the indices of the top-n predicted relevance
    indices <- order(-pred_relevance)
    topn_indices <- indices[1:n]
    
    # Calculate the top-n accuracy for the current user
    topn_accuracy <- sum(true_response[topn_indices]) / n
    
    # Calculate the top-n F1 score for the current user
    true_positives <- sum(true_response[topn_indices]*pred_response[topn_indices])
    predicted_positives <- sum(pred_response[topn_indices])
    actual_positives <- sum(true_response[topn_indices])
    
    precision <- true_positives / n
    recall <- true_positives / actual_positives
    
    f1_score <- 2 * precision * recall / (precision + recall)
    
    ordered_true_probabilities = p_true[i, indices]
    ndcg = NDCG(ordered_true_probabilities)
    
    topn_accuracies[i] <- topn_accuracy
    topn_f1_scores[i] <- f1_score
    ndcgs[i] <- ndcg
  }
  
  # Calculate the average top-n accuracy and top-n F1 score across all users
  avg_topn_accuracy <- mean(topn_accuracies, na.rm = TRUE)
  avg_topn_f1_score <- mean(topn_f1_scores, na.rm = TRUE)
  avg_ndcg <- mean(ndcgs, na.rm = TRUE)
  
  return(c(avg_topn_accuracy, avg_topn_f1_score, avg_ndcg))
}


# Function that fits rgmirt model using cross validation and returns the average RMSE and average bias
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
  
  # calculate to which fold each person and item belongs
  folds.persons <-  cut(seq(1,nrow(data)),breaks=k[1],labels=FALSE)
  folds.items <- cut(seq(1,ncol(data)),breaks=k[2],labels=FALSE)

  rmses = c()
  biases = c()
  accuracies = c()
  f1s = c()
  acc10 = acc20 = f110 = f120 = ndcg = c()
  
  # loop though cross validation folds
  for (i in 1:k[1])
  {
    for (j in 1:k[2])
    {
      # save the persons and items that belong to the hold out test set
      per.ind <- which(folds.persons==i)
      item.ind <- which(folds.items==j)
      
      # create train and test folds
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
      # do the same for columns
      if(sum(emptycols)>0){
        train = train[,!emptycols]
        item.ind = item.ind[!item.ind %in% which(emptycols)]
      }
      
      # fit model
      fit <- mirt(train, model, '2PL', method=method, GenRandomPars=randompars, lambda=lambda, pars=pars, technical = list(NCYCLES=1))
      # save parameters
      itempars <- coef(fit, simplify = TRUE)$items
      a <- itempars[,1:(ncol(itempars)-3)]
      d <- itempars[, ncol(itempars)-2]
      theta <- fscores(fit)
      
      # calculate true and predicted probablities
      p.pred <- predict.mirt(a, d, theta)
      p.true <- predict.mirt(a_true, as.vector(d_true), theta_true)
      save(itempars, theta, file=paste('~/rgmirt/parameters/est/par', lambda1, lambda2, ndim, iteration, sparsity,i,j, '.csv', sep='_'))
      # if there are rows without responses, they have no theta estimates, so we have to fill in some empty rows to make the dimensions correct
      if(sum(emptyrows)>0){
        p.pred = insertEmptyRows(p.pred, which(emptyrows))
      }
      if(sum(emptycols)>0){
        p.pred = insertEmptyCols(p.pred, which(emptycols))
      }
 
      # calculate metrics and save them
      rmses <- c(rmses, RMSE(est=p.pred[per.ind, item.ind], par=p.true[per.ind, item.ind]))
      biases <- c(biases, bias(est=p.pred[per.ind, item.ind], par=p.true[per.ind, item.ind]))
      top10 <- acc_f1_at_n(p_true=p.true[per.ind, item.ind], p_pred=p.pred[per.ind, item.ind], n=10)

      top20 <- acc_f1_at_n(p_true=p.true[per.ind, item.ind], p_pred=p.pred[per.ind, item.ind], n=20)
      acc10 <- c(acc10, top10[1])
      f110 <- c(f110, top10[2])
      acc20 <- c(acc20, top20[1])
      f120 <- c(f120, top20[2])
      ndcg <- c(ndcg, top10[3])
      
    }
  }
  return(c(mean(rmses), mean(biases), mean(acc10), mean(acc20), mean(f110), mean(f120), mean(ndcg)))
}

# read appropriate data file
data1 = read.csv(paste('~/rgmirt/data/', ndim, iteration, sparsity, '.csv', sep='_'), row.names=1)

# load true parameter values
load(file=paste('~/rgmirt/parameters/true/true', ndim, iteration, sparsity, '.RData', sep='_'))
      
# initialize paramters
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

# fit model
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

# save results
fileConn<-file(paste0('~/rgmirt/results/', lambda1,'_',lambda2,'_',ndim, '_',iteration,'_', sparsity, '.txt'))
writeLines(c(as.character(metrics[1]), as.character(metrics[2]), as.character(metrics[3]), as.character(metrics[4]), as.character(metrics[5]), as.character(metrics[6]), as.character(metrics[7])), fileConn)
close(fileConn)



