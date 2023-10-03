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
  
  top10 <- acc_f1_at_n(p_true=p.true[per.ind, item.ind], p_pred=p.pred[per.ind, item.ind], n=10)
  
  top20 <- acc_f1_at_n(p_true=p.true[per.ind, item.ind], p_pred=p.pred[per.ind, item.ind], n=20)
  acc10 <- top10[1]
  f110 <- top10[2]
  acc20 <- top20[1]
  f120 <- top20[2]
  ndcg <- top10[3]
  
  return(c(rmse, bias, acc10, acc20, f110, f120, ndcg))
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
writeLines(c(as.character(metrics[1]), as.character(metrics[2]), as.character(metrics[3]), as.character(metrics[4]), as.character(metrics[5]), as.character(metrics[6]), as.character(metrics[7])), fileConn)
close(fileConn)


