# Read command line arguments
args = commandArgs(trailingOnly=TRUE)
lambda1 = as.numeric(args[1]) # l1 penalty
lambda2 = as.numeric(args[2]) # l2 penalty
ndim = as.numeric(args[3])    # number of latent dimensions
fold1 = as.numeric(args[4])   # cross falidation fold over users
fold2 = as.numeric(args[5])   # cross valdiation fold over items

# load adapted version of mirt library (rgmirt)
library(Matrix)
library(recosystem)  

RMSE <- function(par, est)
{
  dif = par-est
  dif = dif[!is.na(dif)]
  sqrt(sum((dif)^2) / length(dif))
}

bias <- function(par, est)
{
  dif = par-est
  dif = dif[!is.na(dif)]
  mean(na.omit(dif))
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
    if (i %% 100==0) print(i)
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
  topn_accuracies[is.na(topn_accuracies)] = 0
  topn_f1_scores[is.na(topn_f1_scores)] = 0
  avg_topn_accuracy <- mean(topn_accuracies, na.rm = TRUE)
  avg_topn_f1_score <- mean(topn_f1_scores, na.rm = TRUE)
  avg_ndcg <- mean(ndcgs, na.rm = TRUE)
  
  return(c(avg_topn_accuracy, avg_topn_f1_score, avg_ndcg))
}

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
  new[,-cs] <- m
  return(new)
}

cv.nmf <- function(data,          # matrix of responses
                   ndim,         # model to be passed to mirt, can be number of dimensions
                   k = c(3,2),    # number of folds for 2d cross validation
                   lambda1,        # l1 and l2 penalty
                   lambda2,
                   fold1,
                   fold2)    # true parameter values         
{
  
  # calculate to whcih fold each person and item belongs
  folds.persons <-  cut(seq(1,nrow(data)),breaks=k[1],labels=FALSE)
  folds.items <- cut(seq(1,ncol(data)),breaks=k[2],labels=FALSE)
  
  rmses = biases = acc10 = acc20 = f110 = f120 = ndcg = c()

      
      # save the persons and items that belong to the hold out test set
      per.ind <- which(folds.persons==fold1)
      item.ind <- which(folds.items==fold2)
      
      true <- train <- data + 1
      #train[train==0] = -1
      #true[true==0] = -1
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
      # convert data to right format for recosystem package 
      data_sparse = train
      data_sparse[is.na(data_sparse)] <- 0
      sparse_matrix = as(as.matrix(data_sparse), 'dgTMatrix')
      ml.trn <- data_matrix(sparse_matrix)
      
      # Train model
      defaultW <- getOption("warn")
      options(warn = -1)
      r <- Reco()
      r$train(ml.trn, opts=list(dim=ndim, verbose=F, costp_l1=lambda1, costp_l2=lambda2, costq_l1=lambda1, costq_l2=lambda2))
      options(warn = defaultW)
      fit <- r$output(out_memory(), out_memory())
      options(warn = defaultW)
      
      
      # calculate true and predicted probablities
      p.pred <- fit$P %*% t(fit$Q)
      
      # if there are rows without responses, they have no theta estimates, so we have to fill in some empty rows to make the dimensions correct
      if(sum(emptyrows)>0){
        p.pred = insertEmptyRows(p.pred, which(emptyrows))
      }
      if(sum(emptycols)>0){
        p.pred = insertEmptyCols(p.pred, which(emptycols))
      }
      rmse <- RMSE(est=p.pred[per.ind, item.ind], par=true[per.ind, item.ind])
      bias <- bias(est=p.pred[per.ind, item.ind], par=true[per.ind, item.ind])
      
      true=true-1
      true[is.na(true)] = 0  # Count NA items as irrelevant
      p.pred=p.pred-1
      top10 <- acc_f1_at_n(p_true=true[per.ind, item.ind], p_pred=p.pred[per.ind, item.ind], n=10)
      acc10 <- top10[1]
      f110 <- top10[2]
      ndcg <- top10[3]
      
    
  return(c(rmse, bias, acc10, f110, ndcg))
}


# read datate
data1 = read.csv('~/rgmirt/movielens.csv', row.names=1)

# fit model
print('fitting model...')
start = Sys.time()
metrics <- cv.nmf(data1, 
                  ndim, 
                  k = c(3,2), 
                  lambda1,
                  lambda2,
                  fold1=fold1,
                  fold2=fold2)
end = Sys.time()
time = end-start
print(time)

# save results
fileConn<-file(paste0('~/mlresults/nmf_', lambda1,'_',lambda2,'_',ndim, '_', fold1, '_', fold2, '.txt'))
writeLines(c(as.character(metrics[1]), as.character(metrics[2]), as.character(metrics[3]), as.character(metrics[4]), as.character(metrics[5])), fileConn)
close(fileConn)


