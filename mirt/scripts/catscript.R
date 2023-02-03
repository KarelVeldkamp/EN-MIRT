library(plyr)
library(mirt)

# predict probabilities based on mirt parameters
predict.mirt <- function(a, d, theta)
{ 
  prob <- t(exp(a %*% t(theta) +d) / (1 + exp(a %*% t(theta)+d)))
  
  return(prob)
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

iteminfomat <- function(theta, params){
  "
  calculates items informatin matrix given estimated theta level and item parameters
  "
  dim = length(theta) # latent dimension
  a = params[1:dim] # slopes
  d = params[dim+1] # intercept
  
  # get probability of a correct response
  P = predict.mirt(a, d, theta) 
  Q = 1-P
  
  # calculate infomration matrix
  info_mat = matrix(NA, dim, dim)
  for (r in 1:dim){
    for (s in 1:dim){
      info_mat[r,s] = a[r]*a[s]*Q*P
    }
  }
  return(info_mat)
}



simulate_cat = function(model, n_pers, n_item, design_matrix){
  "
  Simulate computerized adaptive testing given an estimated MIRT model
  model: mirt model object
  n_pers: number of persons to simulate 
  design matrix: matrix where the first column indicates the selection criterium and the second column speifies the 
                 exploration parameter
  "
  # store test attributes
  #n_item = length(model@Model$itemtype)
  n_dim = model@Model$nfact
  params = coef(model, simplify=T)$items[, 1:(n_dim+1)]
  
  # simulate random theta values and calculate the true probabilities for each person/item
  theta_true = matrix(rnorm(n_dim*n_pers), ncol=n_dim)
  
  # sample responses for each participant on each item
  true_probabilities = predict.mirt(a=params[,1:n_dim], d=params[,n_dim+1], theta=theta_true)
  responses = apply(true_probabilities, c(1,2), function(prob){rbinom(1,1,prob)})
  
  # create variables to store the average error of estimation of theta and the probability of liking an item
  avg_rmses_theta = avg_probabilities = avg_ndcg = matrix(NA, nrow=nrow(design_matrix), ncol=n_item)
  
  # loop though rows of the design matrix
  for (c in 1:nrow(design_matrix)){
    # initialize matrices for storing theta estimation errors and reponse probabilities
    error_mat_theta = matrix(NA, nrow=n_pers, ncol=n_item)
    NDCG_mat = matrix(NA, nrow=n_pers, ncol=n_item)
    
    # loop though participants to simulate CAT
    for (p in 1:n_pers){
      # print progress
      cat('\r', round(100*((c-1)*n_pers+p)/(n_pers*nrow(design_matrix)),2), '%      ')
      
      # initialize theta estimates and keep track of which items are answered
      theta_est = rep(0, n_dim)
      theta_estimates = matrix(NA, nrow=n_item, ncol=n_dim)
      answered = rep(FALSE, nrow(params))
     
      a = design_matrix[c, 'exploration'] 
      # loop though the items
      for (i in 1:n_item){
	a = a -.1
        a = max(c(0,a))
        # choose next item
        iteminfos = apply(params, 1, function(x){iteminfomat(t(theta_est), x)}, simplify=F)
        order = select_next_item(params, answered, t(theta_est), 
                                        iteminfos, exploration=a)
        next_item_ix = order[1]
        
        ordered_true_probabilities = true_probabilities[p, order]
        NDCG_mat[p, i] = NDCG(ordered_true_probabilities)
       
        # update response probabilities and response pattern
        answered[next_item_ix] = TRUE
        resp_patt = responses[p,]
        resp_patt[!answered] = NA
        
        # reestimate theta and store
        theta_est = fscores(model, response.pattern = resp_patt, method = 'MAP')[1:n_dim]
        theta_estimates[i,] = theta_est
        
      }
      # predict response probabilities given estimated theta
      preds = predict.mirt(params[, 1:n_dim], params[,n_dim+1], theta_estimates)
      # calculate RMSE of theta estimates and store
      rmses_theta = apply(theta_estimates, 1, function(est){RMSE(theta_true[p,], est)})
      error_mat_theta[p,] = rmses_theta
    }
    # average probabilities and errors over participants
    avg_probabilities[c, ] = colMeans(prob_mat)
    avg_rmses_theta[c, ] = colMeans(error_mat_theta)
    avg_dcg[c, ] = colMeans(DCG_mat)
    avg_ndcg[c, ] = colMeans(NDCG_mat)
  }
  
  # return results
  return(list(avg_rmses_theta, avg_probabilities, avg_dcg, avg_ndcg))
}

RMSE = function(true, est){
  "Helper function that calculates the RMSE"
  sqrt(mean((true-est)^2))
}

select_next_item <- function(params, answered, theta_est, iteminfos, exploration){
  "
  Function that samples next item according to criterium
  fit: mirt object
  which_not_answered: boolean vector representing whether each item was answered
  criterium: selection criterium: random, probability or determinant, mixed
  theta_est: current theta estimate for subject
  iteminfos: list of item information matrices 
  a: exploration parameter
  
  returns: index of the next item
  "
  n_dim = ncol(params)-1

  testinfo = Reduce('+', iteminfos[answered])
  if (length(testinfo) < 1){ testinfo = matrix(0, nrow=n_dim, ncol=n_dim)}
  
  # get standardised determinants for unanswered items
  dets = unlist(lapply(iteminfos, function(x){det(testinfo+x)}))
  dets = (dets - mean(dets)) / sd(dets)
  dets[answered] = -Inf
  
  # get standardised predictions for unanswered items
  pred_unanswered = predict.mirt(params[, 1:n_dim], params[,n_dim+1], theta_est)
  pred_unanswered = (pred_unanswered - mean(pred_unanswered)) / sd(pred_unanswered)
  pred_unanswered[answered] = -Inf
  crit = exploration * dets + (1-exploration) * pred_unanswered
  next_item_ix = which.max(crit)
  
  order = order(crit, decreasing=T)

  return(order)
}

args = commandArgs(trailingOnly=TRUE)
exploration = as.numeric(args[1])
n = as.numeric(args[2])

# load estimates movielens paramters
load('movielens_parameters.RData')

design_matrix = data.frame('exploration'= exploration)
results = simulate_cat(model, 1, 100,  design_matrix=design_matrix)

write(results[[1]], file=paste0('~/catresults/cat_results_theta_', exploration, '_', n, '.txt'))
write(results[[4]], file=paste0('~/catresults/cat_results_ndcg_', exploration, '_', n, '.txt'))
