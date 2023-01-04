simulate_cat = function(fit, n_pers, design_matrix, plot=TRUE){
  "
  Simulate computerized adaptive testing given an estimated MIRT model
  fit: mirt object
  n_pers: number of persons to simulate 
  criterium: selection criterium, random, probability or determinant
  "
  # store test attributes
  n_item = length(fit@Model$itemtype)
  items = lapply(1:n_item, function(x){mirt::extract.item(fit, x)})
  n_dim = fit@Model$nfact
  
  # simulate random theta values and calculate the true probabilities for each person/item
  theta_true = matrix(rnorm(n_dim*n_pers), ncol=n_dim)
  params = coef(fit, simplify=TRUE)$items
  
  # sample responses for each participant on each item
  true_probabilities = predict.mirt(a=params[,1:ndim], d=params[,ndim+1], theta=theta_true)
  responses = apply(true_probabilities, c(1,2), function(prob){rbinom(1,1,prob)})
  
  avg_rmses_theta = avg_probabilities = matrix(NA, nrow=nrow(design_matrix), ncol=n_item)
  for (c in 1:nrow(design_matrix)){
    # initialize matrix for storing theta estimation errors
    error_mat_theta = matrix(NA, nrow=n_pers, ncol=n_item)
    prob_mat = matrix(NA, nrow=n_pers, ncol=n_item)
    # loop though participants to simulate CAT
    for (p in 1:n_pers){
      cat('\r', round(100*((c-1)*n_pers+p)/(n_pers*nrow(design_matrix)),2), '%      ')
      # initialize thea estimates and tracking variables
      theta_est = c(0,0)
      theta_estimates = matrix(NA, nrow=n_item, ncol=n_dim)
      answered = rep(FALSE, n_item)
      # loop though the number of items
      for (i in 1:n_item){
        # choose and administer next item
        next_item_ix = select_next_item(fit, answered, design_matrix[c, 'criterium'], t(theta_est), 
                                        iteminfos, a=design_matrix[c, 'a'])
        prob_mat[p, i] = true_probabilities[p, next_item_ix]
        answered[next_item_ix] = TRUE
        resp_patt = responses[p,]
        resp_patt[!answered] = NA
        
        # reestimate theta and add to 
        theta_est = fscores(fit, response.pattern = resp_patt)[1:n_dim]
        
        theta_estimates[i,] = theta_est
      }
      params = coef(fit, simplify=TRUE)$items
      preds = predict.mirt(params[, 1:n_dim], params[,ndim+1], theta_estimates)
      
      rmses_theta = apply(theta_estimates, 1, function(est){RMSE(theta_true[p,], est)})
      error_mat_theta[p,] = rmses_theta
    }
    avg_probabilities[c, ] = colMeans(prob_mat)
    avg_rmses_theta[c, ] = colMeans(error_mat_theta)
  }
  if(plot){
    plot(avg_rmses_theta[1, ], type='l', ylab ='RMSE(Theta)', xlab='Number of items', main='Theta')
    if (nrow(design_matrix)>1){
      for (row in 2:nrow(design_matrix)){
        lines(avg_rmses_theta[row, ], col=row)
      }
    }
    
    legend('topright', legend=paste(design_matrix[,1], design_matrix[,2]), fill= 1:nrow(design_matrix))
    
    plot(avg_probabilities[1, ], type='l', ylim=c(0,.8), ylab ='Probability correct', xlab='Number of items', main='Probability')
    if (nrow(design_matrix)>1){
      for (row in 2:nrow(design_matrix)){
        lines(avg_probabilities[row, ], col=row)
      }
    }
    
    legend('topright', legend=paste(design_matrix[,1], design_matrix[,2]), fill= 1:nrow(design_matrix))
  }
  
  return(list(avg_rmses_theta, avg_probabilities))
}

RMSE = function(true, est){
  sqrt(mean((true-est)^2))
}

select_next_item <- function(fit, answered, criterium, theta_est, iteminfos, a){
  "
  Function that samples next item according to criterium
  fit: mirt object
  which_not_answered: boolean vector representing whether each item was answered
  criterium: selection criterium: random, probability or determinant, mixed
  theta_est: current theta estimate for subject
  iteminfos: list of item information matrices 
  
  returns: index of the next item
  "
  if (criterium == 'random'){
    # sample random unanswered item
    next_item_ix = sample(which(!answered), 1)
  } else if(criterium == 'probability'){
    # ignore items that are already answered
    params = coef(fit, simplify=TRUE)$items
    pred_unanswered = predict.mirt(params[, 1:n_dim], params[,ndim+1], theta_est)
    pred_unanswered[answered] = NA
    # choose item with max probability
    next_item_ix = which.max(pred_unanswered)
  } else if(criterium=='determinant'){
    testinfo = Reduce('+', iteminfos[answered])
    if (length(testinfo) < 1){ testinfo = matrix(0, nrow=n_dim, ncol=n_dim)}
    # calculate hypothetical determinant of information matrix for each item to be added
    dets = lapply(iteminfos[!answered], function(x){det(testinfo+x)})
    next_item_ix = which(!answered)[which.max(dets)]
  } else if(criterium=='mixed'){
    testinfo = Reduce('+', iteminfos[answered])
    if (length(testinfo) < 1){ testinfo = matrix(0, nrow=n_dim, ncol=n_dim)}
    
    # get standardised determinants for unanswered items
    dets = unlist(lapply(iteminfos, function(x){det(testinfo+x)}))
    dets = (dets - mean(dets)) / sd(dets)
    dets[answered] = -Inf
    
    # get standardised predictions for unanswered items
    params = coef(fit, simplify=TRUE)$items
    pred_unanswered = predict.mirt(params[, 1:n_dim], params[,ndim+1], theta_est)
    pred_unanswered = (pred_unanswered - mean(pred_unanswered)) / sd(pred_unanswered)
    pred_unanswered[answered] = -Inf
    crit = a * dets + (1-a) * pred_unanswered
    next_item_ix = which.max(crit)
  } else stop('Invalid criterium')
  return(next_item_ix)
}


design_matrix = data.frame('criterium'=c('determinant', 'probability', 'mixed', 'mixed', 'mixed', 'mixed'), 'a'= c(1, 0, .2, .4, .6, .8))
results = simulate_cat(fit, 1000, design_matrix=design_matrix)

library(plyr)
theta_est_df = reshape::melt(results[[1]])
colnames(theta_est_df) = c('criterium', 'nitems', 'value')
theta_est_df$criterium = mapvalues(theta_est_df$criterium, from=1:6, to = c(1,0, .2, .4, .6, .8))


prob_df = reshape::melt(results[[2]])
colnames(prob_df) = c('criterium', 'nitems', 'value')
prob_df$criterium = mapvalues(prob_df$criterium, from=1:6, to = c(1,0, .2, .4, .6, .8))



ggplot(theta_est_df, aes(x=nitems, y=value)) +
  geom_line(aes(col =factor(criterium)), size=.7) +
  scale_color_viridis(discrete = T, option='G', name='exploration parameter') +
  xlab("Number of items") + 
  ylab("RMSE(Theta)") +
  ggtitle('Average error of theta estimation after n items')


ggplot(prob_df, aes(x=nitems, y=value)) +
  geom_smooth(aes(color =factor(criterium))) +
  scale_color_viridis(discrete = T, option='G', name='exploration parameter') + 
  ggtitle('Average probability correct after n items')
