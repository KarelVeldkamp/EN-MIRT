# Data
these are the results from our simulation studies and the application to the movielens dataset:

### CAT_simulation_results.csv
what: indicates which metric is given in this row (ndcg=the NDCG of the ranking, theta= the RMSE of theta)
method: Which method was determined to select the next item (Exploit or Explore-Exploit)
X1, X2 ... X300: the value of the metrix at timepoints 1 though 300. 

### MIRT_simulation_results.csv
lambda1: l1 penalty used during estimation
lambda2: l2 penalty used during estimation
ndim: number of latent dimensions
iteration: which of iterations/datasets produced this result (twenty datasets were simulated for each combination of ndim and sparsity)
sparsity: the sparsity of the dataset
rmse: root mean square error of the estimated probabilities
bias: bias of the estimated probabilities

### Movielens_results.csv
lambda1: l1 penalty used during estimation
lambda2: l2 penalty used during estimation
ndim: number of latent dimensions
fold1: which fold this metric was calculated on (fold over users)
fold2: which fold this metric was calculated on (fold over items)
rmse: root mean square error of the predicted response probabilities
bias: bias of the predicted response probabilities

