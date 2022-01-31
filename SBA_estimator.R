# This Estimator comes from Job market paper of Xiaoying Lan
# Link https://drive.google.com/file/d/19ru2U36whsvIEbjGpWwG0CLtX5eS8Uyg/view?usp=sharing

# This function provides estimation for logistic regression via gradient method
Gradient_logistic_estimation <- function(X1, Y1, scaling = 100, gb = 1){
	beta_s <- matrix(0,nrow = dim(X1)[2], ncol = 1)
	for (jj in 1: (dim(X1)[1] / scaling)){
		X1_beta_s = pmin(X1 %*% beta_s, 700) # limit the number inserted in exp()
		beta_s <- beta_s - 1 / jj^gb / dim(X1)[1]* t(X1) %*%
		  ( exp(X1_beta_s) /(1+exp(X1_beta_s)) - Y1)
	}
	return(beta_s)
}

SBA_predict <- function(regsult, regor)
  # This function is used to predict probability of any feature using SBA_estimator
  # regsult is the SBA_estimator
  # regor is the matrix or vector of regressors
{
  pp = length(regsult[[1]][,1])
  qq = length(regsult[[2]])
  beta_pre = append(1, regsult[[1]][,1])
  index_predict = matrix(regor, ncol = pp+1) %*% matrix(beta_pre, nrow = pp +1)
  index_reg <- cbind(index_predict^0,index_predict^1) # get polynomial estimator
  for (ii in 2:(qq-1)){
    index_reg <- cbind(index_reg, index_predict^ii)
  }
  sba_pre <- pmin(index_reg %*% regsult[[2]],700)
  sba_pre <- exp(sba_pre) / (1 + exp(sba_pre))
  return(sba_pre)
}

# estimation and inference
SBA_estimator <- function(X, Y, beta_0 = NULL, q = 3, l = 100, ga = 0.9, iter_times = NULL, iter_times2 = 100, sieve_w = 10)
  # X is regressor
  # Y is binary outcome
  # beta_0 is the initial value of guess, default as 1
  # q is the number of polynomials
  # l is the learning rate gamma_1
  # ga is the learning rate gamma
  # iter_times is the iterate times of updating beta, default as the number of observation
  # itera_times2 is the iterate times of updating error function, default 200
  # sieve is the tuning parameter of magnitude of  polynomials
  # It returns the coefficients (first coefficient is eliminated and default as 1),  standard error and coeffcients of error distribution.
{ library(Matrix)
	library(expm)
	library(pracma)
  library(gsubfn)
# estimation  
  N <- dim(X)[1]
  p <- dim(X)[2]
  if (is.null(beta_0)){
		beta_0 <- rep(1, p)
		}
  if (is.null(iter_times)){
		iter_times <- N
    }
  beta_t <- beta_0 # beta_t is used for iteration
  dim(beta_t) <- c(p, 1)
  beta_ave <- matrix(0, nrow = p, ncol = 1)
  # estimation
  for (j in 1:iter_times){
  	X_trans <- X %*% beta_t
  	X_reg <- cbind(X_trans^0, X_trans^1) # get polynomial estimator
  	for (ii in 2:q){
  		X_reg <- cbind(X_reg, X_trans^ii)
  	}
#  	Iterate_fun <- function(q_search, ind_var, de_var) sum(log(1 + exp(ind_var %*% as.matrix(q_search)))) - t(de_var) %*% ind_var %*% as.matrix(q_search)
#  	Iterate_result <- fminsearch(Iterate_fun, c(1:(q+1)), ind_var = X_reg, de_var = Y)
  	coeff_q <- Gradient_logistic_estimation(X_reg, Y, scaling = iter_times2)  # get coefficient of error distribution
  	g_trans <- X_reg %*% coeff_q
  	g_trans <- pmin(g_trans, 700) # limit the number inserted in exp()
  	g_hat <- exp(g_trans)/(1+exp(g_trans))
  	beta_t <- beta_t - l / j^ga  * solve(t(X) %*% X) %*% t(X) %*% (g_hat - Y) #  BGD
  	beta_t1 <- beta_t
  	beta_ave <- ((j - 1) * beta_ave  + beta_t1) /j
  }
# inference: calculate standard error
  X_trans <- X %*% beta_ave /sieve_w
  X_reg <- cbind(X_trans^0, X_trans^1)
  for (ii in 2:q){
    X_reg <- cbind(X_reg, X_trans^ii)
  }
  # Iterate_fun <- function(q_search) 
    # sum(log(1 + exp( pmin(t(X_reg) * q_search, 700)  ))) - sum(t(Y) %*% X_reg * q_search)
  # Iterate_result <- fminsearch(Iterate_fun, rep(1, q+1) ) 
  # coeff_q <- Iterate_result$xmin
  # dim(coeff_q) <- c(q+1,1)
  coeff_q <- Gradient_logistic_estimation(X_reg, Y, scaling = iter_times2)
  g_trans <- X_reg %*% coeff_q
  g_trans <- pmin(g_trans, 700) # limit the number inserted in exp()
  g_hat <- exp(g_trans)/(1+exp(g_trans)) # get g_hat
  
  coeff_prime_q = coeff_q[2]
  for (j in 3:(q+1)){
    coeff_prime_q <- cbind(coeff_prime_q, (j-1) * coeff_q[j])
  }
  dim(coeff_prime_q) <- c(q,1)
  g_prime_trans <- X_reg[,1:q] %*% coeff_prime_q / sieve_w
  g_prime_hat <- exp(g_trans)/(1+exp(g_trans))^2 * g_prime_trans # get g_prime_hat
  
  sigma_11 <- matrix(0, nrow = p , ncol = p )
  sigma_22 <- matrix(0, nrow = p , ncol = p )
  sigma_33 <- matrix(0, nrow = p , ncol = p )
  sigma_trans <- matrix(0 , nrow = q+1 ,ncol =p )
  
  for (j in 1:N){
    sigma_11 <- sigma_11 + g_hat[j] * (1 - g_hat[j]) * 
      matrix(X[j,], nrow = p) %*% matrix(X[j,], ncol = p)
  }
  sigma_11 <- sigma_11 / N 
  
  for (j in 1:N){
    sigma_trans <- sigma_trans + g_prime_hat[j] * 
      matrix(X_reg[j,], nrow = q+1) %*% matrix(X[j,], ncol = p)
  }
  simga_trans <- sigma_trans / N
  
  for (j in 1:N){
    sigma_22 <- sigma_22 + g_prime_hat[j] *
      matrix(X[j,], nrow = p) %*% matrix(X[j,], ncol = p) - 
      N * matrix(X[j,], nrow = p) %*% matrix(X_reg[j,], ncol = q+1) %*% 
      solve(t(X_reg) %*% X_reg) %*% sigma_trans
  }
  sigma_22 <- sigma_22 / N   
  
  for (j in 1:N){                  # this is sigma for test
    sigma_33 <- sigma_33 + g_prime_hat[j]  * 
      matrix(X[j,], nrow = p) %*% matrix(X[j,], ncol = p)
  }
  sigma_33 <- sigma_33 / N 
  
  
  sigma_sba <- t(solve(sigma_22[2:p,2:p])) %*% sigma_11[2:p,2:p] %*% solve(sigma_22[2:p,2:p])
  sigma_sba_test <- t(solve(sigma_33[2:p,2:p])) %*% sigma_11[2:p,2:p] %*% solve(sigma_33[2:p,2:p])
  st_sba <- (diag(sigma_sba))^0.5 / N^0.5
  st_sba_test <- (diag(sigma_sba_test))^0.5 / N^0.5
  
  # dim(st_sba) <- c(p-1, 1)
  beta_ave_oneless <- matrix(beta_ave[2:p], nrow = p-1) /beta_ave[1]
  pvalue_sba <- beta_ave_oneless / st_sba *beta_ave[1]
  report_list1 <- cbind(beta_ave_oneless, st_sba , st_sba_test )
  report_list2 <- coeff_q
  colnames(report_list1) <- c('coeff','st_error','test_st_error')
  colnames(report_list2) <- c('coeff_polynomial')
  report_list <- list(report_list1, report_list2)
# plot
  plot(X_trans/ beta_ave[1] * sieve_w, g_hat, col="red", xlab="index", ylab="probability_predict")
  points(X_trans/ beta_ave[1] * sieve_w, pnorm(X_trans/ beta_ave[1] * sieve_w), col = "blue")
  legend("topleft", legend=c("SBA_estimator", "Normal_distribution"),
         col=c("red", "blue"), lty=1:2, cex=0.8)
  
  return(report_list) #standardize the first coefficient into 1
  # return(coeff_q)
}

# simulation
# generate binary choice model (normal distribution)
R <- diag(8)
mu <- rep(0, dim(R)[1])
n <- 5000
tt <- rnorm(n)
dim(tt) <- c(n,1)
xx <- mvtnorm::rmvnorm(n, mean = mu, sigma = R)
yy <- (xx %*% matrix(c(1,2,4,5,-1,-2,-4,-5),nrow=8) + tt) >0
# get parameter value, standard error 
binary_sba <- SBA_estimator(xx,yy,iter_times = 1000, iter_times2 = 100 )
SBA_predict(binary_sba, xx[4980:5000,])
