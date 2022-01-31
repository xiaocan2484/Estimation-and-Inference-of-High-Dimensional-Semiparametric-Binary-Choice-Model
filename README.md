# Estimation-and-Inference-of-High-Dimensional-Semiparametric-Binary-Choice-Model
This is code for my job market paper: Estimation and Inference of High Dimensional Semiparametric Binary Choice Model

https://drive.google.com/file/d/19ru2U36whsvIEbjGpWwG0CLtX5eS8Uyg/view

Abstract:
Binary choice models can be easily estimated (using, e.g. maximum likelihood estimation) when the distribution of the latent errors is known, as in Logit or Probit. In contrast, most estimators with unknown error distribution (e.g., maximum score, maximum rank correlation, or Klein-Spady) are computationally difficult or numerically unstable, making estimation impractical with more than a few regressors. This paper proposes an estimator that is convex at each iteration, and so is numerically well behaved even with many regressors and large sample sizes. The proposed estimator, which is root-n consistent and asymptotically normal, is based on batch gradient descent, while using a sieve to estimate the unknown error distribution function. In high dimensional setting, the estimator is square root of p/n consistent when p/n goes to 0  and asymptotic normal when p^2/n goes to 0 , where p is the number of regressors and n is the number of observations. An application of proposed estimator to predicting bankruptcy is provided.

I offer code of Matlab version and R version

We can use machine learning to estimate error distribution such as polynomials, random forrest, neural network.
