\name{RcppArmadillo-STPPMx-functions}
\alias{replace_cpp}
\alias{dmvt_cpp}
\alias{centering_cpp}
\alias{partitions_draw_cpp}
\alias{partitions_draw_cpp_fast}
\alias{predict_cpp_new}

\title{C+ STPPMX functions}
\description{
  Main MCMC Neal's algorithm, predictive function and other useful function implemented in C++ using Armadillo.
}
\usage{
replace_cpp(x)
dmvt_cpp(x, v, mu, Sig)
partitions_draw_cpp(y, X_sim, X_reg, M, alpha, c, v, s, sigma, beta, m_0, Phi, w, time_effects, space_effects, mixed_effects)
partitions_draw_cpp_fast(y, X_sim, X_reg, M, alpha, c, v, s, sigma, beta, m_0, Phi, w, time_effects, space_effects, mixed_effects, sizes_mem, sums_mem, products_mem)
predict_cpp_new(obs, station, time, X_train_sim, x_futur_sim, x_futur_reg, y_test, beta_list , sigma_list, space_effects_list, time_effects_list, mixed_effects_list, clusters_sizes_list, partition_list, alpha_vec, nb_density, nb_draw, w, Phi, v, s, m_0)
}
\arguments{
  \item{x}{a numeric vector}
  \item{v}{Student's t degrees of freedom}
  \item{mu}{Student's t mean vector}
  \item{Sig}{Student's t covariance matrix}

  \item{y}{vector of observations}
  \item{X_sim}{matrix with rows being observations and columns being covariates used in similarity}
  \item{X_reg}{matrix with rows being observations and columns being covariates used in regression (without intercept columns of 1's)}
  \item{M}{Number of Neal's algorithm auxiliary parameters}
  \item{alpha}{Dirichlet process concentration parameter}
  \item{c}{vector of current partition with cluster indexes}
  \item{v}{gamma prior shape parameter}
  \item{s}{gamma prior scale parameter}
  \item{sigma}{vector of current clusters sigma variance parameters}
  \item{beta}{matrix of current clusters beta parameters as columns}
  \item{m_0}{multivariate normal prior mean parameter}
  \item{Phi}{Inverse Wishart prior scale matrix parameter}
  \item{w}{Inverse Wishart prior degree of freedom}
  \item{time_effects}{vector of random time effects (vector of 0's if model doesn't use)}
  \item{space_effects}{vector of random space effects (vector of 0's if model doesn't use)}
  \item{mixed_effects}{vector of random spatio-temporal effects (vector of 0's if model doesn't use)}
  \item{sizes_mem}{vectors of clusters sizes kept in memory at each MCMC iteration if fast = TRUE}
  \item{sums_mem}{matrix of clusters related sums kept in memory at each MCMC iteration if fast = TRUE}
  \item{products_mem}{array of clusters related matrix products kept in memory at each MCMC iteration if fast = TRUE}

  \item{obs}{number indexing which test observation to predict}
  \item{station}{which station is observation from}
  \item{time}{which time point is the observation at}
  \item{X_train_sim}{training similarity matrix}
  \item{x_futur_sim}{similarity covariates vector of the test observation to predict}
  \item{x_futur_reg}{regression covariates vector of the test observation to predict}
  \item{y_test}{test observed value that we want to predict}
  \item{beta_list}{List with all the MCMC iterations beta matrices}
  \item{sigma_list}{List with all the MCMC iterations sigma vectors}
  \item{time_effects_list}{List with all the MCMC iterations time effects vectors}
  \item{space_effects_list}{List with all the MCMC iterations space effects vectors}
  \item{mixed_effects_list}{List with all the MCMC iterations mixed effects vectors}
  \item{clusters_sizes_list}{List with all the MCMC iterations clusters sizes vectors}
  \item{partition_list}{List with all the MCMC iterations cluster indexes vectors}
  \item{alpha_vec}{vector with all the MCMC iterations alpha parameters}
  \item{nb_density}{number of saved MCMC iterations density parameters (equivalent to number of iterations after burn_in)}
  \item{nb_draw}{number of draws to approximate the CRPS}
}
\value{
  \code{replace_cpp(x)} replace the number vector by its ordered index from 1:N, N being that x vector length.

  \code{dmvt_cpp(x, v, mu, Sig)}  Evaluate the multivariate student t density with mean \code{mu}, covariance matrix \code{Sig} and degree of freedom \code{v}
  at \code{x}.

  \code{partitions_draw_cpp(y, X_sim, X_reg, M, alpha, c, v, s, sigma, beta, m_0, Phi, w, time_effects, space_effects, mixed_effects)} Code used to draw random partition according to the modified Neal algorithm, used in the STPPMx model, then updates the associated parameters (citer memoire). Returns a list of the updated partition and parameters.

  \code{partitions_draw_cpp_fast(y, X_sim, X_reg, M, alpha, c, v, s, sigma, beta, m_0, Phi, w, time_effects, space_effects, mixed_effects, sizes_mem, sums_mem, products_mem)} Fast version of \code{partitions_draw_cpp}, based on using more memory to skip some heavy marix computation updates that are redundant when iterating over the partition.

  \code{predict_cpp_new(obs, station, time, X_train_sim, x_futur_sim, x_futur_reg, y_test, beta_list , sigma_list, space_effects_list, time_effects_list, mixed_effects_list, clusters_sizes_list, partition_list, alpha_vec, nb_density, nb_draw, w, Phi, v, s, m_0)}
} Prediction function used to obtain approximate CRPS and AE from an STPPMX model return object and future dataset that we want to predict.
\details{
Main MCMC Neal's algorithm and predictive function implemented in C++ using Armadillo with some other useful C++ coded function, brought to R with
  RcppArmadillo for faster implementation. All of those functions were designed for my master's thesis, which I have referenced below. Therefore, all details concerning priors, appproximations, random effects models, etc. are in that paper.
}
\references{
Reference to my paper: Dounas, R. (2025). Prévisions météorologiques par régression non-paramétrique bayésienne basée sur le modèle PPMx spatio-temporel. [Master's thesis, Université de Montréal]. Papyrus. URL: https://hdl.handle.net/1866/41129.

Reference to Armadillo library: Conrad Sanderson and Ryan Curtin. Armadillo: An Efficient Framework for Numerical Linear Algebra. International Conference on Computer and Automation Engineering, 2025. https://doi.org/10.1109/ICCAE64891.2025.10980539.
}
\author{Rayane Dounas}
