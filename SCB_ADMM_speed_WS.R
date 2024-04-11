#' bi-ADMM: a Biclustering Algorithm for the General Model (faster version with one step iteration)
#'
#' Same algorithm as \code{biADMM}. Call python code to speed up the running time.
#' @param X The data matrix to be clustered. The rows are the samples, and the columns are the features.
#' @param nu1 A regularization parameter for row shrinkage
#' @param nu2 A regularization parameter for column shrinkage
#' @param gamma_1 A regularization parameter for row shrinkage
#' @param gamma_2 A regularization parameter for column shrinkage
#' @param m m-nearest-neighbors in the weight function
#' @param phi The parameter phi in the weight function
#' @param niters Iteraion times
#' @param tol Stopping criterion
#' @param output When output = 1, print the results at each iteration. No print when output equals other value.
#'
#'
#' @return A list of results, containing matrix of A, v, z, lambda1, and lambda2
#' @export
#'
#' @examples
#' # generate dataset
#' set.seed(123)
#' X = data_gen(n = 100, p = 80 , true_p = 40)
#' # set parameters
#' nu1 = nu2 = nu3 = gamma_1 = gamma_2 = gamma_3 = 0.1
#' m = 5
#' phi = 0.5
#' # biADMM algorithm
#' res2 = sparse.biADMM.speed.WS(X, nu1, nu2, gamma_1, gamma_2,
#'  m, phi, niter = 10, tol = 0.0001, output = 0)
#' res1 = sparse.biADMM.speed.WS(X, A= NULL, nu1, nu2, nu3,v=NULL ,z=NULL ,g=NULL ,gamma_1, gamma_2, gamma_3,
#'                              feature_weight = feature_weight, m , phi,niter = 2000,tol = tol,output = 0)
#' dim(res2$A)

sparse.biADMM.speed.WS = function(X, A= NULL, nu1, nu2, nu3,v=NULL ,z=NULL ,g=NULL ,
                               gamma_1, gamma_2, gamma_3,
                               feature_weight = rep(1,ncol(X)),
                               m = 5, phi=0.5,niter = 1000,tol = 0.1,output = 1){

  path <- paste("./R", "SBC_ADMM.python.WS.py", sep="/") #niters=1
  source_python(path)
  # head(A)
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  k_row <- m; k_col <- m
  w_row <- kernel_weights(t(X), phi)
  w_col <- kernel_weights(X, phi)
  w_row <- knn_weights(w_row, k_row, n)
  w_col <- knn_weights(w_col, k_col, p)
  w_row <- w_row/sum(w_row)
  w_col <- w_col/sum(w_col)
  w_row <- w_row/sqrt(p)
  w_col <- w_col/sqrt(n)
  
  w_l <- w_row
  u_k <- w_col
  w_l <- matrix(w_l, length(w_l),1)
  u_k <- matrix(u_k, length(u_k),1)
  
  feature_weight <- feature_weight/ sum(feature_weight)/sqrt(p)
  feature_weight <- matrix(feature_weight,length(feature_weight),1)
  #res <- SCB_ADMM_python_AR(X,A, nu1, nu2, nu3, gamma_1, gamma_2,gamma_3, 
  #                       w_l, u_k, feature_weight, niter, tol, output = output)
  
  res <- SCB_ADMM_python_WS(X,A,v, z, g, nu1, nu2, nu3, gamma_1, gamma_2, gamma_3,
                            w_l, u_k, feature_weight,  tol, niter,output = output)
  
  result <- list(A = res[[1]], v = res[[2]], z = res[[3]], g = res[[4]],
                 lambda_1 = res[[5]], lambda_2 = res[[6]], lambda_3 = res[[7]],
                 gamma_1 =res[[8]],gamma_2 =res[[9]],gamma_3=res[[10]],iters = res[[11]])
  return(result)
}