##########################################
## Chetverikov, Kim, and Wilhelm (2018) ##
##             functions                ##
##########################################

# kernel function (Epanechnikov)
K <- function(u, h) {0.75*(1-(u/h)^2)/h  * ( abs(u/h) <= 1 )}

# functions 
K_ixh <- function(xhgrid, X) {
  xxgrid <- expand.grid(X, X)
  res    <- colSums(matrix( sign( xxgrid[,2] - xxgrid[, 1] )*K( (xxgrid[,2] - xhgrid[1]) , xhgrid[2] )*K( (xxgrid[,1] - xhgrid[1]), xhgrid[2] ), nrow = length(X)))
  return(res)
}

Y_ind <- function(ygrid, y) {(y <= ygrid)}

T_fn <- function(k, y, Yindex) colSums( k*y[, Yindex] )/sqrt( colSums( k^2 ) ) 

b_fn <- function(e, k, y, f, Yindex) max( colSums( e*k*(y[, Yindex] - f[, Yindex]) )/sqrt( colSums( k^2 ) ) )

Tb_fn <- function(Yindex, k, f, y, e) {
  val   <- max(sapply( Yindex, b_fn, simplify = TRUE, e = e, k = k, f = f, y=y ) )
  return(val)
}

CKW.test <- function(Y, X, brep = 200, num.grid = 100, u = 2/3, delta = 1/2, alpha = 0.95) {
  
  # grid setup
  n     <- length(Y) # sample size
  hmin  <- 1/(n)^(delta) # mininum of bandwidth values
  l     <- seq(0, round(log(1/hmin)/log(1/u)), 1) # number of bandwidth values
  hgrid <- (u)^l #bandwidth - geometric progress
  xgrid <- X # X grid
  ygrid <- seq(min(Y), max(Y), length.out = num.grid) # Y grid
  
  e  <- matrix(rnorm(n*brep), ncol = brep) # multipliers for bootstrap from standard normal distribution
  
  xhgrid <- expand.grid(xgrid, hgrid)
  K_i    <- 2*matrix( apply(xhgrid, 1, K_ixh, X), nrow = n )
  yy     <- sapply( ygrid, Y_ind, simplify = TRUE, Y)
  Yindex <- seq(1, num.grid, 1)
  
  T_st   <- (sapply( Yindex, T_fn, simplify = TRUE, k = K_i, y = yy ))
  
  T_stat <- max(T_st)
  
  T_st_ind <- rowSums(T_st == T_stat) > 0
  
  # GM bootstrap critical value
  xxgrid    <- expand.grid(X, X)
  mat_denom <- matrix(K( (xxgrid[,2] - xxgrid[,1]), 1/sqrt(n) ), nrow = n)
  np_denom  <- colSums(mat_denom)
  fits      <- t(mat_denom)%*%yy/np_denom
  
  boot_dist <- foreach(j = 1:brep, .export = "fits") %do% Tb_fn(Yindex, K_i, fits, yy, e[,j])
  boot_dist <- unlist(boot_dist)
  
  bw   <- sum(xhgrid[,2]*T_st_ind)
  cval <- quantile(boot_dist, alpha)
  rej  <- (T_stat > cval)
  result <- c(rej, T_stat, cval, bw)
  names(result) <- c("rej", "T_stat", "cval", "bw")
  
  return(result)
}