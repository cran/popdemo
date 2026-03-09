################################################################################
#' Project population dynamics 
#'
#' @description
#' Analyse long-term dynamics of a stochastic population matrix projection model.
#'
#' @param A a list of matrices. \code{stoch} uses \code{\link{project}} to 
#' perform a stochastic' projection where the matrix varies with each timestep. 
#' The sequence of matrices is determined using \code{Aseq}. Matrices must be 
#' square, non-negative and numeric, and all matrices must have the same dimension. 
#' @param what Deprecated: later versions will always return all
#' values of growth, log growth and variance in growth. This argument determines 
#' what should be returned. A character vector with possible entries
#' "lambda" (to calcualate stochastic growth), "log_lambda" (to calcualate logarithm 
#' of stochastic growth), "var" (to calculate variance in log stochastic growth) 
#' and/or "all" (to calculate all values).
#' @param vector (optional) a numeric vector describing the age/stage 
#' distribution used to calculate the projection. If \code{vector} is not 
#' specified, a random vector is generated. Long-term stochastic dynamics should 
#' usually be the same for any vector, although if all the matrices in A are 
#' reducible (see \code{\link{isIrreducible}}), that may not be the case.
#' @param Aseq the sequence of matrices in a stochastic projection. 
#' \code{Aseq} may be either:
#' \itemize{
#'  \item "unif" (default), which results in every matrix in \code{A} having an 
#'  equal, random chance of being chosen at each timestep.
#'  \item a square, nonnegative left-stochastic matrix describing a 
#'  first-order markov chain used to choose the matrices. This should have the 
#'  same dimension as the number of matrices in \code{A}. 
#'  \item a numeric vector giving a specific sequence which corresponds to the
#'  matrices in \code{A}.
#'  \item a character vector giving a specific sequence which corresponds to the
#'  names of the matrices in \code{A}.
#' }
#' @param Astart (optional) in a stochastic projection, the matrix with which to
#' initialise the projection (either numeric, corresponding to the matrices in 
#' \code{A}, or character, corresponding to the names of matrices in \code{A}). 
#' When \code{Astart = NULL}, a random initial matrix is chosen.
#' @param iterations the number of projection intervals. The default is 1e+5.
#' @param discard the number of initial projection intervals to discard, to 
#' discount near-term effects arising from the choice of vector. The default is
#' 1e+3
#' @param PREcheck many functions in \code{popdemo} first check Primitivity, 
#' Reducibility and/or Ergodicity of matrices, with associated warnings and/or 
#' errors if a matrix breaks any assumptions. Set \code{PREcheck=FALSE} if you
#' want to bypass these checks.
#'
#' @details 
#' Calculates stochastic growth and its variance for a given stochastic population 
#' matrix projection model.
#'
#' @return 
#' A numeric vector with three possible elements: "lambda" (the geometric mean stochastic population
#' growth rate) "log_lambda" (the arithmetic mean of the logarithm of the stochastic population growth rate) 
#' and "var" (the variance of the logarithm of the stochastic population growth rate). Values
#' returned depend on what's passed to \code{what}.
#'
#' @examples
#'   # load the Polar bear data
#'   ( data(Pbear) )
#'
#'   # Find the stochastic info for a time series with uniform probability of each
#'   # matrix
#'   ( all_unif <- stoch(Pbear, what = "all", Aseq = "unif") )
#'
#'   # Find the stochastic growth for a time series with uniform probability of each
#'   # matrix
#'   ( lambda_unif <- stoch(Pbear, what = "lambda", Aseq = "unif") )
#'
#'   # Find the variance in stochastic growth for a time series with uniform 
#'   # probability of each matrix
#'   ( var_unif <- stoch(Pbear, what = "var", Aseq = "unif") )
#'                  
#'   # Find stochastic growth and its variance for a time series with a sequence of
#'   # matrices where "bad" years happen with probability q
#'   q <- 0.5
#'   prob_seq <- c(rep(1-q,3)/3, rep(q,2)/2)
#'   Pbear_seq <- matrix(rep(prob_seq,5), 5, 5)
#'   ( all_q <- stoch(Pbear, what = "all", Aseq = Pbear_seq) )
#'   
#' @concept stochastic growth
#' @concept variance
#' @concept projection
#' @concept project
#' @concept population
#'
#' @export stoch
#' @import methods
#'
stoch<-
function (A, what = "all", Aseq = "unif", vector = NULL, 
          Astart = NULL, iterations = 1e+4, discard = 1e+3, 
          PREcheck = FALSE){
## check structure of A and rearrange into an array
if(!any((is.list(A) & all(sapply(A, is.matrix))), 
        (is.array(A) & length(dim(A)) == 3)) ){
    stop("A must be a list of matrices")
}
if(is.list(A) & length(A) == 1) stop("A must contain more than one matrix")
if (is.list(A) & length(A) > 1) {
    numA <- length(A)
    alldim <- sapply(A, dim)
    if (!diff(range(alldim))==0) {
        stop("all matrices in A must be square and have the same dimension as each other")
    }
    dimA <- mean(alldim)
    L <- A
    A <- numeric(dimA * dimA * numA); dim(A) <- c(dimA, dimA, numA)
    for(i in 1:numA){
        A[,,i] <- L[[i]]
    }
    dimnames(A)[[1]] <- dimnames(L[[1]])[[1]]
    dimnames(A)[[2]] <- dimnames(L[[1]])[[2]]
    dimnames(A)[[3]] <- names(L)
}
## extract some info about matrices
order <- dim(A)[1]
##what should be calculated?
ifelse("lambda" %in% what, growth <- TRUE, growth <- FALSE)
ifelse("log_lambda" %in% what, log_growth <- TRUE, log_growth <- FALSE)
ifelse("var" %in% what, variance <- TRUE, variance <- FALSE)
if("all" %in% what){
    growth <- TRUE
    log_growth <- TRUE
    variance <- TRUE
}
if(!growth & !log_growth & !variance) stop('"what" does not contain the right information')
if (!missing(what)) {
    warning(
        "Argument 'what' is deprecated and later versions will always return all\n
        values of growth, log growth and variance in growth."
    )
}## check vector
if(!is.null(vector) & length(vector) != order){
    stop("vector must be equal to dimension of A")
}
if(is.null(vector)) vector <- stats::runif(order)
vector <- vector / sum(vector)
##project the model
pr <- project(A = A, vector = vector, time = iterations, standard.A = FALSE,
              standard.vec = FALSE, return.vec = FALSE, 
              Aseq = Aseq, PREcheck = PREcheck)
##calculate the stuff
#work out per-timestep growth
gr <- pr[(discard:iterations) + 1] / pr[discard:iterations]
# find the per-timestep mean growth (stochastic growth)
if (growth) gr_mean <- exp(mean(log(gr)))
# fin the per-timestep mean log growth
if (log_growth) log_gr_mean <- mean(log(gr))
# find the per-timestep variance in growth
if (variance) gr_var <- stats::var(log(gr))
if (growth & log_growth & variance) {
    final <- data.frame(lambda = gr_mean, log_lambda = log_gr_mean, var = gr_var)
}
if (growth & log_growth & !variance) {
    final <- data.frame(lambda = gr_mean, log_lambda = log_gr_mean)
}
if (growth & !log_growth & variance) {
    final <- data.frame(lambda = gr_mean, var = gr_var)
}
if (!growth & log_growth & variance) {
    final <- data.frame(log_lambda = log_gr_mean, var = gr_var)
}
if (growth & !log_growth & !variance) final <- data.frame(lambda = gr_mean)
if (!growth & log_growth & !variance) final <- data.frame(log_lambda = log_gr_mean)
if (!growth & !log_growth & variance) final <- data.frame(var = gr_var)
return(final)
}

# decompose growth into lambda vs transient
# if(decomp){
#     #find the per-timestep growth decomposition
#     gr_a <- apply(A[, , MC], 3, eigs, what = "lambda")
#     gr_t <- gr / gr_a
#     gr_decomp <- cbind(asymptotic = gr_a, transient = gr_t)
#     #find the mean growth components
#     if(growth) gr_decomp_mean <- colMeans(gr_decomp)
#     #find the vcv of the decomposed components
#     if(variance) gr_decomp_var <- var(gr_decomp)
#     if(growth & variance) final <- list(lambda = gr_mean, 
#                                         var = gr_var, 
#                                         lambda_decomp = gr_decomp_mean,
#                                         var_decomp = vcv)
#     if(growth & !variance) final <- gr_mean
#     if(!growth & variance) final <- gr_var
# }
