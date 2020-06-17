#' Estimating the smoothing parameter.
#' 
#' First generate the \code{nb.groups} groups of sites, each of which are to function as a validation set.
#' Estimate the smoothing parameter \code{alpha} by minimizing the Mean Square Error calculated by \code{\link{CValpha}}.
#' \deqn{\frac{1}{K}\sum\limits_{k=1}^{K}\sum\limits_{i=1}^{n-1} \frac{\left(\xi_i^{(k)}-\E[\xi_i^{(k)}]\right)^2}{\E[\xi_i^{(k)}]},}{1/K \sum_{k=1}^{K} \sum_{i=1}^{n-1} ( \xi_i^(k) - E[\xi_i^(k)] )^2 / E[\xi_i^(k)],}
#' where \eqn{K} is the number of groups of sites used in the cross validation, 
#' \eqn{\xi^{(k)}}{\xi^(k)} is the \eqn{k}'th validation SFS, and 
#' \eqn{E[\xi^{(k)}]}{E[\xi^(k)]} is the expected SFS based on the spline \eqn{\Lambda(t)} 
#' inferred from the \eqn{K-1} training sets, not including set \eqn{k}.
#' 
#' The spline is given as
#' \deqn{\Lambda = a_i + b_i (t-t_i) + c_i (t-t_i)^2 + d_i (t-t_i)^3 for t_i\leq t < t_{i+1},}{\Lambda = 
#' a_i + b_i (t-t_i) + c_i (t-t_i)^2 + d_i (t-t_i)^3 for t_i\le t < t_{i+1},}
#' where \eqn{\Lambda(t) = a_m + b_m (t-t_m)} for \eqn{t_m \leq t}{t_m \le t}, \eqn{\Lambda'(0)=b_0=1} and \eqn{\Lambda(0)=0}.
#' 
#' For given time points \eqn{t_i} for \eqn{i=0,1,\ldots,m} and  values \eqn{a_0=\Lambda(0)=0, \Lambda(t_1) = a_1, \ldots, \Lambda(t_m)=a_m},
#' the spline is completly specified. Hence, the score function is minimized wrt. \eqn{a_0, a_1, \ldots, a_m} given some time points.
#' A reparameterisation is used \eqn{\psi_i=a_{i}-a_{i-1}} for \eqn{i=1,2,\ldots,m}.
#'
#' Given a spline \eqn{\Lambda(t)} and the smoothing parameter \eqn{\alpha}, 
#' the score function is
#' \deqn{S(\Lambda)= (1-\alpha)\sum\limits_{i=1}^{n-1}\frac{\left(\E[\xi_i]-\xi_i\right)^2}{\E[\xi_i]}+
#' \alpha\int\limits_{0}^{\infty}\left(\Lambda''(t)\right)^2dt.}{S(\Lambda) = (1-\alpha) * 
#' \sum_{i=1}^{n-1} ( E[\xi_i] - \xi_i )^2 / E[\xi_i] +
#' \alpha * intergral_{0}^{\infty} ( \Lambda''(t) )^2 dt.}
#'
#' The probability of observing a site with \eqn{i} derived alleles is evaluated using the formula by Polanski and Kimmel (2003)
#' \deqn{p_i(\Lambda)=\frac{\sum\limits_{j=2}^ne_j (\Lambda)W_{ij}}{\sum\limits_{j=2}^n e_j(\Lambda)V_j}}{p_i(\Lambda) = ( \sum_{j=2}^n e_j(\Lambda) W_{ij} ) / ( \sum_{j=2}^n e_j(\Lambda) V_j )}
#' for \eqn{i=1,2,\ldots,n-1} for unfolded SFS (\code{is.folded=F}). For the folded SFS the probability of observing a site with minor 
#' allele count \eqn{i} among \eqn{n} sequences is 
#' \deqn{p_i^F(\Lambda) = p_i(\Lambda) + p_{n-i}(\Lambda)}
#' for \eqn{i=1,\ldots,n/2} if \eqn{n} is even and for \eqn{n} odd
#' \deqn{p_i^F(\Lambda) = p_i(\Lambda) + p_{n-i}(\Lambda)}
#' for \eqn{i=1,\ldots,(n-1)/2} and \eqn{p_{(n+1)/2}^F(\Lambda) = p_{(n+1)/2}(\Lambda)}.
#' 
#' Here \eqn{e_j(\Lambda)} is evaluated using \code{\link{Ejfct}}.
#' 
#' The expected SFS \eqn{E[\xi^{(k)}]}{E[\xi^(k)]} is estimated by the total number of segregating sites in group \eqn{k} 
#' multiplied by the above probabilities.
#' 
#' The roughness penalty is evaluated by a linear equation system, implemented in \code{\link{SplineFct}}.
#' 
#' Minimizing this function w.r.t. \code{alpha} will provide an estimate of the smooting parameter.
#' 
#' In order to parallel the cross validation procedure, the \code{\link{CValpha}} function can be redefined to call \code{\link{helpCValpha}} in parallel on groups.
#' 
#' The minimum of the MSE is found by using \code{\link[stats]{optimize}} with tolerance \code{tol}. 
#' 
#' @param rawSFS the observed SFS. 
#' @param n.samples the number of sequences used to generate the observed SFS.
#' @param knots the time points \eqn{t_i} for \eqn{i=0,\ldots,m} where \eqn{t_0=0}.
#' @param theta.prm the parameters to minize over where \code{theta.prm[i]}\eqn{=psi[i]=a_i-a_{i-1}} for \eqn{i=1,2,\ldots,m}, where \eqn{\Lambda(t_i)=a_i}.
#' @param nb.groups the number of training sets to devid the observed SFS into. must be larger than 2.
#' @param Wmat,Vmat the Polanski and Kimmel matrices.
#' @param is.folded is the SFS folded(True) or unfolded (False).
#' @param nb.time the number of time points used to estimate \eqn{e_j(\Lambda)}
#' @param tol: the desired accurarcy in \code{\link{optimize}}.
#' @param loc.solv the local solver. May be either \code{COBYLA} or \code{LBFGS}.
#' @param conv.lim the desired accuracy of the paraemters.
#' @inheritParams nloptr::auglag
#' 
#' @return a list with two elements. The result from calling optimize.
#' \describe{
#'   \item{minimum}{The estimated smoothing parameter}
#'   \item{value}{The value of the MSE given the estimated smoothing parameter}
#' }
#'
#' @seealso \code{\link[stats]{optimize}} \code{\link{CValpha}}
#' 
#' @references 
#' Polanski, A. and Kimmel, M. (2003), "New Explicit Expressions for Relative Frequencies of Single-Nucleotide Polymorphisms With Application 
#' to Statistical Inference on Population Growth." Genetics, 165, 427-436.
#' 
#' @examples 
#' 
#' rawSFS <- c(1925,867,532,372,280,221,181,151,129,111,98,87,77,70,63,58,53,49,45,42,39,36,34,32,30,28,27,25,24,23,22,21,20,19,18,17,17,16,15,15,14,14,13,13,12,12,12,11,11)
#' n.samples <- 50
#' Wmat <- WFunction(n.samples)
#' Vmat <- VFunction(n.samples)
#' 
#' ## depending on your system this may take a minut.
#' estimateAlpha(rawSFS,n.samples,seq(0,1,len=5),seq(0,1,len=5)[-1],5,Wmat,Vmat,F,tol=0.1,conv.lim=0.1)
#' 
#' @export

estimateAlpha <-
function(rawSFS,n.samples,knots,theta.prm,nb.groups,Wmat,Vmat,is.folded,nb.time=500,
                           tol=10^(-2),conv.lim=10^{-3},
                           loc.solv="COBYLA",
                           control=list(maxeval=20000, ftol_abs=10^(-5),stopval=0)){

  ## generate the training sets.
  SFS.groups <- matrix(0,ncol=nb.groups,nrow=length(rawSFS))

  rand.sample <- sample(rep(1:length(rawSFS),rawSFS))

  group.size <- floor(sum(rawSFS)/nb.groups)
  rest.size <- sum(rawSFS) %% nb.groups
  group.size <- c(rep(group.size+1,rest.size),rep(group.size,nb.groups-rest.size))
  group.index <- c(0,cumsum(group.size))

  print("up to the group step")
  for(i in 1:nb.groups){

    sub.SFS <- rand.sample[(group.index[i]+1):group.index[i+1]]
    SFS.groups[,i] <- as.vector(tabulate(sub.SFS,nbins=length(rawSFS)))
  }

  print("about to call CValpha in optimise")
  ## estimate alpha given number of knots.
  min.alpha <- optimise(CValpha,lower=0,upper=1,tol=tol,
                        SFS.groups=SFS.groups,nb.groups=nb.groups,
                        knots=knots,theta.prm=theta.prm,n.samples=n.samples,
                        Wmat=Wmat,Vmat=Vmat,is.folded=is.folded,nb.time=nb.time,conv.lim=conv.lim,
                        loc.solv=loc.solv,control=control)
  return(min.alpha)
}
