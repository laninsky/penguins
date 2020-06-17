
#' The Mean Square Error (MSE) of all training sets .
#' 
#' The MSE is evaluated as
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
#' This is performed using \code{\link{estimateAlpha}}.
#' In order to parallel the validation procedure, this function can be redefined to call \code{\link{helpCValpha}} in parallel on groups.
#'    
#' @param SFS.groups The \code{nb.groups} nummber of training sets produced from the observed SFS.
#' @param nb.groups The number of training sets used for estimating \code{alpha}.
#' @param knots The time points used for the estimation of the spline using \code{\link{estimateCubSFS}}
#' @param theta.prm The initialization of the pararmeter to estimate. This is the differences 
#'   \eqn{a_i-a_{i-1}} where \eqn{a_i} is the value of the spline at time point \eqn{t_i}.
#' @param alpha The smoothing parameter to minimize over.
#' @param n.samples The number of sequences.
#' @param Wmat,Vmat The Polanski and Kimmel matrices 
#' @param is.folded Is the SFS folded(TRUE) or unfolded (FALSE).
#' @param loc.solv The local solver. May be either \code{COBYLA} or \code{LBFGS}.
#' @inheritParams nloptr::auglag
#'  
#' @return An integer. The MSE over all validation sets given in \code{SFS.groups} given \code{alpha}.
#' 
#' @seealso \code{\link{estimateAlpha}} \code{\link{helpCValpha}} 
#' 
#' @references 
#' Polanski, A. and Kimmel, M. (2003), "New Explicit Expressions for Relative Frequencies of Single-Nucleotide Polymorphisms With Application 
#' to Statistical Inference on Population Growth." Genetics, 165, 427-436.
#' 
#' @export

CValpha <-
function(SFS.groups,nb.groups,knots,theta.prm,alpha,n.samples,Wmat,Vmat,is.folded,
         nb.time=500,conv.lim=10^{-3},loc.solv="COBYLA",
                    control=list(maxeval=20000, ftol_abs=10^(-5),stopval=0)){
  print("CValpha reporting for duty")

  fit <- mclapply(split(1:nb.groups,as.factor(1:nb.groups)),helpCValpha,
                  SFS.groups=SFS.groups,n.samples=n.samples,knots=knots,theta.prm=theta.prm,alpha=alpha,
                  Wmat=Wmat,Vmat=Vmat,is.folded=is.folded,nb.time=nb.time,conv.lim=conv.lim,loc.solv=loc.solv,control=control, mc.cores=numcores)

  mean(unlist(fit))
}
