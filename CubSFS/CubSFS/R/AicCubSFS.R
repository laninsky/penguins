#' Calculate the Akaike Information Criteria (AIC).
#' 
#' On the basis of an result from the \code{\link{estimate.CubSFS}} function, calculate the Akaike Information Criteria as
#' \deqn{2\times2m-2\log\hat{L},}{2*2m - 2*log L,}
#' where \eqn{m} is the number of time points used in the \code{\link{estimateCubSFS}} call, and \eqn{\hat{L}}{L} is the probability of 
#' \code{rawSFS} under the multinomial distribution with probability vector given by the estimated model. 
#' The probability vector is either based on the result of calling \code{\link{estimateCubSFS}} on \code{rawSFS}, 
#' or by the pointwise median of the result from the bootstraped samples.
#' 
#' @param rawSFS the observed SFS
#' @param n.samples the sample size, i.e. the number of sequences
#' @param SFS.res the result from \code{\link{estimateCubSFS}} based on the observed SFS
#' @param boot.res a list of the results from \code{\link{estimateCubSFS}} on each bootstrapped SFS
#' @param is.folded indicates if the SFS is folded (TRUE) or unfolded (FALSE)
#' @param nb.time number of time points for the fine partitioning of the time \code{knots} in \code{SFS.res}, Default=500.
#' 
#' @details 
#' 
#' The probability vector is evaluated by the formula provided by Polanski and Kimmel (2003) 
#' 
#' \deqn{p_i(\Lambda)=\frac{\sum\limits_{j=2}^ne_j (\Lambda)W_{ij}}{\sum\limits_{j=2}^n e_j(\Lambda)V_j}}{p_i(\Lambda)=(\sum_{j=2}^n e_j(\Lambda) W_{ij})/(\sum_{j=2}^n e_j(\Lambda) V_j)}
#' for \eqn{i=1,2,\ldots,n-1} for unfolded SFS (\code{is.folded=F}). For the folded SFS the probability of observing a site with minor 
#' allele count \eqn{i} among \eqn{n} sequences is 
#' \deqn{p_i^F(\Lambda)=p_i(\Lambda)+p_{n-i}(\Lambda)}
#' for \eqn{i=1,\ldots,n/2} if \eqn{n} is even and for \eqn{n} odd
#' \deqn{p_i^F(\Lambda)=p_i(\Lambda)+p_{n-i}(\Lambda)}
#' for \eqn{i=1,\ldots,(n-1)/2} and \eqn{p_{(n+1)/2}^{F}(\Lambda)=p_{(n+1)/2}(\Lambda)}
#' 
#' When the AIC is based on the call of \code{\link{estimateCubSFS}} on \code{rawSFS} 
#' \deqn{e_j(\Lambda)=\int_0^{\infty}exp(-{j\choose 2}\Lambda(t))dt}{e_j(\Lambda)=intergral_0^\infty exp( -choose(j,2) \Lambda(t) ) dt} 
#' is estimated by assuming piece wise straight lines between any two points. 
#' The \code{nb.time} points used to specify the straight lines are a fine partitioning of the time points used in the 
#' \code{estiamteCubSFS} call.
#' 
#' When the AIC is based on a list of bootstrap results, the point wise median of the resulting splines is generated. A total \code{nb.time} 
#' points are used between 0 and \code{max(knots)} from \code{SFS.res}. Then \eqn{p_i(\Lambda)} is evaluated by assuming straight lines betweeen
#' successive points of the piece wise median.
#'
#' @return a vector with two entries.
#' The AIC based on the Spline generated from the observed SFS, and 
#' the AIC based on the point wise median of the splines based on the bootstrap samples.
#' 
#' @references 
#' Polanski, A. and Kimmel, M. (2003), "New Explicit Expressions for Relative Frequencies of Single-Nucleotide Polymorphisms With Application 
#' to Statistical Inference on Population Growth." Genetics, 165, 427-436.
#' 
#' @seealso \code{\link{estimateCubSFS}}
#' 
#' @export

AicCubSFS <-
function(n.samples,SFS.res,boot.res=NULL,is.folded=F,nb.time=500){
  
  Wmat <- WFunction(n.samples)
  Vmat <- VFunction(n.samples)

  evec.spl <- EjFct(cumsum(SFS.res$ResSpline$par),SFS.res$knots,n.samples,nb.time=nb.time)

  probSFS <- (Wmat%*%evec.spl/sum(Vmat*evec.spl))
  if(is.folded){
    UFprobSFS <- probSFS
    if(length(UFprobSFS) %% 2 == 0){
      probSFS <- (UFprobSFS +rev(UFprobSFS))[1:(length(UFprobSFS)/2)]
    } else{
      probSFS <- (UFprobSFS +rev(UFprobSFS))[1:(floor(length(UFprobSFS)/2))]
      probSFS <- c(probSFS,UFprobSFS[(length(UFprobSFS)+1)/2])
    }
  }

  log.prob.data <- dmultinom(SFS.res$SFS,prob=probSFS,log=T)

  SFS.AIC <- 2*2*(length(SFS.res$knots)-1)-2*log.prob.data

  if(!is.null(boot.res)){
    
    if(length(which(unlist(lapply(boot.res,function(x) identical(x$knots,SFS.res$knots)))==F))!=0){
      stop("The time points of the bootstrap samples must equal that of the observed SFS")
    } else {
      
      tt <- seq(0,max(SFS.res$knots),len=nb.time)
      
      Lambda <-lapply(boot.res,function(boot){
        if(boot$ResSpline$convergence!=3 &
           all( boot$ResSpline$par-diff(boot$knots) < 10^(-5) )){
          return(NULL)
        } else {
            Lambda <- SplineFct(tt,cumsum(boot$ResSpline$par),boot$knots)$spline$value
            return(Lambda)
        }
      })
      
      LambdaMat <- matrix(unlist(Lambda),nrow=nb.time)
      
      med.Lam <- as.vector(apply(LambdaMat,1,function(x)quantile(x,probs=0.5)))
      
      slp.prm <- diff(med.Lam)/diff(tt) 
    
      x1 <- 1/outer(choose(2:n.samples,2),slp.prm)
      x2 <- exp(-outer(choose(2:n.samples,2),med.Lam[-length(med.Lam)]))
      x3 <- exp(-outer(choose(2:n.samples,2),med.Lam[-length(med.Lam)]+slp.prm*diff(tt)))
      evecMat <- x1*(x2-x3)
    
      ## the constant population size after the knots
    
      b_n <- median(unlist(lapply(boot.res,function(boot) {
                                        if(boot$ResSpline$convergence!=3 &
                                           all( boot$ResSpline$par-diff(boot$knots) < 10^(-5) )){
                                          return(NULL)
                                        } else {
                                          bb_n <- SplineFct(tt,cumsum(boot$ResSpline$par),boot$knots)$bb.spl[length(boot$knots)]
                                          return(bb_n)
                                        }
                                          
      })))
          
      x4 <- exp(-choose(2:n.samples,2)*med.Lam[length(tt)])/(choose(2:n.samples,2)*b_n) ##med.Lam[length(tt)]=a_m
      evec.med <- rowSums(evecMat)+x4
    
      prob.SFS <- (Wmat%*%evec.med/sum(Vmat*evec.med))
      if(is.folded){
        UFprobSFS <- probSFS
        if(length(UFprobSFS) %% 2 == 0){
          prob.SFS <- (UFprobSFS +rev(UFprobSFS))[1:(length(UFprobSFS)/2)]
        } else{
          prob.SFS <- (UFprobSFS +rev(UFprobSFS))[1:(floor(length(UFprobSFS)/2))]
          prob.SFS <- c(probSFS,UFprobSFS[(length(UFprobSFS)+1)/2])
        }
      }
      if(is.folded){
        log.prob.data <- dmultinom(SFS.res$SFS,prob=probSFS,log=T)
      } else {
        log.prob.data <- dmultinom(SFS.res$SFS,prob=prob.SFS,log=T)
      }  
      boot.AIC <- 2*2*(length(SFS.res$knots)-1)-2*log.prob.data
      
      res <- c(SFS.AIC,boot.AIC)
      names(res) <- c("AIC from SFS","AIC from Bootstrap")
    }
  } else {
    res <- c(SFS.AIC)
    names(res) <- c("AIC from SFS")
  }
  return(res)

}
