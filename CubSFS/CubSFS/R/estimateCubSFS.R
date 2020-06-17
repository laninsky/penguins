
#' Non-parametric estimation of population size changes from the site frequency spectrum (SFS)
#'
#' Based on the SFS, a cubic spline, imposed on the cumulative intensity of coalescence, is describing the changes of population size back in time.
#' Given some time points (in coalescent units), the cubic spline is determined my minizing the score function
#' \deqn{S(\Lambda)= (1-\alpha)\sum\limits_{i=1}^{n-1}\frac{\left(\E[\xi_i]-\xi_i\right)^2}{\E[\xi_i]}+\alpha\int\limits_{0}^{\infty}\left(\Lambda''(t)\right)^2dt.}{S(\Lambda) = (1-\alpha) * \sum_{i=1}^{n-1} ( E[\xi_i] - \xi_i )^2 / E[\xi_i] + \alpha * integral_0^\infty ( \Lambda''(t) )^2 dt.}
#' Here \eqn{\Lambda(t)} is the cumulative intensity a time \eqn{t} (in coalescent units) of coalescence, \eqn{\xi} is the observed SFS, \eqn{E[\xi]} is the expected SFS 
#' calculated based on the estimated \eqn{\Lambda}, \eqn{\alpha} is a smoothing parameter, and 
#' \eqn{\int_{0}^{\infty} (\Lambda''(t))^2 dt}{integral_0^\infty (\Lambda''(t))^2 dt} is the penalising term.
#'  
#' The probability of observing a site with \eqn{i} derived alleles is evaluated using the formula by Polanski and Kimmel (2003)
#' \deqn{p_i(\Lambda)=\frac{\sum\limits_{j=2}^ne_j (\Lambda)W_{ij}}{\sum\limits_{j=2}^n e_j(\Lambda)V_j}}{p_i(\Lambda) = ( \sum_{j=2}^n e_j(\Lambda) W_{ij} ) / ( \sum_{j=2}^n e_j(\Lambda) V_j )}
#' for \eqn{i=1,2,\ldots,n-1} for unfolded SFS (\code{is.folded=F}). For the folded SFS the probability of observing a site with minor 
#' allele count \eqn{i} among \eqn{n} sequences is 
#' \deqn{p_i^F(\Lambda) = p_i(\Lambda) + p_{n-i}(\Lambda)}
#' for \eqn{i=1,\ldots,n/2} if \eqn{n} is even and for \eqn{n} odd
#' \deqn{p_i^F(\Lambda) = p_i(\Lambda) + p_{n-i}(\Lambda)}
#' for \eqn{i=1,\ldots,(n-1)/2} and \eqn{p_{(n+1)/2}^F(\Lambda) = p_{(n+1)/2}(\Lambda)}
#' 
#' Here \eqn{\Lambda} is an increasing cubic spline, 
#' \deqn{\Lambda(t)=a_i+b_i(t-t_i)+c_i(t-t_i)^2+d_i(t-t_i)^3}{\Lambda(t)=a_i + b_i*(t-t_i) + c_i*(t-t_i)^2 + d_i*(t-t_i)^3} 
#' for \eqn{t_{i-1} \le t < t_i}, such that \eqn{\Lambda(0)=0}, \eqn{\Lambda'(0)=\lambda(0)=1}, and 
#' \eqn{\Lambda(t)=a_m+b_m(t-t_m)} for \eqn{t_m \le t}, where \eqn{\lambda(t) = \Lambda'(t)} is the coealecent rate.
#' 
#' Given an estimated cubic spline, 
#' \deqn{e_j(\Lambda)=\int_0^{\infty}exp(-{j\choose 2}\Lambda(t))dt}{e_j(\Lambda) = intergral_0^\infty exp( -choose(j,2) \Lambda(t) ) dt} 
#' is evaluated by assuming piece wise straight lines between any two points of a fine partitioning of the \code{knots}. The number of time points 
#' in the fine partitioning is \code{nb.time}.
#'    
#' \code{alpha} is estimated using cross validation with \code{nb.groups} equaly sized groups of sites. 
#' Each group represent a validation set to be used exactly ones. The cross validation procedure is implemented in \code{\link{estimateAlpha}}.
#' 
#' The time points (\code{knots}), on which to base the cubic spline, must be in coalescent units. The first time point must be 0. If not specified a total of \code{n.knots} 
#' (default=15) time points will be destributed on the interval from 0 to \code{t_m} (default=2) following
#' \deqn{t_i = \frac{1}{\kappa}\left\{ \exp\left[ \frac{i}{m-1}\log(1+\kappa t_m )\right] -1\right\}}{t_i = 1/\kappa ( exp[ i/(m-1) log(1+\kappa t_m ) ] -1 )}
#' for \eqn{i=0,1,\ldots,m-1}. Here \eqn{\kappa} (default=1) determines the clustering of points near 0.
#' 
#' The minimization of the score function w.r.t. \eqn{0<\phi_i=a_{i}-a_{i-1}} for \eqn{i=1,2,\ldots,m} is performed using the augmented 
#' Lagrangian algorithm implemented in the \code{nloptr} package, to ensure the property of \eqn{\Lambda} to be increasing.
#' 
#' Given the results from \code{estimateCubSFS}, say \code{res}, the call \code{SplineFct(tt,cumsum(res$par),res$konts)} will provide
#' the information on the spline, such as the value \eqn{\Lambda(t)} in the time points given by \code{tt}, and the estimate of the \eqn{b_i}'s, 
#' the \eqn{c_i}'s and the \eqn{d_i}'s. The coalenscent rate \eqn{\lambda(t)} is then evaluated using
#' \deqn{\lambda(t) = \Lambda'(t) = b_i + 2c_i(t-t_i) + 3d_i(t-t_i)^2.}{\lambda(t) = \Lambda'(t) = b_i + 2*c_i*(t-t_i) + 3*d_i*(t-t_i)^2.} 
#'
#' Given the present poulation size \eqn{N(0)} the population size back in generation \eqn{r=2N(0) t}{r=2N(0)*t} is 
#' \deqn{2N(r)=\frac{2N(0)}{\lambda(r/(2N(0)))}.}{2N(r) = 2N(0) / \lambda( r/(2N(0)) ).}
#' 
#' The function PlotCubSFS provides plots of the changes in population size, given a sequence of time points measured in years back in time.
#' 
#' @param rawSFS the observed SFS. 
#' @param n.samples the number of sequences used to generate the observed SFS.
#' @param knots defined time points \eqn{t_i} for \eqn{i=0,1,\ldots,m} in coalescent units. If NULL(default) \code{n.knots} and \code{t_m} are used to create the knots.
#' @param n.knots the total number of time points (inkluding 0). 
#'   If \code{knots} are specificed \code{n.knots}=\code{length(knots)}.
#' @param t_m the value of the last time point in coalescent units. 
#'   If \code{knots} are specificed \code{t_m}=\code{knots[n.knots]}.
#' @param kappa if \code{knots} is NULL, then \code{kappa} determines the clustering of the time points close to 0.
#' @param nb.groups the number of groups used in the cross validation estimation of the smoothing parameter.
#' @param is.folded is the SFS folded(TRUE) or unfolded (FALSE)
#' @param loc.solv the local solver. May be either \code{COBYLA} or \code{LBFGS}.
#' @inheritParams nloptr::auglag
#' 
#' @return A list: The results from \code{\link[nloptr]{auglag}}. inkluding the observed SFS, the optimal value of alpha and 
#'   the value of the knots.
#'   \describe{
#'     \item{ResSpline}{the result from \code{\link[nloptr]{auglag}}:
#'       \describe{
#'         \item{par}{The estimated parameter \eqn{\phi_i=a_{i}-a_{i-1}}. The optimal solution found so far.}
#'         \item{value}{The value of the score function corresponding to \code{par}.}
#'         \item{iter}{the number of iterations used.}
#'         \item{global_solver}{The global NLOPT solver used.}
#'         \item{local_solver}{The local NLOPT solver used.}
#'         \item{convergence}{An interger code indicating successfull completion (3), running out of iterations (5), 
#'                              other successfull completion (>0), or error number (<0).}
#'         \item{message}{Character string produced by NLOPTR and giving additional information on \code{convergence}.}
#'      }}   
#'    \item{SFS}{the observed SFS}
#'    \item{alpha}{the result from \code{\link[stats]{optimize}}}
#'      \describe{
#'        \item{minimum}{the estimate of the smoothing parameter \eqn{\alpha}.}
#'        \item{objective}{the value of the mean square error used to estimate the smoothing parameter.}
#'      }
#'    \item{knots}{the time points used for estimation \eqn{t_0,t_1,\ldots,t_m} in coalescent units}
#'    \item{expSFS}{the expected SFS under the estimated model.}
#'    \item{CoalRate}{a dataframe with three columns. time (\code{CoalTime}), coalescent rate (\code{CoalRate}) and 
#'    the integrated coalescent rate (\code{IntCoalRate}).}
#'    \item{CoalRatePlot}{a plot of the coalescent rate given in \code{CoalRate}}
#'    \item{SFSPlot}{the plot of the observed SFS and the expected SFS under the estimated model.}
#'    \item{qqPlot}{a plot of the expected SFS against the observed SFS.}
#'  }
#'    
#'    
#' @seealso \code{\link{estimateAlpha}} \code{\link[nloptr]{auglag}} \code{\link{SplineFct}} \code{\link{AicCubSFS}} \code{\link{PlotCubSFS}}
#' 
#' @references 
#' Polanski, A. and Kimmel, M. (2003), "New Explicit Expressions for Relative Frequencies of Single-Nucleotide Polymorphisms With Application 
#' to Statistical Inference on Population Growth." Genetics, 165, 427-436.
#' 
#' @export

estimateCubSFS <-
function(rawSFS,n.samples,knots=NULL,
                   n.knots=ifelse(is.null(knots),15,length(knots)),
                   t_m=ifelse(is.null(knots),2,knots[length(knots)]),
                   alpha=NULL,kappa=1,nb.groups=5,is.folded=FALSE, 
                   nb.time=500, tol=10^(-2),
                   conv.lim=10^{-3},
                   loc.solv="COBYLA",
                   control=list(maxeval=20000, ftol_abs=10^(-5),stopval=0)){


  if(n.knots<4){stop("Number of time points must be more than 3")}
  if(nb.groups<3){stop("Number of groups used for cross validation must me more than 2")}

  print("Doing Wmat step now")
  Wmat <- matrix( unlist( mclapply(1:(n.samples-1),WijFunction,n=n.samples,mc.cores=numcores) ),byrow=TRUE,nrow=(n.samples-1))
  print("Doing Vmat step now")
  Vmat <- as.vector(sapply(2:n.samples,VjFunction,n=n.samples))

  print("Doing knots step now")
  if(is.null(knots)){
    knots <- 1/kappa*( exp( (0:(n.knots-1))/(n.knots-1)*log(1+kappa*t_m) ) -1) ## in coalescent units
  } else {
    if(knots[1]!=0){stop("The first time point must be 0")}
  }

  if(knots[n.knots]>10000){warning("time points must be in coalescent units")}

  theta.prm <- diff(knots)

  print("Doing alpha step now")
  if(is.null(alpha)){
    alpha <- estimateAlpha(rawSFS,n.samples,knots,theta.prm,nb.groups,Wmat,Vmat,is.folded,
                         nb.time=nb.time,tol=tol,loc.solv=loc.solv,control=control)
  } else {
    alphaL <- list()
    alphaL$minimum <- alpha
    alphaL$objective <- NA
    alpha <- alphaL
  }
  
  res <- MinScore(rawSFS,n.samples,knots,theta.prm,alpha$minimum,Wmat,Vmat,is.folded,
                  nb.time=nb.time,conv.lim=conv.lim,loc.solv=loc.solv,control=control)

  evec.spl <- EjFct(cumsum(res$par),knots,n.samples)
  
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
  
  expSFS <- sum(rawSFS)*probSFS
  
  tt <- seq(0,max(knots),len=1000)
  
  Lambda <- SplineFct(tt,cumsum(res$par),knots)
  
  ii <- findInterval(tt,knots)
  lambda <- Lambda$bb.spl[ii]+2*Lambda$cc.spl[ii]*(tt-knots[ii])+
    3*Lambda$dd.spl[ii]*(tt-knots[ii])^2
  
  resMat <- data.frame(CoalTime=tt,CoalRate=lambda,IntCoalRate=Lambda$spline$value)
  
  plot(tt,lambda,type="l",xlab="time in coalescent units",ylab="coalescent rate")
  plot.lambda <- recordPlot()
  
  plot(1:length(rawSFS),rawSFS,type="l",xlab="",ylab="SFS")
  points(1:length(expSFS),expSFS,col="red",pch=20)
  legend("topright",legend=c("observed SFS","expected SFS"),col=c("black","red"),
         lty=c(1,NA),pch=c(NA,20),bty="n")
  plot.SFS <- recordPlot()
  
  plot(expSFS,rawSFS,pch=20,xlab="expected SFS",ylab="observed SFS")
  abline(a=0,b=1,lty=2)
  plot.qq <- recordPlot()
  
  return(list(ResSpline=res,SFS=rawSFS,alpha=alpha,knots=knots,expSFS=expSFS,CoalRate=resMat,CoalRatePlot=plot.lambda,SFSPlot=plot.SFS,qqPlot=plot.qq))
}

