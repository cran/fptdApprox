\encoding{latin1}
\name{Integration.Steps}
\alias{Integration.Steps}
\title{Subintervals and Integration Steps To Approximate First-Passage-Time Densities}
\description{
  According to the First-Passage-Time Location (FPTL) function and the arguments in the function call, this function  
  calculates suitable subintervals and integration steps in order to approximate the first-passage-time (f.p.t.) density. 
}
\usage{
Integration.Steps(sfptl, variableStep = TRUE, from.t0 = FALSE,
                  to.T = FALSE, n = 250, p = 0.2)
}
\arguments{
  \item{sfptl}{an object of class \dQuote{summary.fptl}, a result of applying the \code{\link{summary.fptl}} method to 
  an object of class \dQuote{fptl}.}
  \item{variableStep}{a logical value indicating whether a variable integration step is used.}
  \item{from.t0}{a logical value indicating whether the approximation should be calculated from the lower end of the
interval considered, \eqn{t_0}{t0}, specified in the object used as \code{sfptl} argument.}
  \item{to.T}{a logical value indicating whether the approximation should be calculated to the upper end of the
interval considered, \eqn{T}, specified in the object used as \code{sfptl} argument.}
	\item{n}{Number of points used to determine the integration step in subintervals where the f.p.t. density is not
constant according to the FPTL function, that is, in \eqn{[t_i^*, t_{max,i}^+]}{[t[i]*, tmax[i]^+]} 
from interesting instants provided by the FPTL function and contained in the \code{sfptl} object.}
  \item{p}{Ratio of n used to determine the integration step in subintervals where the f.p.t. density function is 
near constant according to the FPTL function, that is, in \eqn{[t_{max,i}^+, t_{i+1}^*]}{[tmax[i]^+, t[i+1]*]}
from interesting instants provided by the FPTL function and contained in the \code{sfptl} object.}
}
\details{
Based on the information provided by the FPTL function contained in the \code{sfptl} object, this function computes and returns 
suitable subintervals and integration steps in order to approximate the density 
function of the f.p.t. variable according to the other arguments in the function call. \cr

If \code{variableStep = TRUE}, the function computes the variable integration step
\deqn{h_i = \frac{t_{max,i}^{-} - t_i^*}{n} \ ,}{h[i] = (tmax[i]^- - t[i]*)/n ,} in subintervals 
\eqn{[t_i^*, t_{max,i}^+]}{[t[i]*, tmax[i]^+]}, \eqn{i=1, \ldots, m}{i=1,..., m}, and 
\deqn{h_i^*= \frac{t_{i+1}^*-t_{max,i}^+}{[np+0.5]} \ ,}{h[i]* = (t[i+1]* - tmax[i]^+)/[np+0.5] ,} in subintervals 
\eqn{[t_{max,i}^+, t_{i+1}^*]}{[tmax[i]^+, t[i+1]*]}, \eqn{i=1, \ldots, m-1}{i=1,..., m-1}, where \eqn{[x]} is the integer part of \eqn{x}. \cr

If \code{variableStep = FALSE} the function determines the fixed integration step \deqn{h = min \left\{ h_i \ , \ i=1, \ldots, m \right\} .}{h = min {h[i] , i=1,..., m} .}

If \code{from.t0 = TRUE}, the function also computes a suitable integration step in subinterval \eqn{[t_0, t_1^*]}{[t0, t[1]*]},
once again depending on the \code{variable.Step} argument. If \code{variableStep = TRUE} it considers the integration step 
\deqn{h_{0} = \frac{t_1^* - t_0}{[np+0.5]} \ .}{h[0] = (t[1]* - t0)/[np+0.5] .}
Otherwise, the fixed integration step defined above is considered. \cr

If \code{to.T = TRUE}, the function computes a suitable integration step in subinterval \eqn{[t_{max,m}^{+}, T]}{[tmax[m]^+, T]},
once again depending on \code{variable.Step} argument. If \code{variableStep = TRUE} it considers the integration step 
\deqn{h_{m} = \frac{T - t_{max,m}^{+}}{[np+0.5]} \ .}{h[m] = (T - tmax[m]^+)/[np+0.5] .}
In other case, the fixed integration step defined above is considered. \cr
}
\value{
  A matrix of subintervals and integration steps that we must consider in order to approximate 
	the f.p.t. density according to the information contained in the \code{sfptl} object and 
	the arguments in the function call. A final column of binary values indicates the subintervals 
	in which we must check whether it is possible to avoid applying the algorithm.
}
\references{
Román, P., Serrano, J. J., Torres, F. (2008) First-passage-time location function: Application to determine
first-passage-time densities in diffusion processes. \emph{Comput. Stat. Data Anal.}, \bold{52}, 4132--4146.
}
\author{Patricia Román-Román, Juan J. Serrano-Pérez and Francisco Torres-Ruiz.}
\seealso{
\code{\link{Approx.fpt.density}} to approximate f.p.t. densities from objects of class \dQuote{summary.fptl} and create objects of class \dQuote{fpt.density}.

\code{\link{summary.fptl}} to locate the f.p.t. variable and create objects of class \dQuote{summary.fptl} from objects of class 
\dQuote{fptl}.

\code{\link{FPTL}} to evaluate the FPTL function and create objects of class \dQuote{fptl}.
}
\examples{
## Continuing the summary.fptl(.) example:
\dontshow{Lognormal <- diffproc(c("m*x","sigma^2*x^2","dnorm((log(x)-(log(y)+(m-sigma^2/2)*(t-s)))/(sigma*sqrt(t-s)),0,1)/(sigma*sqrt(t-s)*x)", "plnorm(x,log(y)+(m-sigma^2/2)*(t-s),sigma*sqrt(t-s))")) ; 
b <- "4.5 + 4*t^2 + 7*t*sqrt(t)*sin(6*sqrt(t))" ; y <- FPTL(dp = Lognormal, t0 = 0, T = 18, x0 = 1, S = b, list(m = 0.48, sigma = 0.07)) ; yy <- summary(y)}
Integration.Steps(yy)
Integration.Steps(yy, from.t0 = TRUE)
Integration.Steps(yy, to.T = TRUE, n = 100, p = 0.25)
}

