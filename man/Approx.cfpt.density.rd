\encoding{latin1}
\name{Approx.cfpt.density}
\alias{Approx.cfpt.density}
\title{Approximating First-Passage-Time Densities for Conditioned Problems}
\description{
  \code{Approx.cfpt.density} computes values of the approximate first-passage-time (f.p.t.) density, for a conditioned problem, 
  from an object of class \dQuote{summary.fptl} that contains the information provided by First-Passage-Time Location (FPTL) function.
}
\usage{
Approx.cfpt.density(sfptl, variableStep = TRUE, from.t0 = FALSE, 
                   to.T = FALSE, skip = TRUE, n = 250, p = 0.2, 
                   alpha = 1, tol = 1e-03, it.max = 50000L)                      
}
\arguments{
  \item{sfptl}{an object of class \dQuote{summary.fptl}, a result of applying the \code{\link{summary.fptl}} method to an object
  of class \dQuote{fptl}.}
  \item{variableStep}{a logical value indicating whether a variable integration step is used.}
  \item{from.t0}{a logical value indicating whether the approximation should be calculated from the lower end of the
interval considered, \eqn{t_0}{t0}, specified in the \code{sfptl} object.}
  \item{to.T}{a logical value indicating whether the approximation should be calculated to the upper end of the
interval considered, \eqn{T}, specified in the \code{sfptl} object.}
  \item{skip}{a logical value indicating whether the intervals at which the FPTL function is near zero could be
avoided.}
  \item{n}{Number of points used to determine the integration step in subintervals \eqn{[t_i^*, t_{max,i}^+]}{[t[i]*, tmax[i]^+]} 
\eqn{i=1, \ldots, m}{i=1,..., m}, from interesting instants provided by the FPTL function.}
	\item{p}{Ratio of n used to determine the integration step in subintervals  
\eqn{[t_{max,i}^+, t_{i+1}^*]}{[tmax[i]^+, t[i+1]*]}, \eqn{i=1, \ldots, m}{i=1,..., m}, \eqn{[t_0, t_1^*]}{[t0, t[1]*]} and 
\eqn{[t_{max,m}^{+}, T]}{[tmax[m]^+, T]}.}
  \item{alpha}{Parameter used to determine the integration step in subintervals 
\eqn{[t_{max,i}^+, t_{i+1}^*]}{[tmax[i]^+, t[i+1]*]}, \eqn{i=1, \ldots, m}{i=1,..., m}, \eqn{[t_0, t_1^*]}{[t0, t[1]*]} and 
\eqn{[t_{max,m}^{+}, T]}{[tmax[m]^+, T]}, in order to reduce the computational cost of approximating the f.p.t. density function 
in those cases where \eqn{t_{i+1}^* - t_{max,i}^+ >> t_{max,i}^{-} - \thinspace t_i^*}{t[i+1]* - tmax[i]^+ >> tmax[i]^- - t[i]*}, 
for some \eqn{i}, \eqn{t_1^* - t_0 >> t_{max,1}^{-} - \thinspace t_1^*}{t[1]* - t0 >> tmax[1]^- - t[1]*} or 
\eqn{T - t_{max,m}^+ >> t_{max,m}^{-} - \thinspace t_m^*}{T - tmax[m]^+ >> tmax[m]^- - t[m]*}, respectively.}  
  \item{tol}{If the cumulative integral of the approximation is greater than or equal to 1 - tol the algorithm is stopped.}
  \item{it.max}{If the number of iterations required for the approximation process is greater than it.max, the function asks for permission to continue.}
}
\details{
For a diffusion process \eqn{\{X(t), t_0 \leq t \leq T \}}{{X(t), t0 \le t \le T}}, the f.p.t. variable,   
conditioned to \eqn{X(t_0) = x_0}{X(t0) = x0}, through a continuous boundary \eqn{S(t)} is defined as
\ifelse{latex}{\deqn{T_{S(t), x_0} = \left\{
\begin{array}{lll}
\mbox{Inf} \ \{ t \geq t_0 \ : \ X(t) > S(t) \negthinspace \mid \negthinspace X(t_0)=x_0 \} & & \mbox{if \ } x_0 < S(t_0) \vspace{5pt} \\
\mbox{Inf} \ \{ t \geq t_0 \ : \ X(t) < S(t) \negthinspace \mid \negthinspace X(t_0)=x_0 \} & & \mbox{if \ } x_0 > S(t_0)
\end{array}
 \right.}}{\deqn{T = Inf{ t \ge t0 : X(t) > S(t) | X(t0) = x0 } ,} 
if \eqn{x_0 < S(t_0)}{x0 < S(t0)}, and
\deqn{T = Inf{ t \ge t0 : X(t) < S(t) | X(t0) = x0 } ,} 
if \eqn{x_0 > S(t_0)}{x0 > S(t0)}. }

Its density function is the solution to a Volterra integral equation of the second kind. The kernel of this equation depends
on the infinitesimal moments of the process, the transition probability density function and the boundary. 

Nevertheless, and apart from some particular processes and boundaries,
closed-form solutions for the integral equation are not available. For this reason, 
in the cases without explicit solutions, numerical procedures are required. That is the situation 
considered here and the numerical procedure implemented by the \code{Approx.fpt.density} function is the one 
proposed by Buonocore et al. (1987), based on the composite trapezoid method. \cr

The \code{Approx.cfpt.density} function computes efficiently the approximate f.p.t. density by
using the information provided by the FPTL function contained in the 
\code{sfptl} object. See the function \code{\link{summary.fptl}} for details. \cr

By default the function does not compute the approximate f.p.t. density 
from the time instant \eqn{t_0}{t0}, but from a more suitable time instant \eqn{t_1^*}{t[1]*}
provided by the FPTL function. It also uses a variable integration step. \cr

The function makes an internal call to \code{\link{Integration.Steps}} function in order to determine the subintervals and 
integration steps to be used in the application of the numerical algorithm according to the \code{variableStep}, 
\code{from.t0}, \code{to.T}, \code{n}, \code{p} and \code{alpha} arguments. \cr

In addition, if \code{skip = TRUE}, the function checks the approximate density value for each \eqn{t_{max,i}^+}{tmax[i]^+}, 
and, if it is almost 0, the application of the numerical algorithm in the subinterval \eqn{[t_{max,i}^+, t_{i+1}^*]}{[tmax[i]^+, t[i+1]*]} 
is avoided, and then continued from instant \eqn{t_{i+1}^*}{t[i+1]*} considering a zero value of the approximate 
density. \cr

Similarly, if \code{to.T = FALSE}, the function checks the cumulative value of the integral for each \eqn{t_{max,i}^+}{tmax[i]^+} 
provided by the FPTL function and, if it is greater than or equal to 1 - tol, the numerical algorithm is stopped. In any case, 
the algorithm is stopped in the final \eqn{t_{max,i}^+}{tmax[i]^+}, and if the cumulative value of the integral is less than 1 - tol
the function issues a warning. \cr
}
\value{
The \code{Approx.cfpt.density} function computes and returns an object of class \dQuote{fpt.density}. It 
is a three-component list:
  \item{x}{a sequence of suitable time instants in \eqn{[t_0, \ T]}{[t0, T]} according to the arguments in the function call.}
  \item{y}{the approximate conditioned f.p.t. density function values on the x sequence.} \cr
  \item{y.x0}{NULL (for consistency with the object of class \dQuote{fpt.density} that produces the \code{Approx.fpt.density} function).}
  
It also includes six additional attributes: 
\tabular{rl}{
\code{Call} \tab the unevaluated function call, substituting each name in this call by its value when \cr 
\code{} \tab the latter has length 1. \cr
\code{Steps} \tab matrix of subintervals and integration steps to consider for computing \cr 
\code{} \tab the approximate conditioned f.p.t. density. \cr 
\code{cumIntegral} \tab vector of the values of the cumulative integral of the \cr
\code{} \tab approximation for each subinterval considered. \cr
\code{skips} \tab a list that contains, for each subinterval, the value 1 if the application of the \cr
\code{} \tab numerical algorithm has been avoided or integer(0) otherwise. \cr
\code{CPUTime} \tab matrix of user and system times, by columns, required to approximate \cr
\code{} \tab the density for each subinterval considered, by rows. \cr
\code{summary.fptl} \tab the object used as \code{sfptl} argument in the function call. \cr
}                                                                                          

\code{x} is the vector result of the concatenation of the sequences of equally spaced values in the suitable subintervals 
determined by the \code{\link{Integration.Steps}} function. 
}
\references{
Buonocore, A., Nobile, A.G. and Ricciardi, L.M. (1987) A new integral equation for the evaluation of
first-passage-time probability densities. \emph{Adv. Appl. Probab.}, \bold{19}, 784--800.

Román, P., Serrano, J. J., Torres, F. (2008) First-passage-time location function: Application to determine
first-passage-time densities in diffusion processes. \emph{Comput. Stat. Data Anal.}, \bold{52}, 4132--4146.

P. Román-Román, J.J. Serrano-Pérez, F. Torres-Ruiz. (2012) An R package for an efficient approximation of first-passage-time 
densities for diffusion processes based on the FPTL function. \emph{Applied Mathematics and Computation}, \bold{218}, 8408--8428.
}
\author{Patricia Román-Román, Juan J. Serrano-Pérez and Francisco Torres-Ruiz.}
\seealso{
\code{\link{summary.fptl}} to locate the f.p.t. variable and create objects of class \dQuote{summary.fptl}.

\code{\link{is.fpt.density}} to test for objects of class \dQuote{fpt.density}.

\code{\link{print.fpt.density}} to show objects of class \dQuote{fpt.density}. 

\code{\link{report.fpt.density}} to generate a report.

\code{\link{plot.fpt.density}} for graphical display.

\code{\link{FPTL}} to evaluate the FPTL function and create objects of class \dQuote{fptl}.
}
\examples{
## Continuing the summary.fptl(.) example:
\dontshow{Lognormal <- diffproc(c("m*x","sigma^2*x^2","dnorm((log(x)-(log(y)+(m-sigma^2/2)*(t-s)))/(sigma*sqrt(t-s)),0,1)/(sigma*sqrt(t-s)*x)", "plnorm(x,log(y)+(m-sigma^2/2)*(t-s),sigma*sqrt(t-s))")) ; 
b <- "4.5 + 4*t^2 + 7*t*sqrt(t)*sin(6*sqrt(t))" ; y <- FPTL(dp = Lognormal, t0 = 0, T = 18, x0 = 1, S = b, list(m = 0.48, sigma = 0.07)) ; yy <- summary(y);
LognormalFEx <- diffproc(c("`h(t)`*x", "sigma^2*x^2", "dnorm((log(x)-(log(y)+`H(s,t)`-(sigma^2/2)*(t - s)))/(sigma*sqrt(t-s)),0,1)/(sigma*sqrt(t-s)*x)", "plnorm(x,log(y)+ `H(s,t)`-(sigma^2/2)*(t-s),sigma*sqrt(t-s))"));
z <- FPTL(dp = LognormalFEx, t0 = 1, T = 10, x0 = 1, S = 15, env = list(sigma=0.1, `h(t)` = "t/4", `H(s,t)` = "(t^2-s^2)/8")); zz <- summary(z)}
## Making an efficient approximation of the f.p.t. density 
## (optimal variable integration steps and small computational cost)
yyy <- Approx.cfpt.density(yy)
yyy

zzz <- Approx.cfpt.density(zz)
zzz

## Making a less efficient approximation of the f.p.t. density 
## (optimal fixed integration step but high computational cost related to 
##  the efficient approximation)
\dontrun{
yyy1 <- Approx.cfpt.density(yy, variableStep = FALSE, from.t0 = TRUE, to.T = 
                         TRUE, skip = FALSE)
yyy1}
}
\keyword{classes}
\keyword{list}
\keyword{methods}
\keyword{print}
