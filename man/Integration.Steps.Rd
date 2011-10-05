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
                  to.T = FALSE, n = 250, p = 0.2, alpha = 1)
}
\arguments{
  \item{sfptl}{an object of class \dQuote{summary.fptl}, a result of applying the \code{\link{summary.fptl}} method to 
  an object of class \dQuote{fptl}.}
  \item{variableStep}{a logical value indicating whether a variable integration step is used.}
  \item{from.t0}{a logical value indicating whether the approximation should be calculated from the lower end of the
interval considered, \eqn{t_0}{t0}, specified in the object used as \code{sfptl} argument.}
  \item{to.T}{a logical value indicating whether the approximation should be calculated to the upper end of the
interval considered, \eqn{T}, specified in the object used as \code{sfptl} argument.}
	\item{n}{Number of points used to determine the integration step in subintervals of \eqn{[t_0, \thinspace T]}{[t0, T]]}
where the f.p.t. density is not constant according to the FPTL function, that is, in \eqn{[t_i^*, t_{max,i}^+]}{[t[i]*, tmax[i]^+]}, 
\eqn{i=1, \ldots, m}{i=1,..., m}, from interesting instants provided by the FPTL function and contained in the \code{sfptl} object.}
  \item{p}{Ratio of n used to determine the integration step in subintervals of \eqn{[t_0, \thinspace T]}{[t0, T]]} 
where the f.p.t. density function is near constant according to the FPTL function, that is, in 
\eqn{[t_{max,i}^+, t_{i+1}^*]}{[tmax[i]^+, t[i+1]*]}, \eqn{i=1, \ldots, m}{i=1,..., m}, \eqn{[t_0, t_1^*]}{[t0, t[1]*]} and 
\eqn{[t_{max,m}^{+}, T]}{[tmax[m]^+, T]}.}
  \item{alpha}{Parameter used to determine the integration step in subintervals of \eqn{[t_0, \thinspace T]}{[t0, T]]}
where the f.p.t. density function is near constant according to the FPTL function, that is, in 
\eqn{[t_{max,i}^+, t_{i+1}^*]}{[tmax[i]^+, t[i+1]*]}, \eqn{i=1, \ldots, m}{i=1,..., m}, \eqn{[t_0, t_1^*]}{[t0, t[1]*]} and 
\eqn{[t_{max,m}^{+}, T]}{[tmax[m]^+, T]}, in order to reduce the computational cost of approximating the f.p.t. density function 
in those cases where \eqn{t_{i+1}^* - t_{max,i}^+ >> t_{max,i}^{-} - \thinspace t_i^*}{t[i+1]* - tmax[i]^+ >> tmax[i]^- - t[i]*}, 
for some \eqn{i}, \eqn{t_1^* - t_0 >> t_{max,1}^{-} - \thinspace t_1^*}{t[1]* - t0 >> tmax[1]^- - t[1]*} or 
\eqn{T - t_{max,m}^+ >> t_{max,m}^{-} - \thinspace t_m^*}{T - tmax[m]^+ >> tmax[m]^- - t[m]*}, respectively.}
}
\details{
Based on the information provided by the FPTL function contained in the \code{sfptl} object, this function computes and returns 
suitable subintervals and integration steps in order to approximate the density 
function of the f.p.t. variable according to the other arguments in the function call. \cr

In what follows, \eqn{[x]} is the integer part of \eqn{x}.

\itemize{
\item If \code{variableStep = TRUE}, the function computes the integration steps 
\deqn{h_i = \frac{t_{max,i}^{-} - \thinspace t_i^*}{n} \ ,}{h[i] = (tmax[i]^- - t[i]*)/n ,} in subintervals 
\eqn{[t_i^*, t_{max,i}^+]}{[t[i]*, tmax[i]^+]}, \eqn{i=1, \ldots, m}{i=1,..., m}. \cr

Then, the upper endpoints of subintervals \eqn{[t_i^*, t_{max,i}^+]}{[t[i]*, tmax[i]^+]} are readjusted according 
to the integration steps \eqn{h_i}{h[i]}, \eqn{i=1, \ldots, m}{i=1,..., m}. \cr

Next, in subintervals \eqn{[t_{max,i}^+, t_{i+1}^*]}{[tmax[i]^+, t[i+1]*]}, \eqn{i=1, \ldots, m-1}{i=1,..., m-1}, 
the function determines the integration steps
\deqn{h_i^* = \frac{t_{i+1}^* - t_{max,i}^+}{max \{ 1, \thinspace n_i \}} \ ,}{h[i]* = (t[i+1]* - tmax[i]^+)/max{1, n[i]}}
where
\ifelse{latex}{\deqn{n_i = \left\{ \begin{array}{ccc}
          \left[ n \thinspace p \thinspace k_{i} \right] &  \mbox{if} & k_i \leq 1 \vspace{0.25cm} \\
          \left[ n \thinspace p \thinspace k_{i}^{alpha} \right] & \mbox{if} & k_i > 1
        \end{array}
\right.} and}{\deqn{n[i] = [n p k[i]] ,} if \eqn{k[i] \le 1}, and \deqn{n[i] = [n p k[i]^alpha] ,} if \eqn{k[i] > 1}, with}
\deqn{k_i = \frac{t_{i+1}^* - \thinspace t_{max,i}^+}{t_{max,i}^{-} - \thinspace t_i^*} \ .}{k[i] = (t[i+1]* - tmax[i]^+)/(tmax[i]^- - t[i]*) .} \cr
 
\item If \code{variableStep = FALSE} the function determines the fixed integration step 
\deqn{h = min \left\{ h_i \ , \ i=1, \ldots, m \right\} \thinspace ,}{h = min {h[i] , i=1,..., m}, } 
and the endpoints of subintervals \eqn{[t_i^*, t_{max,i}^+]}{[t[i]*, tmax[i]^+]} are readjusted according 
to this integration step. \cr \cr \cr
} 

If \code{from.t0 = TRUE}, the function also computes a suitable integration step in subinterval \eqn{[t_0, t_1^*]}{[t0, t[1]*]},
once again depending on the \code{variable.Step} argument. 
\itemize{
\item If \code{variableStep = TRUE}, it considers the integration step 
\deqn{h_{0} = \frac{t_1^* - t_0}{max \{ 1, n_0 \} } \ ,}{h0 = (t[1]* - t0)/max{1, n0} , }
where
\ifelse{latex}{\deqn{n_0 = \left\{ \begin{array}{ccc}
          \left[ n \thinspace p \thinspace k_{0} \right] &  \mbox{if} & k_0 \leq 1 \vspace{0.25cm} \\
          \left[ n \thinspace p \thinspace k_{0}^{alpha} \right] & \mbox{if} & k_0 > 1
        \end{array}
 \right.} and}{\deqn{n0 = [n p k0] ,} if \eqn{k0 \le 1}, and \deqn{n[i] = [n p k0^alpha] ,} if \eqn{k0 > 1}, with}
\deqn{k_0 = \frac{t_1^* - \thinspace t_0}{t_{max,1}^{-} - \thinspace t_1^*} \ .}{k0 = (t[1]* - t0)/(tmax[1]^- - t[1]*) .} \cr

\item If \code{variableStep = FALSE}, the fixed integration step \eqn{h} defined above is considered. \cr \cr \cr
}

And, if \code{to.T = TRUE}, the function computes a suitable integration step in subinterval \eqn{[t_{max,m}^{+}, T]}{[tmax[m]^+, T]},
once again depending on \code{variable.Step} argument. 
\itemize{
\item If \code{variableStep = TRUE} it considers the integration step 
\deqn{h_m = \frac{T - \thinspace t_{max,m}^+}{max \{ 1, n_m \} } \ ,}{h[m] = (T - tmax[m]^+)/max{1, n[m]} , }
where
\ifelse{latex}{\deqn{n_m = \left\{ \begin{array}{ccc}
          \left[ n \thinspace p \thinspace k_{m} \right] &  \mbox{if} & k_m \leq 1 \vspace{0.25cm} \\
          \left[ n \thinspace p \thinspace k_{m}^{alpha} \right] & \mbox{if} & k_m > 1
        \end{array}
 \right.} and}{\deqn{n[m] = [n p k[m]] ,} if \eqn{k[m] \le 1}, and \deqn{n[m] = [n p k[m]^alpha] ,} if \eqn{k[m] > 1}, with}
\deqn{k_m = \frac{T - \thinspace t_{max,m}^+}{t_{max,m}^{-} - \thinspace t_m^*} \ .}{k[m] = (T - tmax[m]^+)/(tmax[m]^{-} - t[m]^*) .} \cr

\item If \code{variableStep = FALSE}, the fixed integration step \eqn{h} defined above is considered. \cr
}

In all cases, \eqn{p} and \eqn{alpha} are ignored if \code{variableStep = FALSE}. Otherwise, 
\eqn{p \geq 0.1}{p \ge 0.1} and \eqn{0.75 \leq alpha \leq 1}{0.75 \le alpha \le 1} are recommended; 
if not, some integration steps can be excessively large.
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
