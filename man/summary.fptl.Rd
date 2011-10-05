\encoding{latin1}
\name{summary.fptl}
\alias{summary.fptl}
\alias{is.summary.fptl}
\alias{print.summary.fptl}
\title{Locating the First-Passage-Time Variable}
\description{
\code{summary.fptl} summary method for class \dQuote{fptl}.

\code{is.summary.fptl} tests if its argument is an object of class \dQuote{summary.fptl}.

\code{print.summary.fptl} shows an object of class \dQuote{summary.fptl}.
}
\usage{
\method{summary}{fptl}(object, zeroSlope = 0.01, p0.tol = 8, k = 3, \dots)

is.summary.fptl(obj)

\method{print}{summary.fptl}(x, \dots)
}
\arguments{
  \item{object}{an object of class \sQuote{fptl}, a result of a call to \code{\link{FPTL}}.}
  \item{obj}{an \R object to be tested.}
  \item{x}{an object of class \sQuote{summary.fptl}, a result of a call to \code{summary.fptl}.}
  \item{zeroSlope}{maximum slope required to consider that a growing function is constant.}
  \item{p0.tol}{controls where the First-Passage-Time Location function begins to increase significantly.}
  \item{k}{controls whether the First-Passage-Time Location function decreases very slowly.}
  \item{\dots}{other arguments passed to functions.}
}
\details{
The \code{summary.fptl} function extracts the information contained in \code{object} about the
location of the variation range of the first-passage-time (f.p.t.) variable. \cr

It makes an internal call to \code{\link{growth.intervals}} function in order to determine the time instants 
\eqn{t_i, \ i=1, \ldots, m}{t[i], i=1,...,m}, from which the First-Passage-Time Location (FPTL)
function starts growing, and its local maximums \eqn{t_{max,i}}{tmax[i]}. For this, \code{zeroSlope} argument
is considered. \cr

If there is no growth subinterval, the execution of the function \code{summary.fptl} is stopped and an error is reported.
Otherwise, for each of the subintervals \eqn{I_{i} = [t_{i},t_{i+1}]}{I[i] = [t[i], t[i+1]]} the function determines:
\itemize{
\item The first time instant \eqn{t_i^* \in [t_i, t_{max,i}]}{t[i]* in [t[i], tmax[i]]} at which the function is bigger than
or equal to \deqn{p_i^* = p_i + 10^{-p0.tol}(p_{max,i} - p_i) \ ,}{p[i]* = p[i] + 10^{-p0.tol}(pmax[i] - p[i]) ,}
where \eqn{p_i = FPTL(t_i)}{p[i] = FPTL(t[i])} and \eqn{p_{max,i} = FPTL(t_{max,i}) \ . }{pmax[i] = FPTL(tmax[i]) .} \cr \cr
\eqn{10^{-p0.tol}} is the ratio of the global increase of the function in the growth subinterval \eqn{[t_i, t_{max,i}]}{[t[i], tmax[i]]} that
should be reached to consider that it begins to increase significantly.
\item The first time instant \eqn{t_{max,i}^{-} \in [t_i, t_{max,i}]}{tmax[i]^- in [t[i], tmax[i]]} at which the FPTL function is
bigger than or equal to \deqn{p_{max,i}^{-} = p_{max,i} \big( 1 - 0.05(p_{max,i} - p_i) \big) \ .}{pmax[i]^- = pmax[i](1 - 0.05(pmax[i] - p[i])) .}
\item The last time instant \eqn{t_{max,i}^{+} \in \big[ t_{max,i}, \thinspace T_i \big]}{tmax[i]^+ in [tmax[i], T[i]]} at which the FPTL function is bigger
than or equal to \deqn{p_{max,i}^{+} = max \left\{ 1 - (1 - p_{max,i}^2)^{(1+q)/2}, FPTL(T_i) \right\},}{pmax[i]^+
= max { 1 - (1 - pmax[i]^2)^{(1+q)/2}, FPTL(T[i])} ,}
where \deqn{T_i = min \big\{ t_{max,i} + k \thinspace (t_{max,i}-t_i^*)(1 - p_{max,i}),
\thinspace t_{i+1} \big\}}{T[i] = min(tmax[i] + k (tmax[i] - t[i]*)(1 - pmax[i]), t[i+1])} and
\deqn{q = \displaystyle{\frac{p_{max,i} - p_i}{p_{max,i}}} \ . }{q = (pmax[i] - p[i])/pmax[i] .}
}

\code{print.summary.fptl} displays an object of class \dQuote{summary.fptl} for immediate understanding of the information it contains.
}
\value{
The \code{summary.fptl} function computes and returns an object of class \dQuote{summary.fptl}. It is a matrix whose
columns correspond to \eqn{t_i}{t[i]}, \eqn{t_i^*}{t[i]*}, \eqn{t_{max,i}^{-}}{tmax[i]^-}, \eqn{t_{max,i}}{tmax[i]}
and \eqn{t_{max,i}^{+}}{tmax[i]^+} values.
It also includes four additional attributes:
\item{Call}{the unevaluated function call.}
\item{FPTLValues}{matrix of values of the FPTL function for the interesting time instants in each subinterval \eqn{I_i}{I[i]}, by rows.
Its columns correspond to \eqn{FPTL(t_i)}{FPTL(t[i])}, \eqn{FPTL(t_i^*)}{FPTL(t[i]*)}, \eqn{FPTL(t_{max,i}^{-})}{FPTL(tmax[i]^-)},
\eqn{FPTL(t_{max,i})}{FPTL(tmax[i])} and \eqn{FPTL(t_{max,i}^{+})}{FPTL(tmax[i]^+)} values.}
\item{FPTLCall}{the unevaluated FPTL function call that resulted in the object used
as \code{object} argument in the function call.}
\item{dp}{the object used as \code{dp} argument in the FPTL function call that resulted in the \code{object} argument
in the function call. \cr \cr}

\code{is.summary.fptl} returns \code{TRUE} or \code{FALSE} depending on whether its argument is an object of class
\dQuote{summary.fptl} or not.}
\references{
Román, P., Serrano, J. J., Torres, F. (2008) First-passage-time location function: Application to determine
first-passage-time densities in diffusion processes. \emph{Comput. Stat. Data Anal.}, \bold{52}, 4132--4146.
}
\author{Patricia Román-Román, Juan J. Serrano-Pérez and Francisco Torres-Ruiz.}
\seealso{
\code{\link{Approx.fpt.density}} to approximate f.p.t. densities from objects of class \dQuote{summary.fptl} 
and create objects of class \dQuote{fpt.density}.

\code{\link{FPTL}} to evaluate the FPTL function and create objects of class \dQuote{fptl}.

\code{\link{report.summary.fptl}} to generate a report.

\code{\link{growth.intervals}} to study the growth of the vector of values resulting from the evaluation of a function.
}
\examples{
## Continuing the FPTL(.) example:
\dontshow{Lognormal <- diffproc(c("m*x","sigma^2*x^2","dnorm((log(x)-(log(y)+(m-sigma^2/2)*(t-s)))/(sigma*sqrt(t-s)),0,1)/(sigma*sqrt(t-s)*x)", "plnorm(x,log(y)+(m-sigma^2/2)*(t-s),sigma*sqrt(t-s))")) ;
b <- "4.5 + 4*t^2 + 7*t*sqrt(t)*sin(6*sqrt(t))" ; y <- FPTL(dp = Lognormal, t0 = 0, T = 18, x0 = 1, S = b, list(m = 0.48, sigma = 0.07))}
## Summarizing an object of class fptl
yy <- summary(y)
yy
print(yy, digits=10)
yy1 <- summary(y, zeroSlope = 0.001)
yy1
yy2 <- summary(y, zeroSlope = 0.001, p0.tol = 10)
yy2

## Testing summary.fptl objects
is.summary.fptl(yy)
}
\keyword{classes}
\keyword{array}
\keyword{methods}
\keyword{print}
