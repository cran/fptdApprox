\encoding{latin1}
\name{FPTL}
\alias{FPTL}
\alias{is.fptl}
\alias{print.fptl}
\title{First-Passage-Time Location Function}
\description{
\code{FPTL} computes values of the First-Passage-Time Location (FPTL) function of a diffusion process for a continuous boundary.

\code{is.fptl} tests if its argument is an object of class \dQuote{fptl}.

\code{print.fptl} shows an object of class \dQuote{fptl}.
}
\usage{
FPTL(dp, t0, T, x0, S, env = NULL, n = 10000)

is.fptl(obj)

\method{print}{fptl}(x, \dots)
}
\arguments{
  \item{dp}{an object of class \dQuote{diffproc} defining a family of diffusion processes.}
  \item{obj}{an \R object to be tested.}
  \item{x}{an object of class \dQuote{fptl}, a result of a call to this function.}
  \item{t0, T}{lower and upper limits of the considered time interval. Must be finite.}
  \item{x0}{initial value of process in the time instant specified in the \code{t0} argument.}
  \item{S}{character string with the mathematical expression of the boundary.}
  \item{env}{a list object that contains the unspecified parameters in the \code{dp} object and \code{S}. 
  It is copied into a temporary environment for evaluating the mathematical expressions in the \code{dp} object and \code{S}.}
  \item{n}{number of points spanning the considered time interval.}
  \item{\dots}{additional arguments potentially passed (currently none is considered).}
}
\details{
The mathematical expression of the boundary \code{S} should be a function of \eqn{t} and may include arguments \code{t0},
\code{x0} and the parameters specified in the \code{env} argument. The \code{FPTL} function checks if the mathematical expression 
shows syntax errors and if \R can compute its symbolic derivative with respect to \eqn{t}.
		
The \code{env} argument is a list of tagged values in \code{parameter = value} form for the unspecified parameters in the \code{dp} object
and \code{S} other than \eqn{x}, \eqn{t}, \eqn{y} and \eqn{s}. Defaults to NULL, interpreted as an empty list. \R looks for the objects not found 
into the temporary environment for evaluating \code{dp} object and \code{S} in the parent.frame()
environment. \cr
	
The FPTL function for the problem of the first-passage-time of a diffusion process \eqn{X(t)}, \eqn{t_0 \leq t \leq T}{t0 \le t \le T}, through a 
boundary \eqn{S(t)} is defined as
\ifelse{latex}{\deqn{FPTL(t) = \left\{
\begin{array}{lll}
P[ X(t)>S(t) \thinspace | \thinspace X(t_0)=x_0] = 1 - F(S(t),t|x_0,t_0) & & \mbox{if \ } x_0 < S(t_0) \vspace{5pt} \\
P[ X(t)<S(t) \thinspace | \thinspace X(t_0)=x_0] = F(S(t),t|x_0,t_0) & & \mbox{if \ } x_0 > S(t_0)
\end{array}
 \right. \ , }}{\deqn{FPTL(t) = P[X(t)>S(t)|X(t0)=x0] = 1 - F(S(t),t|x0,t0) ,} 
if \eqn{x_0 < S(t_0)}{x0 < S(t0)}, and
\deqn{FPTL(t) = P[X(t)<S(t)|X(t0)=x0] = F(S(t),t|x0,t0) ,} 
if \eqn{x_0 > S(t_0)}{x0 > S(t0)}, }
where \eqn{F(x,t|y,s)}{F(x,t|y,s)} is the transition probability distribution function of the process.  \cr

Initially, the FPTL function is evaluated at a sequence of n equally spaced values from \code{t0} to \code{T}. 
Then the \code{FPTL} function makes an internal call to the \code{\link{growth.intervals}} function in order to study the growth
of the evaluation vector. Finally, the FPTL function is evaluated at a more adequate sequence of values from \code{t0} to \code{T} according 
to the abovementioned study.
}
\value{
The function \code{FPTL} computes and returns an object of class \dQuote{fptl}. It is a two-component list: 
\item{x}{a sequence of n values from \code{t0} to \code{T}.}
\item{y}{the corresponding values of the FPTL function for the \code{x} sequence.} \cr
It also includes two additional attributes: \cr \cr
\code{dp} \eqn{\qquad \qquad \qquad}{ } the object used as \code{dp} argument in the function call. \cr
\code{Call} \eqn{\qquad \qquad \quad}{ } the unevaluated function call. \cr

\code{is.fptl} returns \code{TRUE} or \code{FALSE} depending on whether its argument is an object of class \dQuote{fptl} or not. \cr

Since n is usually large, the \code{print.fptl} function does not display an object of class \dQuote{fptl} as a list, but in its \sQuote{basic}
structure instead. However, each component can be displayed separately in the usual way.
}
\references{
Román, P., Serrano, J. J., Torres, F. (2008) First-passage-time location function: Application to determine
first-passage-time densities in diffusion processes. \emph{Comput. Stat. Data Anal.}, \bold{52}, 4132--4146.
}
\author{Patricia Román-Román, Juan J. Serrano-Pérez and Francisco Torres-Ruiz.}
\seealso{
\code{\link{diffproc}} about creation of class \dQuote{diffproc} objects.

\code{\link{summary.fptl}} for summaries and \code{\link{plot.fptl}} for graphical display.
}
\examples{
% data(Lognormal)
## Continuing the diffproc(.) example:
\dontshow{Lognormal <- diffproc(c("m*x","sigma^2*x^2","dnorm((log(x)-(log(y)+(m-sigma^2/2)*(t-s)))/(sigma*sqrt(t-s)),0,1)/(sigma*sqrt(t-s)*x)", "plnorm(x,log(y)+(m-sigma^2/2)*(t-s),sigma*sqrt(t-s))"))}
## Specifying a boundary
b <- "4.5 + 4*t^2 + 7*t*sqrt(t)*sin(6*sqrt(t))"

## Computing the FPTL function and creating an object of class fptl
y <- FPTL(dp = Lognormal, t0 = 0, T = 18, x0 = 1, S = b, list(m = 0.48,
          sigma = 0.07))
y

## Testing fptl objects
is.fptl(y)
}
\keyword{classes}
\keyword{list}
\keyword{methods}
\keyword{print}
