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
FPTL(dp, t0, T, x0, S, env = NULL, n = 4000)

is.fptl(obj)

\method{print}{fptl}(x, \dots)
}
\arguments{
  \item{dp}{an object of class \dQuote{diffproc} defining a family of diffusion processes.}
  \item{obj}{an \R object to be tested.}
  \item{x}{an object of class \dQuote{fptl}, a result of a call to this function.}
  \item{t0, T}{lower and upper limits of the considered time interval. Must be finite.}
  \item{x0}{fixed initial value of process in the time instant specified in the \code{t0} argument.}
  \item{S}{numerical value of a constant boundary or character string with the mathematical expression of a time-dependent boundary.}
  \item{env}{a named list of objects of numeric or character type specifying the values of names which occur in
  the mathematical expressions in objects \code{dp} and \code{S}, or of names which occur in the specification of the previous values.
  Defaults to NULL, interpreted as an empty list.} 
  \item{n}{number of points at which the FPTL function is evaluated.}
  \item{\dots}{additional arguments potentially passed (currently none is considered).}
}
\details{
The FPTL function for the problem of the first-passage-time of a diffusion process \eqn{\{X(t), t_0 \leq t \leq T \}}{{X(t), t0 \le t \le T}},  
conditioned to \eqn{X(t_0) = x_0}{X(t0) = x0}, through a continuous boundary \eqn{S(t)} is defined as
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
to the abovementioned study. \cr

The mathematical expression of the boundary \code{S} should be a function of \eqn{t} and may include arguments \code{t0},
\code{x0} and the parameters specified in the \code{env} argument. The \code{FPTL} function checks if the mathematical expression 
shows syntax errors and if \R can compute its symbolic derivative with respect to \eqn{t}. \cr
	
The \code{env} argument is a list of tagged values in \code{name = value} form with \code{name} other than \eqn{x}, \eqn{t}, \eqn{y} and \eqn{s}.
To name the expression of a function \eqn{h(u)} as a character string we can use 
\eqn{\grave{} \thinspace}{}\if{text}{`}\if{html}{\out{&#96;}}\code{h(u)}\eqn{\grave{} \ }{}\if{text}{`}\if{html}{\out{&#96;}} \code{= value} 
if we want to show its dependence on \eqn{u}, or \code{h = value} otherwise. \cr

The \code{env} argument is copied into a temporary environment for evaluating the mathematical expressions in objects \code{dp} and \code{S}. 
\R looks for the objects not found into this temporary environment in the parent.frame() environment. 
}
\value{
The function \code{FPTL} computes and returns an object of class \dQuote{fptl}. It is a two-component list: 
\item{x}{a sequence of n values from \code{t0} to \code{T}.}
\item{y}{the corresponding values of the FPTL function for the \code{x} sequence.} \cr
It also includes three additional attributes:
\tabular{rl}{
\code{Call} \tab the unevaluated function call, substituting each name in this call for its value when \cr 
\code{} \tab the latter has length 1. \cr
\code{dp} \tab the object used as \code{dp} argument in the function call. \cr 
\code{vars} \tab NULL or a list containing the values of names in the function call for those names \cr 
\code{} \tab with values of length greater than 1. \cr
}

\code{is.fptl} returns \code{TRUE} or \code{FALSE} depending on whether its argument is an object of class \dQuote{fptl} or not. \cr

Since n is usually large, the \code{print.fptl} function does not display an object of class \dQuote{fptl} as a list, but in its \sQuote{basic}
structure instead. However, each component can be displayed separately in the usual way.
}
\references{
Román, P., Serrano, J. J., Torres, F. (2008) First-passage-time location function: Application to determine
first-passage-time densities in diffusion processes. \emph{Comput. Stat. Data Anal.}, \bold{52}, 4132--4146.

P. Román-Román, J.J. Serrano-Pérez, F. Torres-Ruiz. (2012) An R package for an efficient approximation of first-passage-time 
densities for diffusion processes based on the FPTL function. \emph{Applied Mathematics and Computation}, \bold{218}, 8408--8428.
}
\author{Patricia Román-Román, Juan J. Serrano-Pérez and Francisco Torres-Ruiz.}
\seealso{
\code{\link{diffproc}} about creation of class \dQuote{diffproc} objects.

\code{\link{summary.fptl}} for summaries and \code{\link{plot.fptl}} for graphical display.
}
\examples{
% data(Lognormal)
## Continuing the diffproc(.) examples:
\dontshow{Lognormal <- diffproc(c("m*x","sigma^2*x^2","dnorm((log(x)-(log(y)+(m-sigma^2/2)*(t-s)))/(sigma*sqrt(t-s)),0,1)/(sigma*sqrt(t-s)*x)", "plnorm(x,log(y)+(m-sigma^2/2)*(t-s),sigma*sqrt(t-s))"));
LognormalFEx <- diffproc(c("`h(t)`*x", "sigma^2*x^2", "dnorm((log(x)-(log(y)+`H(s,t)`-(sigma^2/2)*(t - s)))/(sigma*sqrt(t-s)),0,1)/(sigma*sqrt(t-s)*x)", "plnorm(x,log(y)+ `H(s,t)`-(sigma^2/2)*(t-s),sigma*sqrt(t-s))"))}
## Specifying a boundary
b <- "4.5 + 4*t^2 + 7*t*sqrt(t)*sin(6*sqrt(t))"

## Computing FPTL functions and creating objects of class fptl
y <- FPTL(dp = Lognormal, t0 = 0, T = 18, x0 = 1, S = b, env = list(m = 0.48,
          sigma = 0.07))
y

z <- FPTL(dp = LognormalFEx, t0 = 1, T = 10, x0 = 1, S = 15, env = list(sigma=0.1, 
          `h(t)` = "t/4", `H(s,t)` = "(t^2-s^2)/8"))
z
          
## Testing fptl objects
is.fptl(y)
is.fptl(z)
}
\keyword{classes}
\keyword{list}
\keyword{methods}
\keyword{print}
