\encoding{latin1}
\name{diffproc}
\alias{diffproc}
\alias{is.diffproc}
\alias{as.diffproc}
\alias{print.diffproc}
\title{Diffusion Processes}
\description{
\code{diffproc} creates an object of class \dQuote{diffproc} from the given set of values.
 
\code{as.diffproc} attempts to turn its argument into an object of class \dQuote{diffproc}. 

\code{is.diffproc} tests if its argument is an object of class \dQuote{diffproc}.

\code{print} shows an object of class \dQuote{diffproc}.
}
\usage{
diffproc(text)
as.diffproc(text)
is.diffproc(obj)

\method{print}{diffproc}(x, \dots)
}
\arguments{
  \item{text}{a character vector or list of length four to be coerced, containing the infinitesimal mean, infinitesimal variance, 
  transition probability density function and transition probability distribution function of a diffusion process.}
  \item{obj}{an \R object to be tested.}
  \item{x}{an object of class \dQuote{diffproc}.}
  \item{\dots}{additional arguments potentially passed (currently none is considered).}
}
\details{
The main goal of the \code{diffproc} function is to validate the mathematical expressions in the given character string vector or list.
If no errors happen, the function converts the input into an object of class \dQuote{diffproc}. Otherwise, the function reports the errors.

The mathematical expressions of the infinitesimal mean and variance should be expressions in \eqn{x} and \eqn{t}, \eqn{A_1(x,t)}{A1(x,t)} and \eqn{A_2(x,t)}{A2(x,t)}, respectively.

The mathematical expressions of the transition probability density and distribution functions should be expressions in \eqn{x}, \eqn{t}, \eqn{y} and \eqn{s}, \eqn{f(x,t|y,s)} and \eqn{F(x,t|y,s)}, respectively.

In addition, all mathematical expressions may depend on other parameters.

The function checks if the mathematical expressions show syntax errors and if \R can compute the symbolic derivative with respect 
to \eqn{x} for the infinitesimal variance and transition probability density function. 
}
\note{
The transition probability density functions specified by the \code{dnorm} function should be expressed in terms of the standard 
normal distribution, since only its derivative is considered in \R. 
}
\value{
If possible, this function returns an object of class \dQuote{diffproc} that defines a family of diffusion processes.
It is a four-component list:
  \item{mean}{character of length 1 with the mathematical expression of the infinitesimal mean of the process.}
  \item{var}{character of length 1 with the mathematical expression of the infinitesimal variance of the process.}
  \item{tpdf}{character of length 1 with the mathematical expression of the transition probability density function of the process.}
  \item{tpdF}{character of length 1 with the mathematical expression of the transition probability distribution function of the process.} 
\code{is.diffproc} returns \code{TRUE} or \code{FALSE} depending on whether its argument is an object of class \dQuote{diffproc} or not.
}
\author{Patricia Román-Román, Juan J. Serrano-Pérez and Francisco Torres-Ruiz.}
\seealso{
\code{\link{D}} to know more about the simple functions and expressions that support symbolic derivative in \R, and in 
particular about the functions provided in \R to evaluate probability density functions. 

\code{\link{FPTL}} to locate the first-passage-time of a diffusion process through a continuous boundary.
}
\examples{
## Creating a diffproc object for the lognormal diffusion process
x <- c("m*x","sigma^2*x^2",
       "dnorm((log(x)-(log(y)+(m-sigma^2/2)*(t-s)))/(sigma*sqrt(t-s)),0,1)/
       (sigma*sqrt(t-s)*x)", "plnorm(x,log(y)+(m-sigma^2/2)*(t-s),
       sigma*sqrt(t-s))") 

Lognormal <- diffproc(x)
Lognormal

## Creating a diffproc object for the Ornstein Uhlenbeck diffusion process
x <- c("alpha*x + beta","sigma^2","dnorm((x-(y*exp(alpha*(t-s))-beta*
       (1-exp(alpha*(t-s)))/alpha))/(sigma*sqrt((exp(2*alpha*(t-s)) - 1)/
       (2*alpha))),0,1)/(sigma*sqrt((exp(2*alpha*(t-s)) - 1)/(2*alpha)))",
       "pnorm(x, y*exp(alpha*(t-s)) - beta*(1 - exp(alpha*(t-s)))/alpha,
       sigma*sqrt((exp(2*alpha*(t-s)) - 1)/(2*alpha)))")
			 
OU <- diffproc(x)
OU

## Testing diffproc objects
is.diffproc(Lognormal)
is.diffproc(OU)
}
\keyword{classes}
\keyword{list}
\keyword{methods}
\keyword{print}
