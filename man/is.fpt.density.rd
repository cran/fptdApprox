\encoding{latin1}
\name{is.fpt.density}
\alias{is.fpt.density}
\title{Testing for objects of class \dQuote{fpt.density}}
\description{
  \code{is.fpt.density} tests if its argument is an object of class \dQuote{fpt.density}.
}
\usage{
is.fpt.density(obj)
}
\arguments{
  \item{obj}{an \R object to be tested.}
}
\value{
\code{is.fpt.density} returns \code{TRUE} or \code{FALSE} depending on whether its argument is an object of class \dQuote{fpt.density} or not.

An object of class \dQuote{fpt.density} is a three-component list:
  \item{x}{a sequence of suitable time instants in \eqn{[t_0, \ T]}{[t0, T]} according to the arguments in the function call.}
  \item{y}{the approximate f.p.t. density function values on the x sequence for the unconditioned or conditioned problem at hand.}
  \item{y.x0}{NULL for a conditioned f.p.t. problem or a matrix with the values, by columns, of the approximate f.p.t. densities conditioned to each considered value
   \eqn{x_0}{x0} of the initial ditribution for an unconditioned f.p.t. problem.} \cr
  
It also includes six additional attributes. For more details, see the values of \code{\link{Approx.cfpt.density}} and \code{\link{Approx.fpt.density}} functions.
}
\author{Patricia Román-Román, Juan J. Serrano-Pérez and Francisco Torres-Ruiz.}
\examples{
## Testing fpt.density objects
\dontshow{Lognormal <- diffproc(c("m*x","sigma^2*x^2","dnorm((log(x)-(log(y)+(m-sigma^2/2)*(t-s)))/(sigma*sqrt(t-s)),0,1)/(sigma*sqrt(t-s)*x)", "plnorm(x,log(y)+(m-sigma^2/2)*(t-s),sigma*sqrt(t-s))")) ; 
b <- "4.5 + 4*t^2 + 7*t*sqrt(t)*sin(6*sqrt(t))" ; y <- FPTL(dp = Lognormal, t0 = 0, T = 18, x0 = 1, S = b, list(m = 0.48, sigma = 0.07)) ; yy <- summary(y);
yyy <- Approx.cfpt.density(yy); yyy.cp <- Approx.fpt.density(dp = Lognormal, t0 = 0, T = 18, id = 1, S = "4.5 + 4*t^2 + 7*t*sqrt(t)*sin(6*sqrt(t))", env = list(m = 0.48, sigma = 0.07))}
## Continuing the Approx.cfpt.density example:
is.fpt.density(yyy)

## Continuing the Approx.fpt.density example:
is.fpt.density(yyy.cp)
\dontrun{
is.fpt.density(yyy.ucp)}
}
