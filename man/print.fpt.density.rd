\encoding{latin1}
\name{print.fpt.density}
\alias{print.fpt.density}
\title{Printing First-Passage-Time Densities}
\description{
  \code{print.fpt.density} shows an object of class \dQuote{fpt.density}.
}
\usage{
\method{print}{fpt.density}(x, \dots)
}
\arguments{
  \item{x}{an object of class \dQuote{fpt.density}.}
  \item{\dots}{further arguments passed to \code{\link{print}} and \code{\link{format}} methods.}
}
\value{
Since the length of components of an object of class \dQuote{fpt.density} is usually large, the \code{print.fpt.density} function does not display such object as a list, but in its \sQuote{basic} structure instead. However, each component can be displayed separately in the usual way.
}
\references{
P. Román-Román, J.J. Serrano-Pérez, F. Torres-Ruiz. (2012) An R package for an efficient approximation of first-passage-time densities for diffusion processes based on the FPTL function. \emph{Applied Mathematics and Computation}, \bold{218}, 8408--8428.

P. Román-Román, J.J. Serrano-Pérez, F. Torres-Ruiz. (2014) More general problems on first-passage times for diffusion processes: A new version of the fptdApprox R package. \emph{Applied Mathematics and Computation}, \bold{244}, 432--446.
}
\author{Patricia Román-Román, Juan J. Serrano-Pérez and Francisco Torres-Ruiz.}
\examples{
\dontshow{Lognormal <- diffproc(c("m*x","sigma^2*x^2","dnorm((log(x)-(log(y)+(m-sigma^2/2)*(t-s)))/(sigma*sqrt(t-s)),0,1)/(sigma*sqrt(t-s)*x)", "plnorm(x,log(y)+(m-sigma^2/2)*(t-s),sigma*sqrt(t-s))")) ; 
b <- "4.5 + 4*t^2 + 7*t*sqrt(t)*sin(6*sqrt(t))" ; y <- FPTL(dp = Lognormal, t0 = 0, T = 18, x0 = 1, S = b, list(m = 0.48, sigma = 0.07)) ; yy <- summary(y);
yyy <- Approx.cfpt.density(yy); yyy.cp <- Approx.fpt.density(dp = Lognormal, t0 = 0, T = 18, id = 1, S = "4.5 + 4*t^2 + 7*t*sqrt(t)*sin(6*sqrt(t))", env = list(m = 0.48, sigma = 0.07))}
## Continuing the Approx.cfpt.density example:
yyy
print(yyy, digits=10)

## Continuing the Approx.fpt.density example:
yyy.cp
\dontrun{
yyy.ucp}
}
