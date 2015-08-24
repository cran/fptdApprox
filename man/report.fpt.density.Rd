\encoding{latin1}
\name{report.fpt.density}
\alias{report.fpt.density}
\title{Writing a Report of a fpt.density Object}
\description{
  A report is generated with the information contained in an object of class \dQuote{fpt.density}.
}
\usage{
\method{report}{fpt.density}(obj, report.sfptl = FALSE, tex = FALSE, digits = 8, \dots)
}
\arguments{
  \item{obj}{an object of class \dQuote{fpt.density}, a result of a call to \code{\link{Approx.cfpt.density}} or \code{\link{Approx.fpt.density}} functions.}
  \item{report.sfptl}{logical. If TRUE, also writes a report with the information contained in the object of class \dQuote{summary.fptl} 
  from which the \code{obj} object was created.} 
  \item{tex}{logical, specifies whether to generate formatted LaTeX output (\code{TRUE}) or plain text file (\code{FALSE}).}
  \item{digits}{integer indicating the significant digits to be used.}
  \item{\dots}{additional arguments potentially passed (currently none is considered).}
}
\details{
  The resulting report can be recicled by copy-and-pasting into a document (if \code{tex = FALSE}) 
  or directly included into a LaTeX file (if \code{tex = TRUE}).
}
\value{
  Send the report to the command windows.
}
\references{
P. Román-Román, J.J. Serrano-Pérez, F. Torres-Ruiz. (2012) An R package for an efficient approximation of first-passage-time densities for diffusion processes based on the FPTL function. \emph{Applied Mathematics and Computation}, \bold{218}, 8408--8428.

P. Román-Román, J.J. Serrano-Pérez, F. Torres-Ruiz. (2014) More general problems on first-passage times for diffusion processes: A new version of the fptdApprox R package. \emph{Applied Mathematics and Computation}, \bold{244}, 432--446.
}
\author{Patricia Román-Román, Juan J. Serrano-Pérez and Francisco Torres-Ruiz.}
\seealso{
\code{\link{Approx.cfpt.density}} to approximate the f.p.t. density from objects of class \dQuote{summary.fptl} and to create objects of class \dQuote{fpt.density}.

\code{\link{Approx.fpt.density}} to approximate the f.p.t. density from objects of class \dQuote{dp} and to create objects of class \dQuote{fpt.density}.

\code{\link{is.fpt.density}} to test whether an object is an object of class \dQuote{fpt.density}.
}
\examples{
\dontshow{Lognormal <- diffproc(c("m*x","sigma^2*x^2","dnorm((log(x)-(log(y)+(m-sigma^2/2)*(t-s)))/(sigma*sqrt(t-s)),0,1)/(sigma*sqrt(t-s)*x)", "plnorm(x,log(y)+(m-sigma^2/2)*(t-s),sigma*sqrt(t-s))")) ; 
b <- "4.5 + 4*t^2 + 7*t*sqrt(t)*sin(6*sqrt(t))" ; y <- FPTL(dp = Lognormal, t0 = 0, T = 18, x0 = 1, S = b, list(m = 0.48, sigma = 0.07)) ; yy <- summary(y);
yyy <- Approx.cfpt.density(yy); yyy.cp <- Approx.fpt.density(dp = Lognormal, t0 = 0, T = 18, id = 1, S = "4.5 + 4*t^2 + 7*t*sqrt(t)*sin(6*sqrt(t))", env = list(m = 0.48, sigma = 0.07))}
## Continuing the Approx.cfpt.density(.) example:
report(yyy, digits = 4)
report(yyy, report.sfptl = TRUE, digits = 4)
report(yyy, tex = TRUE, digits = 4)
report(yyy, report.sfptl = TRUE, tex = TRUE, digits = 4)

## Continuing the Approx.fpt.density(.) example:
report(yyy.cp)
\dontrun{
report(yyy.ucp)}
}
\keyword{methods}
