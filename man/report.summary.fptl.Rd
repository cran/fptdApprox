\encoding{latin1}
\name{report.summary.fptl}
\alias{report.summary.fptl}
\title{Writing a Report of a summary.fptl Object}
\description{
  A report is generated with the information contained in an object of class 
  \dQuote{summary.fptl}.
}
\usage{
\method{report}{summary.fptl}(obj, tex = FALSE, digits = 8, heading = TRUE, \dots)
}
\arguments{
  \item{obj}{an object of class \dQuote{summary.fptl}, a result of a call to \code{\link{summary.fptl}} function or the attribute \code{summary.fptl}
  of the value (of class \dQuote{fpt.density}) of the \code{\link{Approx.fpt.density}} function (created through one or several successive internal calls to the 
  \code{summary.fptl} function).}
  \item{tex}{logical, specifies whether to generate formatted LaTeX output (\code{TRUE}) or plain text file (\code{FALSE}).}
  \item{digits}{integer indicating the significant digits to be used.}
  \item{heading}{logical, specifies whether to include information about the diffusion process and the boundary (\code{TRUE}) 
  or not (\code{FALSE}).}
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
\code{\link{summary.fptl}} to create objects of class \dQuote{summary.fptl} from objects of class 
\dQuote{fptl}.

\code{\link{Approx.fpt.density}} to create objects of class \dQuote{fpt.density} with the attribute \code{summary.fptl} of class \dQuote{summary.fptl}

\code{\link{is.summary.fptl}} to test whether an object is of class \dQuote{summary.fptl}.
}
\examples{
## Continuing the summary.fptl(.) example:
\dontshow{Lognormal <- diffproc(c("m*x","sigma^2*x^2","dnorm((log(x)-(log(y)+(m-sigma^2/2)*(t-s)))/(sigma*sqrt(t-s)),0,1)/(sigma*sqrt(t-s)*x)", "plnorm(x,log(y)+(m-sigma^2/2)*(t-s),sigma*sqrt(t-s))")) ; 
b <- "4.5 + 4*t^2 + 7*t*sqrt(t)*sin(6*sqrt(t))" ; y <- FPTL(dp = Lognormal, t0 = 0, T = 18, x0 = 1, S = b, list(m = 0.48, sigma = 0.07)) ; yy <- summary(y);
LognormalFEx <- diffproc(c("`h(t)`*x", "sigma^2*x^2", "dnorm((log(x)-(log(y)+`H(s,t)`-(sigma^2/2)*(t - s)))/(sigma*sqrt(t-s)),0,1)/(sigma*sqrt(t-s)*x)", "plnorm(x,log(y)+ `H(s,t)`-(sigma^2/2)*(t-s),sigma*sqrt(t-s))"));
z <- FPTL(dp = LognormalFEx, t0 = 1, T = 10, x0 = 1, S = 15, env = list(sigma=0.1, `h(t)` = "t/4", `H(s,t)` = "(t^2-s^2)/8")); zz <- summary(z)}
report(yy, digits = 4)
report(yy, tex = TRUE, digits = 4)

report(zz)
report(zz, tex = TRUE)
}
\keyword{methods}
