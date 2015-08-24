\encoding{latin1}
\name{plot.fptl}
\alias{plot.fptl}
\title{Plotting Method for fptl Objects}
\description{
  This function creates a plot of the First-Passage-Time Location (FPTL) function for a first-passage-time problem, displaying the information of interest contained in an object of class \dQuote{fptl} and a corresponding object of class \dQuote{summary.fptl}.
}
\usage{
\method{plot}{fptl}(x, sfptl, from.t0 = TRUE, to.T = TRUE, dp.legend = TRUE,  
     dp.legend.cex = 1, ylab = TRUE, growth.points = TRUE, 
     instants = TRUE, \dots)
}
\arguments{
  \item{x}{an object of class \dQuote{fptl}, a result of a call to \code{\link{FPTL}}.}
  \item{sfptl}{an object of class \dQuote{summary.fptl}, a result of applying the \code{\link{summary}} method to the  \code{x}
  object.}
  \item{from.t0}{a logical value indicating whether the FPTL function should be plotted from the lower end of the
	interval considered, \eqn{t_0}{t0}, specified in the \code{x} object.}
  \item{to.T}{a logical value indicating whether the approximation should be plotted to the upper end of the
	interval considered, \eqn{T}, specified in the \code{x} object.}
  \item{dp.legend}{logical. If TRUE, adds a legend to the plot in order to identify the diffusion process 
  and boundary used in the call to \code{FPTL} function which in turn generated the \code{x} object.}
  \item{dp.legend.cex}{the magnification to be used for legend relative to the current setting of \code{cex}.}
  \item{ylab}{logical. If TRUE, adds a title for the y axis.}
  \item{growth.points}{logical. If TRUE, adds one or more vertical lines and labels to the plot in order to identify 
  the time instants from which the FPTL function starts growing.}
  \item{instants}{logical. If TRUE, draws and identify the other points of interest provided by the FPTL function.}
  \item{\dots}{graphical parameters to set before generating the plot, see \code{\link{par}}.}
}
\details{
	If the \code{sfptl} object is missing, the function makes an internal call to the \code{summary.fptl} function in order to identify the points of interest provided by the FPTL function.

	If the FPTL function shows at least a local maximum and \code{from.t0 = FALSE}, the FPTL function should be plotted from the first time instant from which the function starts growing.
	
	If the FPTL function shows at least a local maximum and \code{to.T = FALSE}, the FPTL function should be plotted to the \eqn{t_{max,m}^{+}}{tmax[m]^+} value related to the last local maximum \eqn{t_{max,m}}{tmax[m]} if the function does not decrease subsequently, or to the local minimum following the last local maximum \eqn{t_{max,m}}{tmax[m]} if the function decreases subsequently.

	If \code{dp.legend = TRUE}, a legend is placed in the top inside of the plot frame.
	
	Additional graphical arguments as \code{cex}, \code{lwd} and \code{ps} can be specified.
}
\references{
P. Román-Román, J.J. Serrano-Pérez, F. Torres-Ruiz. (2012) An R package for an efficient approximation of first-passage-time densities for diffusion processes based on the FPTL function. \emph{Applied Mathematics and Computation}, \bold{218}, 8408--8428.
}
\author{Patricia Román-Román, Juan J. Serrano-Pérez and Francisco Torres-Ruiz.}
\seealso{
\code{\link{FPTL}} to know more about the FPTL function and objects of class \dQuote{fptl}.

\code{\link{summary.fptl}} for extract the information of interest in an object of class \dQuote{fptl}.
}
\examples{
## Continuing the FPTL(.) example:
\dontshow{Lognormal <- diffproc(c("m*x","sigma^2*x^2","dnorm((log(x)-(log(y)+(m-sigma^2/2)*(t-s)))/(sigma*sqrt(t-s)),0,1)/(sigma*sqrt(t-s)*x)", "plnorm(x,log(y)+(m-sigma^2/2)*(t-s),sigma*sqrt(t-s))")) ; 
b <- "4.5 + 4*t^2 + 7*t*sqrt(t)*sin(6*sqrt(t))" ; y <- FPTL(dp = Lognormal, t0 = 0, T = 18, x0 = 1, S = b, list(m = 0.48, sigma = 0.07)); 
LognormalFEx <- diffproc(c("`h(t)`*x", "sigma^2*x^2", "dnorm((log(x)-(log(y)+`H(s,t)`-(sigma^2/2)*(t - s)))/(sigma*sqrt(t-s)),0,1)/(sigma*sqrt(t-s)*x)", "plnorm(x,log(y)+ `H(s,t)`-(sigma^2/2)*(t-s),sigma*sqrt(t-s))"));
z <- FPTL(dp = LognormalFEx, t0 = 1, T = 10, x0 = 1, S = 15, env = list(sigma=0.1, `h(t)` = "t/4", `H(s,t)` = "(t^2-s^2)/8"))}
plot(y)
plot(y, cex.main = 1.4, growth.points = FALSE)
plot(y, cex.main = 1.4, growth.points = FALSE, instants = FALSE)
plot(y, cex.main = 1.4, dp.legend = FALSE, growth.points = FALSE, instants = FALSE)
plot(y, cex = 1.25, cex.main = 1.25)
plot(y, cex = 1.25, cex.main = 1.25, dp.legend.cex = 0.8)

plot(z)
plot(z, from.t0 = FALSE)
plot(z, to.T = FALSE)
plot(z, from.t0 = FALSE, to.T = FALSE)
}
\keyword{methods}
