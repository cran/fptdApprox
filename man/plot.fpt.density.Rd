\encoding{latin1}
\name{plot.fpt.density}
\alias{plot.fpt.density}
\title{Plotting Method for fpt.density Objects}
\description{
  This function creates a plot of the approximate density function for a first-passage-time (f.p.t.) problem. 
  Optionally, it displays the information provided by the First-Passage-Time Location (FPTL) function. 
}
\usage{
\method{plot}{fpt.density}(x, from.t0, to.T, dp.legend = TRUE, dp.legend.cex = 1,  
     ylab = TRUE, growth.points = FALSE, instants = FALSE, \dots)
}
\arguments{
  \item{x}{an object of class \dQuote{fpt.density}, a result of a call to \code{\link{Approx.fpt.density}}.}
  \item{from.t0}{a logical value indicating whether the approximation should be plotted from the lower end of the
	interval considered, \eqn{t_0}{t0}, specified in the \code{x} object.}
  \item{to.T}{a logical value indicating whether the approximation should be plotted to the upper end of the
	interval considered, \eqn{T}, specified in the \code{x} object.}
  \item{dp.legend}{logical. If TRUE, adds a legend to the plot for identifying the diffusion process and 
  boundary used in the call to \code{\link{FPTL}} function which in turn generated the object of class \dQuote{fptl} 
  from whose summary the \code{x} object was finally created.}
  \item{dp.legend.cex}{the magnification to be used for legend relative to the current setting of \code{cex}.}
  \item{ylab}{logical. If TRUE, adds a title for the y axis.}
  \item{growth.points}{logical. If TRUE, adds one or more vertical lines and labels to the plot in order to identify 
  the time instants from which the FPTL function starts growing.}
  \item{instants}{logical. If TRUE, adds vertical lines and labels to the plot in order to identify the other points 
  of interest provided by the FPTL function.}
  \item{\dots}{graphical parameters to set before generating the plot, see \code{\link{par}}.}
}
\details{
	If \code{from.t0} or \code{to.T} arguments are missing the function considers the corresponding arguments used in 
	the call to \code{\link{Approx.fpt.density}} function, which in turn generated the \code{x} object. 

	If the approximate density function was calculated from the lower end of the interval considered and \code{from.t0 = FALSE}, 
	the approximate density function should be plotted from the first point from which the FPTL function starts growing.
	
	If the approximate density function was calculated to the upper end of the interval considered and \code{to.T = FALSE}, 
	the approximate density function should be plotted to the last point of interest provided by the FPTL function.
	
	If \code{dp.legend = TRUE}, a legend is placed in the top inside of the plot frame.
	
	Additional graphical arguments as \code{cex}, \code{lwd} and \code{ps} can be specified.
}
\author{Patricia Román-Román, Juan J. Serrano-Pérez and Francisco Torres-Ruiz}
\seealso{
\code{\link{FPTL}} to know more about the FPTL function and objects of class \dQuote{fptl}.

\code{\link{summary.fptl}} to extract the information contained in objects of class \dQuote{fptl} and create objects of class \dQuote{summary.fptl}.

\code{\link{Approx.fpt.density}} to approximate the f.p.t. density from objects of class \dQuote{summary.fptl} and create objects of class \dQuote{fpt.density}.
}
\examples{
## Continuing the Approx.fpt.density(.) example:
\dontshow{Lognormal <- diffproc(c("m*x","sigma^2*x^2","dnorm((log(x)-(log(y)+(m-sigma^2/2)*(t-s)))/(sigma*sqrt(t-s)),0,1)/(sigma*sqrt(t-s)*x)", "plnorm(x,log(y)+(m-sigma^2/2)*(t-s),sigma*sqrt(t-s))")) ; 
b <- "4.5 + 4*t^2 + 7*t*sqrt(t)*sin(6*sqrt(t))" ; y <- FPTL(dp = Lognormal, t0 = 0, T = 18, x0 = 1, S = b, list(m = 0.48, sigma = 0.07)) ; yy <- summary(y) ; z <- Approx.fpt.density(yy)}
plot(z)
plot(z, from.t0 = TRUE)
plot(z, growth.points = TRUE)
plot(z, growth.points = TRUE, instants = TRUE)
plot(z, dp.legend = FALSE, growth.points = TRUE, instants = TRUE)
plot(z, cex = 1.25, cex.main = 1.15)
plot(z, cex = 1.25, cex.main = 1.15, dp.legend.cex = 0.8, growth.points = TRUE, instants = TRUE)
}
\keyword{methods}
