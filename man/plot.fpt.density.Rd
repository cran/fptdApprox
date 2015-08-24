\encoding{latin1}
\name{plot.fpt.density}
\alias{plot.fpt.density}
\title{Plotting Method for fpt.density Objects}
\description{
  This function creates a plot of the approximate density function for the conditioned or unconditioned f.p.t. problem at hand. 
  Optionally, it displays the information provided by the First-Passage-Time Location (FPTL) function for a conditioned problem.
  Moreover, for an unconditioned problem creates an additional plot with all the approximate f.p.t. 
  densities conditioned to each value \eqn{x_0}{x0} selected from the initial ditribution.
}
\usage{
\method{plot}{fpt.density}(x, from.t0, to.T, dp.legend = TRUE, dp.legend.cex = 1,  
     ylab = TRUE, growth.points = FALSE, instants = FALSE, \dots)
}
\arguments{
  \item{x}{an object of class \dQuote{fpt.density}, a result of a call to \code{\link{Approx.cfpt.density}} or \code{\link{Approx.fpt.density}} functions.}
  \item{from.t0}{a logical value indicating whether the approximation should be plotted from the lower end of the
	interval considered, \eqn{t_0}{t0}, specified in the \code{x} object.}
  \item{to.T}{a logical value indicating whether the approximation should be plotted to the upper end of the
	interval considered, \eqn{T}, specified in the \code{x} object.}
  \item{dp.legend}{logical. If TRUE, adds a legend to the plot for identifying the diffusion process and 
  boundary.}
  \item{dp.legend.cex}{the magnification to be used for legend relative to the current setting of \code{cex}.}
  \item{ylab}{logical. If TRUE, adds a title for the y axis.}
  \item{growth.points}{logical. If TRUE, for a conditioned f.p.t. problem adds one or more vertical lines and labels to the plot of the approximate 
  density function in order to identify the time instants from which the FPTL function starts growing.}
  \item{instants}{logical. If TRUE, for a conditioned f.p.t. problem adds vertical lines and labels to the plot of the approximate density function 
  in order to identify the other points of interest provided by the FPTL function.}
  \item{\dots}{graphical parameters to set before generating the plot, see \code{\link{par}}.}
}
\details{
	If \code{from.t0} or \code{to.T} arguments are missing the function considers the corresponding arguments used in the call to \code{\link{Approx.cfpt.density}} or \code{\link{Approx.fpt.density}} functions, which in turn generated the \code{x} object. 

	If the approximate density functions were calculated from the lower end of the interval considered and \code{from.t0 = FALSE}, the approximate density functions should be plotted from the lower endlimit of the first subinterval in the attribute \code{Steps} 
  of the \code{x} object.
	
	If the approximate density functions were calculated to the upper end of the interval considered and \code{to.T = FALSE}, the approximate density functions should be plotted to the upper endlimit of last subinterval in the attribute \code{Steps} of the \code{x} object.
	
	If \code{dp.legend = TRUE}, a legend is placed in the top inside of each plot frame.
	
	Additional graphical arguments as \code{cex}, \code{lwd} and \code{ps} can be specified.
}
\references{
P. Román-Román, J.J. Serrano-Pérez, F. Torres-Ruiz. (2012) An R package for an efficient approximation of first-passage-time densities for diffusion processes based on the FPTL function. \emph{Applied Mathematics and Computation}, \bold{218}, 8408--8428.

P. Román-Román, J.J. Serrano-Pérez, F. Torres-Ruiz. (2014) More general problems on first-passage times for diffusion processes: A new version of the fptdApprox R package. \emph{Applied Mathematics and Computation}, \bold{244}, 432--446.
}
\author{Patricia Román-Román, Juan J. Serrano-Pérez and Francisco Torres-Ruiz}
\seealso{
\code{\link{FPTL}} to know more about the FPTL function and objects of class \dQuote{fptl}.

\code{\link{summary.fptl}} to extract the information contained in objects of class \dQuote{fptl} and create objects of class \dQuote{summary.fptl}.

\code{\link{Approx.cfpt.density}} to approximate the f.p.t. density from objects of class \dQuote{summary.fptl} and create objects of class \dQuote{fpt.density}.

\code{\link{Approx.fpt.density}} to approximate the f.p.t. density from objects of class \dQuote{dp} and create objects of class \dQuote{fpt.density}.
}
\examples{
\dontshow{Lognormal <- diffproc(c("m*x","sigma^2*x^2","dnorm((log(x)-(log(y)+(m-sigma^2/2)*(t-s)))/(sigma*sqrt(t-s)),0,1)/(sigma*sqrt(t-s)*x)", "plnorm(x,log(y)+(m-sigma^2/2)*(t-s),sigma*sqrt(t-s))")) ; 
b <- "4.5 + 4*t^2 + 7*t*sqrt(t)*sin(6*sqrt(t))" ; y <- FPTL(dp = Lognormal, t0 = 0, T = 18, x0 = 1, S = b, list(m = 0.48, sigma = 0.07)) ; yy <- summary(y);
yyy <- Approx.cfpt.density(yy); yyy.cp <- Approx.fpt.density(dp = Lognormal, t0 = 0, T = 18, id = 1, S = "4.5 + 4*t^2 + 7*t*sqrt(t)*sin(6*sqrt(t))", env = list(m = 0.48, sigma = 0.07))}
## Continuing the Approx.cfpt.density(.) example:
plot(yyy)
plot(yyy, from.t0 = TRUE)
plot(yyy, growth.points = TRUE)
plot(yyy, growth.points = TRUE, instants = TRUE)
plot(yyy, dp.legend = FALSE, growth.points = TRUE, instants = TRUE)
plot(yyy, cex = 1.25, cex.main = 1.15)
plot(yyy, cex = 1.25, cex.main = 1.15, dp.legend.cex = 0.8, growth.points = TRUE, instants = TRUE)

## Continuing the Approx.fpt.density example:
plot(yyy.cp)
\dontrun{
plot(yyy.ucp)}

}
\keyword{methods}
