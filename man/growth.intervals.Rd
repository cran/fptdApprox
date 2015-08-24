\encoding{latin1}
\name{growth.intervals}
\alias{growth.intervals}
\title{Studying the Growth of a Vector}
\description{
	For the vector of values resulting from the evaluation of a function, this function determines the positions of the values between which the vector grows. 
}
\usage{
growth.intervals(x, y, zeroSlope = 0.01)
}
\arguments{
  \item{x}{a vector of values.}
  \item{y}{the vector of the corresponding values of a function for the \code{x} vector.}
  \item{zeroSlope}{maximum slope (in degrees) required to consider that a growing function is constant.}
}
\details{
The \code{growth.intervals} function ignores the pairs of values between which the vector \code{y} starts and stops growing if the slope is less than \code{zeroSlope}. 
}
\value{
	This function returns NULL if the vector \code{y} is not growing. Otherwise, the
	function computes a matrix of indexes with two columns. The first column 
	contains the positions of the values from which the vector \code{y} starts growing and the	
	second column those at which the vector \code{y} stops growing.
}
\author{Patricia Román-Román, Juan J. Serrano-Pérez and Francisco Torres-Ruiz.}
\examples{
u <- seq(0, 5, length = 200)
v <- sin(u)
w <- growth.intervals(u, v)
w

plot(u, v, type = "l", las = 1)
abline(v = u[as.vector(w)])

## Continuing the FPTL(.) examples:
\dontshow{Lognormal <- diffproc(c("m*x","sigma^2*x^2","dnorm((log(x)-(log(y)+(m-sigma^2/2)*(t-s)))/(sigma*sqrt(t-s)),0,1)/(sigma*sqrt(t-s)*x)", "plnorm(x,log(y)+(m-sigma^2/2)*(t-s),sigma*sqrt(t-s))")) ; 
b <- "4.5 + 4*t^2 + 7*t*sqrt(t)*sin(6*sqrt(t))" ; y <- FPTL(dp = Lognormal, t0 = 0, T = 18, x0 = 1, S = b, list(m = 0.48, sigma = 0.07))}
growth.intervals(y$x, y$y)
growth.intervals(y$x, y$y, zeroSlope = 0.001)
}
