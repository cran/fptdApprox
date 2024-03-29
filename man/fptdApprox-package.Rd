\encoding{latin1}
\name{fptdApprox-package}
\alias{fptdApprox-package}
\alias{fptdApprox}
\docType{package}
\title{
Approximation of First-Passage-Time Densities for Diffusion Processes
}
\description{
Efficient approximation of first-passage-time (f.p.t.) densities for diffusion processes based on the First-Passage-Time Location (FPTL) function.

For a complete list of functions, use library(help=\dQuote{fptdApprox}).
}
\details{
\tabular{ll}{
Package: \tab fptdApprox\cr
Type: \tab Package\cr
Version: \tab 2.4\cr
Date: \tab 2022-08-21\cr
License: \tab GPL-2 \cr
LazyLoad: \tab yes\cr
}
The \pkg{fptdApprox} package allows to approximate efficiently the f.p.t. density for a diffusion process
through a continuous time-dependent boundary in the cases of conditioned and unconditioned f.p.t. problems.

For an unconditioned f.p.t. problem, a step by step study can be performed.
First, the diffusion process under consideration must be defined by using the function 
\code{\link{diffproc}}. Then, the \code{\link{FPTL}} function will be used to calculate the FPTL function for the specified 
process and boundary. The information provided by the function is then extracted by the method \code{\link{summary.fptl}} and used to find the range of 
variation of the f.p.t. variable. Finally, such information is used by the function \code{\link{Approx.cfpt.density}}
to obtain the approximation of the f.p.t. density.

In the general case (conditioned and unconditioned f.p.t. problems) the function \code{\link{Approx.fpt.density}}
allows to obtain directly the approximate f.p.t. density.
}
\author{Patricia Rom�n-Rom�n, Juan J. Serrano-P�rez and Francisco Torres-Ruiz.

Maintainer: Juan J. Serrano-P�rez, \email{jjserra@ugr.es}
}
\references{
Buonocore, A., Nobile, A.G. and Ricciardi, L.M. (1987) A new integral equation for the evaluation of
first-passage-time probability densities. \emph{Adv. Appl. Probab.}, \bold{19}, 784--800.

Rom�n, P., Serrano, J. J., Torres, F. (2008) First-passage-time location function: Application to determine
first-passage-time densities in diffusion processes. \emph{Comput. Stat. Data Anal.}, \bold{52}, 4132--4146.

P. Rom�n-Rom�n, J.J. Serrano-P�rez, F. Torres-Ruiz. (2012) An R package for an efficient approximation of first-passage-time densities for diffusion processes based on the FPTL function. \emph{Applied Mathematics and Computation}, \bold{218}, 8408--8428.

P. Rom�n-Rom�n, J.J. Serrano-P�rez, F. Torres-Ruiz. (2014) More general problems on first-passage times for diffusion processes: A new version of the fptdApprox R package. \emph{Applied Mathematics and Computation}, \bold{244}, 432--446.
}
