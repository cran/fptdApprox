\encoding{latin1}
\name{report}
\alias{report}
\title{Writing a Report of an Object}
\description{
  \code{report} is a function used to generate a report of the results of various statistical methods.  
  The function invokes particular \sQuote{methods} which depend on the \sQuote{class} of the first argument.  
}
\usage{
report(obj, \dots)
}
\arguments{
  \item{obj}{an object for which a report is desired.}
  \item{\dots}{additional arguments affecting the report generated.}  
}
\details{  
  The functions \code{\link{report.summary.fptl}} and \code{\link{report.fpt.density}} are examples
  of particular methods which report objects of class \dQuote{summary.fptl} and \dQuote{fpt.density}.
}
\value{
  The form of the value returned by \code{report} depends on the class of its argument.  
}
\author{Patricia Román-Román, Juan J. Serrano-Pérez and Francisco Torres-Ruiz.}
\seealso{
  \code{\link{report.summary.fptl}}
  
  \code{\link{report.fpt.density}} 
}
\keyword{methods}
