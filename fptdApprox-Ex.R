pkgname <- "fptdApprox"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('fptdApprox')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("Approx.fpt.density")
### * Approx.fpt.density

flush(stderr()); flush(stdout())

### Name: Approx.fpt.density
### Title: Approximating First-Passage-Time Densities
### Aliases: Approx.fpt.density is.fpt.density print.fpt.density
### Keywords: classes list methods print

### ** Examples

## Continuing the summary.fptl(.) example:
## Don't show: 
Lognormal <- diffproc(c("m*x","sigma^2*x^2","dnorm((log(x)-(log(y)+(m-sigma^2/2)*(t-s)))/(sigma*sqrt(t-s)),0,1)/(sigma*sqrt(t-s)*x)", "plnorm(x,log(y)+(m-sigma^2/2)*(t-s),sigma*sqrt(t-s))")) ; 
b <- "4.5 + 4*t^2 + 7*t*sqrt(t)*sin(6*sqrt(t))" ; y <- FPTL(dp = Lognormal, t0 = 0, T = 18, x0 = 1, S = b, list(m = 0.48, sigma = 0.07)) ; yy <- summary(y)
## End Don't show
## Making an efficient approximation of the f.p.t. density 
## (optimal variable integration steps and small computational cost)
z <- Approx.fpt.density(yy)
z
print(z, digits=10)

## Making a less efficient approximation of the f.p.t. density 
## (optimal fixed integration step but high computational cost related to 
##  the efficient approximation)
## Not run: 
##D z1 <- Approx.fpt.density(yy, variableStep = FALSE, from.t0 = TRUE, to.T = 
##D                          TRUE, skip = FALSE)
##D z1
## End(Not run)

## Testing fpt.density objects
is.fpt.density(z)



cleanEx()
nameEx("FPTL")
### * FPTL

flush(stderr()); flush(stdout())

### Name: FPTL
### Title: First-Passage-Time Location Function
### Aliases: FPTL is.fptl print.fptl
### Keywords: classes list methods print

### ** Examples

## Continuing the diffproc(.) example:
## Don't show: 
Lognormal <- diffproc(c("m*x","sigma^2*x^2","dnorm((log(x)-(log(y)+(m-sigma^2/2)*(t-s)))/(sigma*sqrt(t-s)),0,1)/(sigma*sqrt(t-s)*x)", "plnorm(x,log(y)+(m-sigma^2/2)*(t-s),sigma*sqrt(t-s))"))
## End Don't show
## Specifying a boundary
b <- "4.5 + 4*t^2 + 7*t*sqrt(t)*sin(6*sqrt(t))"

## Computing the FPTL function and creating an object of class fptl
y <- FPTL(dp = Lognormal, t0 = 0, T = 18, x0 = 1, S = b, list(m = 0.48,
          sigma = 0.07))
y

## Testing fptl objects
is.fptl(y)



cleanEx()
nameEx("Integration.Steps")
### * Integration.Steps

flush(stderr()); flush(stdout())

### Name: Integration.Steps
### Title: Subintervals and Integration Steps To Approximate
###   First-Passage-Time Densities
### Aliases: Integration.Steps

### ** Examples

## Continuing the summary.fptl(.) example:
## Don't show: 
Lognormal <- diffproc(c("m*x","sigma^2*x^2","dnorm((log(x)-(log(y)+(m-sigma^2/2)*(t-s)))/(sigma*sqrt(t-s)),0,1)/(sigma*sqrt(t-s)*x)", "plnorm(x,log(y)+(m-sigma^2/2)*(t-s),sigma*sqrt(t-s))")) ; 
b <- "4.5 + 4*t^2 + 7*t*sqrt(t)*sin(6*sqrt(t))" ; y <- FPTL(dp = Lognormal, t0 = 0, T = 18, x0 = 1, S = b, list(m = 0.48, sigma = 0.07)) ; yy <- summary(y)
## End Don't show
Integration.Steps(yy)
Integration.Steps(yy, from.t0 = TRUE)
Integration.Steps(yy, to.T = TRUE, n = 100, p = 0.25)



cleanEx()
nameEx("diffproc")
### * diffproc

flush(stderr()); flush(stdout())

### Name: diffproc
### Title: Diffusion Processes
### Aliases: diffproc is.diffproc as.diffproc print.diffproc
### Keywords: classes list methods print

### ** Examples

## Creating a diffproc object for the lognormal diffusion process
x <- c("m*x","sigma^2*x^2",
       "dnorm((log(x)-(log(y)+(m-sigma^2/2)*(t-s)))/(sigma*sqrt(t-s)),0,1)/
       (sigma*sqrt(t-s)*x)", "plnorm(x,log(y)+(m-sigma^2/2)*(t-s),
       sigma*sqrt(t-s))") 

Lognormal <- diffproc(x)
Lognormal

## Creating a diffproc object for the Ornstein Uhlenbeck diffusion process
x <- c("alpha*x + beta","sigma^2","dnorm((x-(y*exp(alpha*(t-s))-beta*
       (1-exp(alpha*(t-s)))/alpha))/(sigma*sqrt((exp(2*alpha*(t-s)) - 1)/
       (2*alpha))),0,1)/(sigma*sqrt((exp(2*alpha*(t-s)) - 1)/(2*alpha)))",
       "pnorm(x, y*exp(alpha*(t-s)) - beta*(1 - exp(alpha*(t-s)))/alpha,
       sigma*sqrt((exp(2*alpha*(t-s)) - 1)/(2*alpha)))")
			 
OU <- diffproc(x)
OU

## Testing diffproc objects
is.diffproc(Lognormal)
is.diffproc(OU)



cleanEx()
nameEx("growth.intervals")
### * growth.intervals

flush(stderr()); flush(stdout())

### Name: growth.intervals
### Title: Studying the Growth of a Vector
### Aliases: growth.intervals

### ** Examples

u <- seq(0, 5, length = 200)
v <- sin(u)
w <- growth.intervals(u, v)
w

plot(u, v, type="l")
abline(v = u[as.vector(w)])

## Continuing the FPTL(.) example:
## Don't show: 
Lognormal <- diffproc(c("m*x","sigma^2*x^2","dnorm((log(x)-(log(y)+(m-sigma^2/2)*(t-s)))/(sigma*sqrt(t-s)),0,1)/(sigma*sqrt(t-s)*x)", "plnorm(x,log(y)+(m-sigma^2/2)*(t-s),sigma*sqrt(t-s))")) ; 
b <- "4.5 + 4*t^2 + 7*t*sqrt(t)*sin(6*sqrt(t))" ; y <- FPTL(dp = Lognormal, t0 = 0, T = 18, x0 = 1, S = b, list(m = 0.48, sigma = 0.07))
## End Don't show
growth.intervals(y$x, y$y)
growth.intervals(y$x, y$y, zeroSlope = 0.001)



cleanEx()
nameEx("plot.fpt.density")
### * plot.fpt.density

flush(stderr()); flush(stdout())

### Name: plot.fpt.density
### Title: Plotting Method for fpt.density Objects
### Aliases: plot.fpt.density
### Keywords: methods

### ** Examples

## Continuing the Approx.fpt.density(.) example:
## Don't show: 
Lognormal <- diffproc(c("m*x","sigma^2*x^2","dnorm((log(x)-(log(y)+(m-sigma^2/2)*(t-s)))/(sigma*sqrt(t-s)),0,1)/(sigma*sqrt(t-s)*x)", "plnorm(x,log(y)+(m-sigma^2/2)*(t-s),sigma*sqrt(t-s))")) ; 
b <- "4.5 + 4*t^2 + 7*t*sqrt(t)*sin(6*sqrt(t))" ; y <- FPTL(dp = Lognormal, t0 = 0, T = 18, x0 = 1, S = b, list(m = 0.48, sigma = 0.07)) ; yy <- summary(y) ; z <- Approx.fpt.density(yy)
## End Don't show
plot(z)
plot(z, cex = 1.25)
plot(z, from.t0 = TRUE, cex = 1.25)
plot(z, cex = 1.25, growth.points = TRUE)
plot(z, cex = 1.25, growth.points = TRUE, instants = TRUE)
plot(z, cex = 1.25, dp.legend = FALSE, growth.points = TRUE, instants = 
     TRUE)
plot(z, cex = 1.8, lwd = 2, dp.legend.cex = 0.7, growth.points = TRUE, 
     instants = TRUE)



cleanEx()
nameEx("plot.fptl")
### * plot.fptl

flush(stderr()); flush(stdout())

### Name: plot.fptl
### Title: Plotting Method for fptl Objects
### Aliases: plot.fptl
### Keywords: methods

### ** Examples

## Continuing the FPTL(.) example:
## Don't show: 
Lognormal <- diffproc(c("m*x","sigma^2*x^2","dnorm((log(x)-(log(y)+(m-sigma^2/2)*(t-s)))/(sigma*sqrt(t-s)),0,1)/(sigma*sqrt(t-s)*x)", "plnorm(x,log(y)+(m-sigma^2/2)*(t-s),sigma*sqrt(t-s))")) ; 
b <- "4.5 + 4*t^2 + 7*t*sqrt(t)*sin(6*sqrt(t))" ; y <- FPTL(dp = Lognormal, t0 = 0, T = 18, x0 = 1, S = b, list(m = 0.48, sigma = 0.07))
## End Don't show
plot(y)
plot(y, cex = 1.25)
plot(y, cex = 1.25, growth.points = FALSE)
plot(y, cex = 1.25, growth.points = FALSE, instants = FALSE)
plot(y, cex = 1.25, dp.legend = FALSE, growth.points = FALSE, instants = 
     FALSE)
plot(y, cex = 1.8, lwd = 2, dp.legend.cex = 0.7)



cleanEx()
nameEx("report.fpt.density")
### * report.fpt.density

flush(stderr()); flush(stdout())

### Name: report.fpt.density
### Title: Writing a Report of a fpt.density Object
### Aliases: report.fpt.density
### Keywords: methods

### ** Examples

## Continuing the Approx.fpt.density(.) example:
## Don't show: 
Lognormal <- diffproc(c("m*x","sigma^2*x^2","dnorm((log(x)-(log(y)+(m-sigma^2/2)*(t-s)))/(sigma*sqrt(t-s)),0,1)/(sigma*sqrt(t-s)*x)", "plnorm(x,log(y)+(m-sigma^2/2)*(t-s),sigma*sqrt(t-s))")) ; 
b <- "4.5 + 4*t^2 + 7*t*sqrt(t)*sin(6*sqrt(t))" ; y <- FPTL(dp = Lognormal, t0 = 0, T = 18, x0 = 1, S = b, list(m = 0.48, sigma = 0.07)) ; yy <- summary(y) ; z <- Approx.fpt.density(yy)
## End Don't show
report(z, digits = 4)
report(z, report.sfptl = TRUE, digits = 4)
report(z, tex = TRUE, digits = 4)
report(z, report.sfptl = TRUE, tex = TRUE, digits = 4)



cleanEx()
nameEx("report.summary.fptl")
### * report.summary.fptl

flush(stderr()); flush(stdout())

### Name: report.summary.fptl
### Title: Writing a Report of a summary.fptl Object
### Aliases: report.summary.fptl
### Keywords: methods

### ** Examples

## Continuing the summary.fptl(.) example:
## Don't show: 
Lognormal <- diffproc(c("m*x","sigma^2*x^2","dnorm((log(x)-(log(y)+(m-sigma^2/2)*(t-s)))/(sigma*sqrt(t-s)),0,1)/(sigma*sqrt(t-s)*x)", "plnorm(x,log(y)+(m-sigma^2/2)*(t-s),sigma*sqrt(t-s))")) ; 
b <- "4.5 + 4*t^2 + 7*t*sqrt(t)*sin(6*sqrt(t))" ; y <- FPTL(dp = Lognormal, t0 = 0, T = 18, x0 = 1, S = b, list(m = 0.48, sigma = 0.07)) ; yy <- summary(y)
## End Don't show
report(yy, digits = 4)
report(yy, tex = TRUE, digits = 4)



cleanEx()
nameEx("summary.fptl")
### * summary.fptl

flush(stderr()); flush(stdout())

### Name: summary.fptl
### Title: Locating the First-Passage-Time Variable
### Aliases: summary.fptl is.summary.fptl print.summary.fptl
### Keywords: classes array methods print

### ** Examples

## Continuing the FPTL(.) example:
## Don't show: 
Lognormal <- diffproc(c("m*x","sigma^2*x^2","dnorm((log(x)-(log(y)+(m-sigma^2/2)*(t-s)))/(sigma*sqrt(t-s)),0,1)/(sigma*sqrt(t-s)*x)", "plnorm(x,log(y)+(m-sigma^2/2)*(t-s),sigma*sqrt(t-s))")) ;
b <- "4.5 + 4*t^2 + 7*t*sqrt(t)*sin(6*sqrt(t))" ; y <- FPTL(dp = Lognormal, t0 = 0, T = 18, x0 = 1, S = b, list(m = 0.48, sigma = 0.07))
## End Don't show
## Summarizing an object of class fptl
yy <- summary(y)
yy
print(yy, digits=10)
yy1 <- summary(y, zeroSlope = 0.001)
yy1
yy2 <- summary(y, zeroSlope = 0.001, p0.tol = 10)
yy2

## Testing summary.fptl objects
is.summary.fptl(yy)



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
