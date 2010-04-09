##########################################################################################
##########################################################################################
###                                                                                    ###
###                                     WIENER PROCESS                                 ###
###                                                                                    ###
### Examples for constant and linear boundaries for which the f.p.t. density is known  ###
###                                                                                    ###
##########################################################################################
##########################################################################################

# Creating the diffproc object "Wiener"
Wiener <- diffproc(c("m","sigma^2","dnorm((x-(y+m*(t-s)))/(sigma*sqrt(t-s)),0,1)/(sigma*sqrt(t-s))",
                     "pnorm(x,y+m*(t-s),sigma*sqrt(t-s))"))
Wiener

# Theoretical expression of the f.p.t. density through a linear boundary A+C*t
gW <- parse(text="abs(A+C*t0-x0)*exp(-(A+C*t-x0-m*(t-t0))^2/(2*sigma^2*(t-t0)))/(sigma*sqrt(2*pi*(t-t0)^3))")

#################################
# CONSTANT BOUNDARY  S=4 (S>x0) #
#################################

# Evaluating the FPTL function
y1W <- FPTL(dp = Wiener, t0 = 0, T = 20, x0 = 0, S = 4, list(m = 1, sigma = 1))
# Displaying graphically the FPTL function
plot(y1W, cex=1.25)

# Extracting and showing the interesting information provided by the FPTL function
yy1W <- summary(y1W)
yy1W
# Reporting the interesting information provided by the FPTL function
report(yy1W)
# Approximating the f.p.t density 
z1W <- Approx.fpt.density(yy1W)
# Reporting information of the approximation process of the f.p.t. density
report(z1W)
# Displaying graphically the approximation of the f.p.t. density
win.graph()
plot(z1W, cex=1.25)

# Overloading the plot of the f.p.t. density
points(z1W$x[-1],eval(gW, list(t0=0,x0=0,m=1,sigma=1,A=4,C=0,t=z1W$x[-1])),type="l",col=2)

##################################################
# LINEAR BOUNDARY  S(t)=10-0.5*t (with S(t0)>x0) #
##################################################

# Evaluating the FPTL function
y2W <- FPTL(dp = Wiener, t0 = 0, T = 20, x0 = 0, S = "10-0.5*t", list(m = 1, sigma = 1))
# Displaying graphically the FPTL function
win.graph()
plot(y2W, cex=1.25)

# Extracting and showing the interesting information provided by the FPTL function
yy2W <- summary(y2W)
yy2W
# Reporting the interesting information provided by the FPTL function
report(yy2W)
# Approximating the f.p.t density 
z2W <- Approx.fpt.density(yy2W)
# Reporting information of the approximation process of the f.p.t. density
report(z2W)
# Displaying graphically the approximation of the f.p.t. density
win.graph()
plot(z2W, cex=1.25)

# Overloading the plot of the f.p.t. density
points(z2W$x[-1],eval(gW, list(t0=0,x0=0,m=1,sigma=1,A=10,C=-0.5,t=z2W$x[-1])),type="l",col=2)

#############################################
# LINEAR BOUNDARY  S=-1+t/2 (with S(t0)<x0) #
#############################################

# Evaluating the FPTL function
y3W <- FPTL(dp = Wiener, t0 = 1, T = 40, x0 = 1, S = "-1+t/2", list(m = 0, sigma = 1))
# Displaying graphically the FPTL function
win.graph()
plot(y3W, cex=1.25)

# Extracting the interesting information provided by the FPTL function
yy3W <- summary(y3W)
# Approximating the f.p.t density 
z3W <- Approx.fpt.density(yy3W)
# Displaying graphically the approximation of the f.p.t. density
win.graph()
plot(z3W, cex=1.25)

# Overloading the plot of the f.p.t. density
points(z3W$x[-1],eval(gW, list(t0=1,x0=1,m=0,sigma=1,A=-1,C=0.5,t=z3W$x[-1])),type="l",col=2)

#################################################################################################
#################################################################################################
###                                                                                           ###
###                                   LOGNORMAL PROCESS                                       ###
###                                                                                           ###
### Examples for:                                                                             ###
###               - Constant and exponential boundaries for which the f.p.t. density is known ###
###               - General boundaries                                                        ###
###                                                                                           ###
#################################################################################################
#################################################################################################

# Creating the diffproc object "Lognormal"
Lognormal <- diffproc(c("m*x","sigma^2*x^2","dnorm((log(x)-(log(y)+(m-sigma^2/2)*(t-s)))/(sigma*sqrt(t-s)),0,1)/
                     (sigma*sqrt(t-s)*x)","plnorm(x,log(y)+(m-sigma^2/2)*(t-s),sigma*sqrt(t-s))"))
Lognormal

# Theoretical expression of the f.p.t. density through an exponential boundary A*exp(B*t)
gL <- expression(abs(log(A)-log(x0)+B*t0)*exp(-(log(A)-log(x0)+B*t-(m-sigma^2/2)*(t-t0))^2/(2*sigma^2*(t-t0)))/
                 (sigma*sqrt(2*pi*(t-t0)^3)))

############################
# CONSTANT BOUNDARY S=2500 #
############################

# Evaluating and showing the FPTL function
y1L <- FPTL(dp = Lognormal, t0 = 0, T = 2, x0 = 1, S = 2500, list(m = 4, sigma = 0.001))
y1L
# Displaying graphically the FPTL function
win.graph()
plot(y1L, cex=1.25, from.t0 = FALSE, to.T = FALSE)

# Extracting and showing the interesting information provided by the FPTL function
yy1L <- summary(y1L)
yy1L
# Approximating the f.p.t density 
z1L <- Approx.fpt.density(yy1L, variableStep = TRUE, from.t0 = FALSE, to.T = FALSE, skip = TRUE)
# Reporting information of the approximation process of the f.p.t. density
report(z1L)
# Displaying graphically the approximation of the f.p.t. density
win.graph()
plot(z1L, cex=1.25)

# Overloading the plot of the f.p.t. density
points(z1L$x[-1],eval(gL, list(t0=0,x0=1,m=4,sigma=0.001,A=2500,B=0,t=z1L$x[-1])),type="l",col=2)

# Displaying graphically the approximation of the f.p.t density together with the information provided by 
# the FPTL function
win.graph()
plot(z1L,cex=1.25,growth.points=T, instants=T)

# Another options that allows to reduce the size of the output plot with good quality.
win.graph()
plot(z1L, cex = 1.8, lwd = 2, dp.legend.cex = 0.7, growth.points = TRUE, instants = TRUE)

############################################
# EXPONENTIAL BOUNDARY  S = 0.5*exp(0.4*t) #
############################################

# Evaluating the FPTL function
y2L <- FPTL(dp = Lognormal, t0 = 0, T = 25, x0 = 1, S = "0.5*exp(0.4*t)", list(m = 0.2, sigma = 0.25))
# Displaying graphically the FPTL function
win.graph()
plot(y2L, cex=1.25)

# Extracting the interesting information provided by the FPTL function
yy2L <- summary(y2L)
# Approximating the f.p.t density 
z2L <- Approx.fpt.density(yy2L)
# Displaying graphically the approximation of the f.p.t. density
win.graph()
plot(z2L, cex=1.25)

# Overloading the plot of the f.p.t. density
points(z2L$x[-1],eval(gL, list(t0=0,x0=1,m=0.2,sigma=0.25,A=0.5,B=0.4,t=z2L$x[-1])),type="l",col=2)

############################################################################
# GENERAL BOUNDARY FOR WHICH THE FPTL FUNCTION HAS TWO GROWTH SUBINTERVALS #
# Application 1 of the paper                                               #
############################################################################
  
# General boundary 
b <- "7+3.2*t+1.4*t*sin(1.75*t)"
# Evaluating the FPTL function
y3L <- FPTL(dp = Lognormal, t0 = 0, T = 10, x0 = 1, S = b, list(m = 0.48, sigma = 0.07))
# Displaying graphically the FPTL function
win.graph()
plot(y3L, cex=1.25)

# Extracting the interesting information provided by the FPTL function
yy3L<- summary(y3L)
# Reporting the interesting information provided by the FPTL function
report(yy3L)
# Reporting the interesting information provided by the FPTL function (in Latex format)
report(yy3L, tex=TRUE)
# Approximating the f.p.t density 
z3L <- Approx.fpt.density(yy3L)
# Displaying the approximation of the f.p.t. density
win.graph()
plot(z3L, cex=1.25)

# Reporting information of the approximation process of the f.p.t. density
report(z3L)
# Reporting information of the approximation process of the f.p.t. density (in Latex format)
report(z3L, tex=TRUE)

############################################################################
# GENERAL BOUNDARY. A STUDY OF THE COMPUTATIONAL COST FOR SEVERAL OPTIONS  #
# FOR THE APPROXIMATION OF THE F.P.T. DENSITY                              #
# Example and Application 2 from the paper                                 #
############################################################################

# General boundary 
b <- "4.5 + 4*t^2 + 7*t*sqrt(t)*sin(6*sqrt(t))"
# Evaluating the FPTL function
y4L <- FPTL(dp = Lognormal, t0 = 0, T = 18, x0 = 1, S = b, list(m = 0.48, sigma = 0.07))
# Displaying graphically the FPTL function
win.graph()
plot(y4L, cex=1.25)
# Extracting the interesting information provided by the FPTL function
yy4L<- summary(y4L)
# Reporting the interesting information provided by the FPTL function
report(yy4L)
# Approximating the f.p.t density 
z4L <- Approx.fpt.density(yy4L)
# Displaying graphically the approximation of the f.p.t. density
win.graph()
plot(z4L, cex=1.25)
# Reporting information of the approximation process of the f.p.t. density
report(z4L)

# Extracting the interesting information provided by the FPTL function with a value for the zeroSlope 
# lower than the default considered
yy4La <- summary(y4L, zeroSlope=0.001)
# Reporting the interesting information provided by the FPTL function (with the previous option for 
# the zeroSlope)
report(yy4La)

# Displaying graphically the FPTL function (with the information provided with the previous option for 
# the zeroSlope)
# Note that in this case the FPTL exhibits 3 growth subintervals
win.graph()
plot(y4L, yy4La, cex=1.25)
# Approximating the f.p.t density 
z4La <- Approx.fpt.density(yy4La)

# Displaying graphically the approximation of the f.p.t. density
# Note that the new approximate f.p.t. density is practically coincident with the previous.
win.graph()
plot(z4La, cex=1.25)


# COMPARATIVE STUDY: Approximating (with several options) the f.p.t density and reporting the process

# n=250, Fixed integration step, Avoiding intervals
Case1 <- Approx.fpt.density(yy4L, n=250, variableStep=FALSE, skip=TRUE)
report(Case1)

# n=250, Fixed integration step, Not Avoiding intervals
Case2 <- Approx.fpt.density(yy4L, n=250, variableStep=FALSE, skip=FALSE)
report(Case2) 

# n=250, Variable integration step, Avoiding intervals
Case3 <- Approx.fpt.density(yy4L, n=250, variableStep=TRUE, skip=TRUE)
report(Case3)

# n=250, Variable integration step, Not Avoiding intervals
Case4 <- Approx.fpt.density(yy4L, n=250, variableStep=TRUE, skip=FALSE)
report(Case4)

# n=400, Fixed integration step, Avoiding intervals
Case5 <- Approx.fpt.density(yy4L, n=400, variableStep=FALSE, skip=TRUE)
report(Case5)

# n=400, Fixed integration step, Not Avoiding intervals
Case6 <- Approx.fpt.density(yy4L, n=400, variableStep=FALSE, skip=FALSE)
report(Case6)

# n=400, Variable integration step, Avoiding intervals
Case7 <- Approx.fpt.density(yy4L, n=400, variableStep=TRUE, skip=TRUE)
report(Case7)

# n=400, Variable integration step, Not Avoiding intervals
Case8 <- Approx.fpt.density(yy4L, n=400, variableStep=TRUE, skip=FALSE)
report(Case8)


# Fixed integration step, Not Avoiding intervals, from t0 and to T
Case9 <- Approx.fpt.density(yy4L, n=250, from.t0=TRUE, to.T=TRUE, variableStep=FALSE, skip=FALSE)
report(Case9, tex=TRUE)
Case10 <- Approx.fpt.density(yy4L, n=400, from.t0=TRUE, to.T=TRUE, variableStep=FALSE, skip=FALSE)
report(Case10, tex=TRUE)


################################################################################################
# CONSTANT BOUNDARY  S=2500                                                                    #
# Case of a very concentrated f.p.t. variable                                                  #
# Example 1 of the paper:                                                                      #
#                                                                                              #
# P. Román, J.J. Serrano, F. Torres (2008)                                                     #
# "First-passage-time location function: Application to determine first-pasage-times densities #
# in diffusion processes"                                                                      #
# Computational Statistics and Data Analysis, 52, 4132-4146                                    #                                                                                              
################################################################################################

# Theoretical expression of the f.p.t. density through an exponential boundary A*exp(B*t)
gCSDA <- expression(abs(log(A)-log(x0)+B*t0)*exp(-(log(A)-log(x0)+B*t-(m-sigma^2/2)*(t-t0))^2/
+                       (2*sigma^2*(t-t0)))/(sigma*sqrt(2*pi*(t-t0)^3)))

# Evaluating the FPTL function
yCSDA <- FPTL(dp = Lognormal, t0 = 0, T = 2, x0 = 1, S = 2500, list(m = 4, sigma = 0.001))
# Displaying graphically the FPTL function
win.graph()
plot(yCSDA)

# Displaying graphically the FPTL function  (with more adequate options)
win.graph()
plot(yCSDA, cex=1.25, from.t0 = FALSE, to.T = FALSE)
 
# Extracting the interesting information provided by the FPTL function
yyCSDA <- summary(yCSDA)
# Approximating the f.p.t density with default options 
zCSDA <- Approx.fpt.density(yyCSDA, variableStep = TRUE, from.t0 = FALSE, to.T = FALSE, skip = TRUE)
# Reporting information of the approximation process of the f.p.t. density
report(zCSDA)
# Displaying the approximation of the f.p.t. density
win.graph()
plot(zCSDA, cex=1.25)
# Overloading the plot of the f.p.t. density
points(zCSDA$x[-1],eval(gCSDA, list(t0=0,x0=1,m=4,sigma=0.001,A=2500,B=0,t=zCSDA$x[-1])),type="l",col=2)


##### COMPARATIVE OF COMPUTATIONAL COSTS #####

# Total iterations = 576
z1CSDA <- Approx.fpt.density(yyCSDA, variableStep = TRUE, from.t0 = TRUE, to.T = TRUE, skip = FALSE)
report(z1CSDA)
# Total iterations = 526
z2CSDA <- Approx.fpt.density(yyCSDA, variableStep = TRUE, from.t0 = TRUE, to.T = FALSE, skip = FALSE)
report(z2CSDA)
#Total iterations = 526
z3CSDA <- Approx.fpt.density(yyCSDA, variableStep = TRUE, from.t0 = FALSE, to.T = TRUE, skip = FALSE)
report(z3CSDA)
# Total iterations = 476
z4CSDA <- Approx.fpt.density(yyCSDA, variableStep = TRUE, from.t0 = FALSE, to.T = FALSE, skip = FALSE)
report(z4CSDA)
# Total iterations = 197070 (it is recommended not run on a PC)
#z5CSDA <- Approx.fpt.density(yyCSDA, variableStep = FALSE, from.t0 = TRUE, to.T = TRUE, skip = FALSE)
#report(z5CSDA)
# Total iterations = 193021 (it is recommended not run on a PC)
#z6CSDA <- Approx.fpt.density(yyCSDA, variableStep = FALSE, from.t0 = TRUE, to.T = FALSE, skip = FALSE)
#report(z6CSDA)
# Total iterations = 4524
z7CSDA <- Approx.fpt.density(yyCSDA, variableStep = FALSE, from.t0 = FALSE, to.T = TRUE, skip = FALSE)
report(z7CSDA)
# Total iterations = 476
z8CSDA <- Approx.fpt.density(yyCSDA, variableStep = FALSE, from.t0 = FALSE, to.T = FALSE, skip = FALSE)
report(z8CSDA)


####################################################
####################################################
###                                              ###
###            GOMPERTZ-TYPE PROCESS             ###
###                                              ###
### Examples for constant and general boundaries ###
###                                              ###
####################################################
#################################################### 

# Creating the diffproc object "Gompertz"

Gompertz <- diffproc(c("m*x*exp(-beta*t)","sigma^2*x^2","dnorm((log(x)-log(y) + (m/beta)*(exp(-beta*t) - exp(-beta*s))   
                     + (t-s)*sigma^2/2)/(sigma*sqrt(t-s)),0,1)/(x*sigma*sqrt(t-s))","plnorm(x,log(y)   
                     - (m/beta)*(exp(-beta*t) - exp(-beta*s)) - (t-s)*sigma^2/2,sigma*sqrt(t-s))"))
Gompertz

################################################################################################
# CONSTANT BOUNDARY  S=exp(9)                                                                  #
# Case 1 of Example 2 of the paper:                                                            #
#                                                                                              #
# P. Román, J.J. Serrano, F. Torres (2008)                                                     #
# "First-passage-time location function: Application to determine first-pasage-times densities #
# in diffusion processes"                                                                      #
# Computational Statistics and Data Analysis, 52, 4132-4146                                    #                                                                                              
################################################################################################

# Evaluating the FPTL function
y1G <- FPTL(dp = Gompertz, t0 = 0, T = 30, x0 = 1, S = exp(9), list(m = 1, beta=0.1, sigma = 0.01))
# Displaying graphically the FPTL function
win.graph()
plot(y1G, cex=1.25)

# Extracting the interesting information provided by the FPTL function
yy1G <- summary(y1G)
# Approximating the f.p.t density 
z1G <- Approx.fpt.density(yy1G)
# Displaying graphically the approximation of the f.p.t. density
win.graph()
plot(z1G, cex=1.25)

# Reporting information of the approximation process of the f.p.t. density
report(z1G, report.sfptl=TRUE)


##########################################################################
# GENERAL BOUNDARY FOR WHICH THE FPTL FUNCTION HAS 6 GROWTH SUBINTERVALS #
##########################################################################

# Evaluating the FPTL function
y3G <- FPTL(dp = Lognormal, t0 = 0, T = 13, x0 = 1, S = "5.25+0.7*(t^2)+2*t*sin(3.75*t)", list(m = 0.48, sigma = 0.07))

# Extracting the interesting information provided by the FPTL function
yy3G <- summary(y3G)

# Approximating (with several options) the f.p.t density 

# Total iterations = 2885
z3.1G <- Approx.fpt.density(yy3G, variableStep = TRUE, from.t0 = FALSE, to.T = TRUE, skip = TRUE)
report(z3.1G)
# Total iterations = 2935
z3.2G <- Approx.fpt.density(yy3G, variableStep = TRUE, from.t0 = TRUE, to.T = TRUE, skip = TRUE)
report(z3.2G)
# Total iterations = 3185
z3.3G <- Approx.fpt.density(yy3G, variableStep = TRUE, from.t0 = TRUE, to.T = TRUE, skip = FALSE)
report(z3.3G)
# Total iterations = 11729
z3.4G <- Approx.fpt.density(yy3G, variableStep = FALSE, from.t0 = TRUE, to.T = TRUE, skip = FALSE)
report(z3.4G)

# Displaying the approximation of the f.p.t. density in the case of lesser computational cost
win.graph()
plot(z3.1G, cex=1.25)