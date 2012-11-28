#############################################################################################
#############################################################################################
###                                                                                       ###
###                                       GOMPERTZ-TYPE PROCESS                           ###
###                                                                                       ###
### Examples for constant and another boundaries associated to some time random variables ###
### related to the growth of a Gompertz-type diffusion process with infinitesimal moments ###
### A1(x,t) = m*x*exp(-beta*t)  and  A2(x,t) = (sigma^2)*(x^2)                            ###
#############################################################################################
#############################################################################################

# Creating the diffproc object "Gompertz"

Gompertz <- diffproc(c("m*x*exp(-beta*t)","sigma^2*x^2","dnorm((log(x)-log(y) + (m/beta)*(exp(-beta*t) - exp(-beta*s))
                     + (t-s)*sigma^2/2)/(sigma*sqrt(t-s)),0,1)/(x*sigma*sqrt(t-s))","plnorm(x,log(y)
                     - (m/beta)*(exp(-beta*t) - exp(-beta*s)) - (t-s)*sigma^2/2,sigma*sqrt(t-s))"))
Gompertz

################################################################################################
# CONSTANT BOUNDARY  S                                                                         #
# Application to real data of the paper:                                                       #
#                                                                                              #
# R. Guti\u00E9rrez, P. Rom\u00E1n, D. Romero, J.J. Serrano, and F. Torres (2008)              #
# "Some time random variables related to a Gompertz-type diffusion process"                    #
# Cybernetics and Systems: An International Journal, 39, 1-13                                  #
################################################################################################

# Evaluating the FPTL function
y1G.FPTL <- FPTL(dp = Gompertz, t0 = 1, T = 30, x0 = 144, S = "2500", list(m = 0.755152, beta=0.183128, sigma = 0.0708605))
# Displaying graphically the FPTL function
win.graph()
plot(y1G.FPTL)

# Extracting the interesting information provided by the FPTL function
y1G.SFPTL <- summary(y1G.FPTL)
# Approximating the f.p.t density
y1G.g <- Approx.fpt.density(y1G.SFPTL)
# Displaying graphically the approximation of the f.p.t. density
win.graph()
plot(y1G.g)

# Reporting information of the approximation process of the f.p.t. density
report(y1G.g, report.sfptl=TRUE)

################################################################################################
# CONSTANT BOUNDARY  S = x0*exp((m/beta)*exp(-beta*t0) - 1)                                    #
# Inflection time in Application to real data of the paper:                                    #
#                                                                                              #
# R. Guti\u00E9rrez, P. Rom\u00E1n, D. Romero, J.J. Serrano, and F. Torres (2008)              #
# "Some time random variables related to a Gompertz-type diffusion process"                    #
# Cybernetics and Systems: An International Journal, 39, 1-13                                  #
################################################################################################

# Evaluating the FPTL function
y2G.FPTL <- FPTL(dp = Gompertz, t0 = 1, T = 30, x0 = 144, S = "x0*exp((m/beta)*exp(-beta*t0) - 1)", list(m = 0.755152, beta=0.183128, sigma = 0.0708605))
# Displaying graphically the FPTL function
win.graph()
plot(y2G.FPTL)

# Extracting the interesting information provided by the FPTL function
y2G.SFPTL <- summary(y2G.FPTL)
# Approximating the f.p.t density
y2G.g <- Approx.fpt.density(y2G.SFPTL)
# Displaying graphically the approximation of the f.p.t. density
win.graph()
plot(y2G.g)

# Reporting information of the approximation process of the f.p.t. density
report(y2G.g, report.sfptl=TRUE)

################################################################################################
# CONSTANT BOUNDARY  S = p*x0*exp((m/beta)*exp(-beta*t0))                                      #
# Time at which the process achieves a percentage 100*p% of the total growth in                # 
# Application to real data of the paper:                                                       #
#                                                                                              #
# R. Guti\u00E9rrez, P. Rom\u00E1n, D. Romero, J.J. Serrano, and F. Torres (2008)              #
# "Some time random variables related to a Gompertz-type diffusion process"                    #
# Cybernetics and Systems: An International Journal, 39, 1-13                                  #
################################################################################################

# Evaluating the FPTL function
y3G.FPTL <- FPTL(dp = Gompertz, t0 = 1, T = 30, x0 = 144, S = "p*x0*exp((m/beta)*exp(-beta*t0))", list(m = 0.755152, beta=0.183128, sigma = 0.0708605, p=0.5))
# Displaying graphically the FPTL function
win.graph()
plot(y3G.FPTL)

# Extracting the interesting information provided by the FPTL function
y3G.SFPTL <- summary(y3G.FPTL)
# Approximating the f.p.t density
y3G.g <- Approx.fpt.density(y3G.SFPTL)
# Displaying graphically the approximation of the f.p.t. density
win.graph()
plot(y3G.g)

# Reporting information of the approximation process of the f.p.t. density
report(y3G.g, report.sfptl=TRUE)
