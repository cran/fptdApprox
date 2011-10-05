Approx.fpt.density <-
function (sfptl, variableStep = TRUE, from.t0 = FALSE, to.T = FALSE, 
    skip = TRUE, n = 250, p = 0.2, alpha = 1, tol = 1e-03, it.max = 50000) 
{
    if (!is.summary.fptl(sfptl)) 
        stop(paste(sQuote("sfptl"), " object is not of class ", 
            shQuote("summary.fptl"), ".", sep = ""))

    if ((!missing(p)) & (!missing(alpha)))
	   ISteps <- Integration.Steps(sfptl, variableStep, from.t0, to.T, n, p, alpha)
    else{
	   if (missing(p) & (!missing(alpha))) ISteps <- Integration.Steps(sfptl, variableStep, from.t0, to.T, n, , alpha)
	   else{
	   	  if ((!missing(p)) & missing(alpha)) ISteps <- Integration.Steps(sfptl, variableStep, from.t0, to.T, n, p)
		  else ISteps <- Integration.Steps(sfptl, variableStep, from.t0, to.T, n)
	   }
    } 

    it <- (ISteps[ ,2] - ISteps[ ,1])/ISteps[ ,3]

    args <- as.list(attr(sfptl, "FPTLCall"))
    if (is.element("env", names(args))){
	if (is.call(args$env)) args <- c(args[3:6], unlist(as.list(args$env)[-1]))
	else if (length(args$env) > 0) args <- c(args[3:6], unlist(args$env)) else args <- args[3:6]
    }
    else args <- args[3:6]
    args <- lapply(args, eval)
	
    if ((sum(it) > it.max) & interactive()){	
    	itbyStep <- format(data.frame(matrix(ISteps[ ,1:2], ncol=2), it), digits=8)
	logic <- as.logical(ISteps[ ,4])	
	if ((skip) & (any(logic))) itbyStep[ ,3][logic] <- paste("0 or", itbyStep[ ,3][logic])
    	row.names(itbyStep) <- paste("Subinterval", 1:nrow(itbyStep))	        
    	names(itbyStep) <- c("Lower end", "Upper end", "Iterations")
    	cat("\nThe table below shows the number of iterations of the approximation process step by step: \n\n")
    	print(itbyStep)
	iterations <- sum(it)

	if ((skip) & (any(logic))) cat("\nIf no interval is avoided, the total number of iterations is", iterations) 
	else cat("\nThe total number of iterations is", iterations)

	if (to.T) cat(".\n")
	else{
		if ((from.t0) & (args$t0 < sfptl[1,2])) itStop <- (cumsum(it)[!ISteps[ ,4]])[-1]		
		else itStop <- cumsum(it)[!ISteps[ ,4]]		
		itStop <- itStop[itStop < iterations]				
		if (length(itStop) >=1){
			cat(", although the algorithm can stop at ")
			if (length(itStop) == 1) cat("iteration ", itStop, ".\n", sep="")
			else cat("iterations ", paste(itStop[-length(itStop)], collapse=", "), " and ", itStop[length(itStop)], ".\n", sep="")
		}
		else cat(".\n")
	}

	repeat{ 
		answer <- readline("Do you want to continue? (y/n) ")
		if ((answer == "n") | (answer == "y")) break
	}   	
	if (answer == "n") stop("Approximation process stopped by the usuary.\n\n")    		  	
    }
        
    g <- alist(x = , y = )	
    attr(g, "Call") <- match.call()		
     
    aux <- function(x, envir) eval(parse(text=paste("substitute(", deparse(x), ", envir)", sep="")))
    if (length(attr(g, "Call")) > 2){	
    	vars <- lapply(attr(g, "Call")[-(1:2)], all.vars)	   
    	logic <- (unlist(lapply(vars, length)) > 0)
    	if (any(logic)){								
		envars <- lapply(vars[logic], mget, envir=parent.frame())				
		attr(g, "Call")[-(1:2)][logic] <- mapply(aux, attr(g, "Call")[-(1:2)][logic], envars, SIMPLIFY = FALSE)		
    	}	
    	logic <- unlist(lapply(lapply(lapply(attr(g, "Call")[-(1:2)], all.names), is.element, el=c("[","[[")), any))
    	if (any(logic)) attr(g, "Call")[-(1:2)][logic] <- lapply(attr(g, "Call")[-(1:2)][logic], eval)
    }
    
    attr(g, "Steps") <- ISteps

    args$S <- eval(parse(text = paste("substitute(", args$S, ", args)", sep="")))
	
    rho <- sign(eval(args$S, list(t = args$t0)) - args$x0)	
    if (rho == 0) 
        stop("S(t0) = x0")	

    if (inherits(try(args$S, silent = TRUE), "try-error")) 
        stop(paste("the mathematical expression of the boundary shows syntax errors.", 
            sep = ""))
    if (inherits(try(D(args$S, "t"), silent = TRUE), 
        "try-error")) 
        stop("R can not compute the symbolic derivative with respect to 't' of the mathematical expression of the boundary")

    cat("\nComputing...\t")

    DB <- deriv(args$S, "t", function.arg = "t")
    A1 <- function(x, t) NULL
    body(A1) <- eval(parse(text = paste("substitute(", attr(sfptl, 
        "dp")$mean, ", args)", sep = "")))
    DA2 <- deriv(eval(parse(text = paste("substitute(", attr(sfptl, 
        "dp")$var, ", args)", sep = ""))), "x", function.arg = c("x", 
        "t"))
    Df <- deriv(eval(parse(text = paste("substitute(", attr(sfptl, 
        "dp")$tpdf, ", args)", sep = ""))), "x", function.arg = c("x", 
        "t", "y", "s"))

    now <- matrix(, nrow = nrow(ISteps) + 1, ncol = 2)
    now[1, ] <- proc.time()[1:2]
    jumps <- logical(nrow(ISteps))
    Integral <- numeric(nrow(ISteps))
    h <- numeric(0)
    F <- 0

    g0 <- ISteps[1, 1]    
    g$x <- ISteps[1, 1] + ISteps[1, 3]
    b <- DB(g$x)
    a2 <- DA2(b, g$x)
    ff <- Df(b, g$x, args$x0, args$t0)
    g$y <- max(0, -rho * (ff * (attr(b, "gradient") - A1(b, g$x) + 
        0.75 * attr(a2, "gradient")) + attr(ff, "gradient") * 
        a2))    
    Integral[1] <- ISteps[1, 3] * g$y/2
    ISteps[1, 1] <- g$x     
    
    for (i in 1:nrow(ISteps)) {
        if (skip & ISteps[i, 4] & (g$y[length(g$y)] < 10^-10)) {
            correction <- h[length(h)] * g$y[length(g$y)]/2
            Integral[i - 1] <- Integral[i - 1] - correction
            F <- F - correction
            g$y[length(g$y)] <- 0
            h <- c(h, ISteps[i, 2] - ISteps[i, 1])
            g$x <- c(g$x, ISteps[i, 2])
            b <- c(b, DB(ISteps[i, 2]))
            g$y <- c(g$y, 0)
            now[i + 1, ] <- proc.time()[1:2]
            Integral[i] <- 0
            jumps[i] <- TRUE
        }
        else {
            u2 <- seq(ISteps[i, 1] + ISteps[i, 3], ISteps[i, 2], by = ISteps[i, 3])
            h2 <- rep(ISteps[i, 3], length(u2))
            b2 <- DB(u2)
            a1 <- A1(b2, u2)
            a2 <- DA2(b2, u2)
            a <- as.vector(attr(b2, "gradient")) - a1 + 0.75 * 
                as.vector(attr(a2, "gradient"))
            f0 <- Df(b2, u2, args$x0, args$t0)

            if (length(a) == 1) 
                a <- rep(a, length(u2))

            if (length(a2) == 1) 
                a2 <- rep(a2, length(u2))

            g$y <- c(g$y, -rho * (f0 * a + as.vector(attr(f0, 
                "gradient")) * a2))

            n <- length(g$x)
            h <- c(h, h2)
            g$x <- c(g$x, u2)

            if (length(b2) == 1) 
                b2 <- rep(b2, length(u2))
            b <- c(b, b2)

            for (k in (1:length(u2))) {
                index <- 1:(n + k - 1)
                f1 <- Df(b2[k], u2[k], b[index], g$x[index])
                if (length(f1) == 1) {
                  Df1 <- attr(f1, "gradient")
                  f1 <- rep(f1, n + k - 1)
                  attr(f1, "gradient") <- Df1
                }
                g$y[n + k] <- max(0, g$y[n + k] + rho * sum(h[index] * 
                  g$y[index] * (f1 * a[k] + as.vector(attr(f1, 
                  "gradient")) * a2[k])))
            }

            Integral[i] <- Integral[i] + sum(h2 * (g$y[n + (1:length(u2))] + 
                g$y[n - 1 + (1:length(u2))]))/2
            F <- F + Integral[i]
            now[i + 1, ] <- proc.time()[1:2]
            if ((!to.T) & (!ISteps[i, 4]) & (F >= (1 - tol))) 
                break
        }
    }
    
    cat("Done.\n")

    if (F < 1 - tol){ 
		cat(paste("\nThe value of the cumulative integral of the approximation is", F, "< 1 - tol.\n"))
		if (to.T) cat("It may be appropriate to extend the considered time interval.\n") 
    else{
        cat("If the value of the cumulative integral is not high and the final stopping instant is less than T, it may be appropriate:")
        cat("\n   - Check if the value of the final stopping instant increases using k argument to summary the fptl class object, or")  
        cat("\n   - Approximate the density again with to.T = TRUE.\n")
    }
    }
    
    g$x <- c(g0, g$x)
    g$y <- c(0, g$y)
    if (i < nrow(ISteps)) {
        Integral <- Integral[1:i]
        jumps <- jumps[1:i]
        now <- now[1:(i + 1), ]
    }
    CPUTime <- apply(now, 2, diff)
    if (!is.matrix(CPUTime)) 
        CPUTime <- matrix(CPUTime, ncol = 2)
    attr(g, "cumIntegral") <- cumsum(Integral)
    attr(g, "skips") <- jumps
    attr(g, "CPUTime") <- CPUTime
    attr(g, "summary.fptl") <- sfptl
    class(g) <- c("fpt.density", "list")    
    return(g)
}

