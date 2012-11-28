report.summary.fptl <-
function (obj, tex = FALSE, digits = 8, heading = TRUE, ...) 
{
    if (!is.summary.fptl(obj)) 
        stop(paste(sQuote("obj"), "is not of class", shQuote("summary.fptl")))
    
    args <- as.list(attr(obj, "Call"))[-1]

    if (is.element("zeroSlope", names(args))) zeroSlope <- eval(args$zeroSlope) else zeroSlope <- 0.01
    if (is.element("p0.tol", names(args))) p0.tol <- eval(args$p0.tol) else p0.tol <- 8
    if (is.element("k", names(args))) k <- eval(args$k) else k <- 3

    argsFPTL <- as.list(attr(obj, "FPTLCall"))

    if (is.element("n", names(argsFPTL))) n <- eval(argsFPTL$n) else n <- 10000
	
    if (is.element("env", names(argsFPTL))){
	if (is.call(argsFPTL$env)) argsFPTL <- c(argsFPTL[3:6], unlist(as.list(argsFPTL$env)[-1]))
	else if (length(argsFPTL$env) > 0) argsFPTL <- c(argsFPTL[3:6], unlist(argsFPTL$env)) else argsFPTL <- argsFPTL[3:6]
    }
    else argsFPTL <- argsFPTL[3:6]	
	
    if (tex) { 
	  noindent <- "\\noindent "
        vskip <- "\\vskip 10pt "
        dollar <- "$"
        arrayrowsep <- " \\\\"
        ldots <- "\\ldots"
        arraycolsep <- " & = & "
        t0.label <- "t_{0}"
        labels1 <- c("t_{", "t_{", "t_{max,", "t_{max,", "t_{max,")
        labels2 <- c("}", "}^{*}", "}^{-}", "}", "}^{+}")
       	  
  	  if (heading){
		Call <- match.call()		
		Call[[1]] <- as.name("report")
		cat("\n\\vskip 20pt \\noindent \\verb$", deparse(Call), "$", sep="")
	  	if (is.name(match.call()$obj)){
			cat("\n\n", vskip, noindent, "The summary.fptl class object {\\it ", match.call()$obj, "} ", sep="")
	  	}
		else cat("\n\n", vskip, noindent, "This summary.fptl class object ", sep="")
		cat("stores the information provided by the First-Passage-Time")
	  	cat("\nLocation (FPTL) function about the location of the variation range of the first-passage-time variable through the boundary")
        	cat("\n$S(t) = ", argsFPTL$S, "$ for the diffusion process $\\{X(t) \\thinspace ; \\thinspace t_0 \\leq t \\leq T \\}$", sep="")
	  	cat("\nwith infinitesimal moments $A_1(x,t) = ", attr(obj, "dp")$mean, "$ and $A_2(x,t) = ", attr(obj, "dp")$var, "$, in the particular case", sep="")
		cat("\n", paste(c("$t_0$","$T$","$x_0$",paste("$", names(argsFPTL[-(1:4)]), "$", sep="")), lapply(argsFPTL[-4], deparse), sep= " = ", collapse = ", "), " and $P(X(", deparse(argsFPTL$t0), ") = ", deparse(argsFPTL$x0), ") = 1$.", sep="")	  	 	
	  	cat("\n\n", vskip, noindent, "It was created with the R expression", sep="")
		Call <- attr(obj, "Call")		
		Call[[1]] <- as.name("summary")	  
	  	cat("\n\\begin{verbatim}\n")
		cat(deparse(Call), fill = 80)
	  	cat("\\end{verbatim}\n")
	  }	  
    }
    else {
        noindent <- ""
        vskip <- "\n"
        dollar <- ""
        arrayrowsep <- ""
        ldots <- "..."
        arraycolsep <- " = "
        t0.label <- "t0"
        labels1 <- c("t[", "t[", "tmax[", "tmax[", "tmax[")
        labels2 <- c("]", "]*    ", "]^-", "]  ", "]^+")
	  if (heading){
	  	if (is.name(match.call()$obj)){
			cat("\nThe summary.fptl class object ", shQuote(match.call()$obj), " ", sep="")
	  	}
		else cat("\nThis summary.fptl class object ", sep="")		
		cat("stores the information provided by the First-Passage-Time Location (FPTL)")
	  	cat("\nfunction about the location of the variation range of the first-passage-time variable through the boundary")
	  	cat("\nS(t) = ", argsFPTL$S, " for the diffusion process {X(t); t0 <= t <= T} with", sep="")
		cat("\ninfinitesimal moments A1(x,t) = ", attr(obj, "dp")$mean, " and A2(x,t) = ", attr(obj, "dp")$var, ",", sep="")
		cat("\nin the particular case ", paste(names(argsFPTL[-4]), lapply(argsFPTL[-4], deparse), sep= " = ", collapse = ", "), " and P(X(", deparse(argsFPTL$t0), ") = ", deparse(argsFPTL$x0), ") = 1.", sep="")
		cat("\n\nIt was created with the R expression\n", sep="")
		Call <- attr(obj, "Call")		
		Call[[1]] <- as.name("summary")
		cat("\n\t", deparse(Call))
	  }
    }

    argsFPTL <- lapply(argsFPTL, eval)

    cat("\n", vskip, noindent, "The FPTL function was evaluated at ", 
        n, " points in ", sep = "")
    cat(dollar, "[", t0.label, ", T] = [", format(argsFPTL$t0, digits=digits), ", ", 
        format(argsFPTL$T, digits=digits), "]", dollar, ".", sep = "")
    cat(" If we consider that it is constant in")
    cat("\nthose growth subintervals where the slope of the function between the endpoints is less than ", 
        dollar, zeroSlope, dollar, " degrees,", sep = "")
    cat("\nit results that the ")

    x <- format(obj, digits=digits)

    if (nrow(obj) > 1) {
        pl <- "s"		
        cat("function starts growing from points:")
        index <- 1:nrow(x)
        labels <- matrix(paste("\t", outer(labels1, index, paste, 
            sep = ""), paste(labels2, t(x), sep = arraycolsep), 
            arrayrowsep, sep = ""), nrow = 5)
        if (tex) {
            cat("\n$$\n\\setlength{\\arraycolsep}{2pt}\n\\begin{array}{rcl}\n")
            cat(labels[1, ], sep = "\n")
            cat("\\end{array}\n$$\n")
            labels[1, ] <- "$$\n\\setlength{\\arraycolsep}{2pt}\n\\begin{array}{rcl}"
            labels <- rbind(paste("\n", vskip, noindent, "For subinterval $I_{", 
                index, "} = [", x[, 1], ", ", c(x[-1, 1], format(argsFPTL$T, digits=digits)), 
                "]$", ":", sep = ""), labels, "\\end{array}\n$$")
        }
        else {
            cat("", labels[1, ], sep = "\n")
            labels[1, ] <- paste("\n", vskip, noindent, "For subinterval I[", 
                index, "] = [", x[, 1], ", ", c(x[-1, 1], format(argsFPTL$T, digits=digits)), 
                "]", ":", sep = "")
        }
        apply(labels, 2, cat, sep = "\n")
        i <- "i"
        if (nrow(obj) > 2) 
            index <- paste(", ", dollar, "i=1,", ldots, ",", 
                nrow(obj), dollar, ", ", sep = "")
        else index <- paste(", ", dollar, "i=1,2", dollar, ", ", 
            sep = "")
    }
    else {
        pl <- ""
        cat("function starts growing from point ", dollar, 
            labels1[1], "1", labels2[1], " = ", x[1, 1], dollar, 
            ".", sep = "")
        if (tex) {
            cat("\n\n", vskip, noindent, "For subinterval $I_{1} = [", 
                x[1, 1], ", ", format(argsFPTL$T, digits=digits), "]$:\n", sep = "")
            cat("$$\n\\setlength{\\arraycolsep}{2pt}\n\\begin{array}{rcl}\n")
            cat(paste("\t", labels1[-1], "1", labels2[-1], arraycolsep, 
                x[1, -1], arrayrowsep, sep = ""), sep = "\n")
            cat("\\end{array}\n$$\n")
        }
        else {
            cat("\n\nFor subinterval I[1] = [", x[1, 1], ", ", format(argsFPTL$T, digits=digits), 
                "]:\n", sep = "")
            cat(paste("\t", labels1[-1], "1", labels2[-1], " = ", 
                x[1, -1], sep = ""), sep = "\n")
        }
        i <- "1"
        index <- " "
    }

    cat("\n", vskip, noindent, "In order to determine the time instant", 
        pl, " ", dollar, labels1[2], i, sub("\\s+$", "", labels2[2]), 
        dollar, index, "it has been considered that the value of the FPTL function is significantly", 
        sep = "")
    cat("\nbigger than the value at ", dollar, labels1[1], i, 
        labels2[1], dollar, " if its difference is over ", dollar, 
        "10^{", - p0.tol, "}", dollar, " times the increment of the function between ", 
        dollar, labels1[1], i, labels2[1], dollar, " and ", dollar, 
        labels1[4], i, sub("\\s+$", "", labels2[4]), dollar, 
        ".", sep = "")
    cat("\n\n")
}
