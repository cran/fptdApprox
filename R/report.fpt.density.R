report.fpt.density <-
function (obj, report.sfptl = FALSE, tex = FALSE, digits = 8, 
    ...) 
{
    if (!is.fpt.density(obj)) 
        stop(paste(sQuote("obj"), "is not of class", shQuote("fpt.density")))    

    args <- as.list(attr(obj, "Call"))

    if (is.element("variableStep", names(args))) 
        variableStep <- as.logical(as.character(args$variableStep))
    else variableStep <- TRUE

    if (is.element("from.t0", names(args))) 
        from.t0 <- as.logical(as.character(args$from.t0))
    else from.t0 <- FALSE

    if (is.element("to.T", names(args))) 
        to.T <- as.logical(as.character(args$to.T))
    else to.T <- FALSE

    if (is.element("skip", names(args))) 
        skip <- as.logical(as.character(args$skip))
    else skip <- TRUE

    if (is.element("n", names(args))) 
        n <- eval(args$n)
    else n <- 250 
   
    if (is.element("tol", names(args))) 
        tol <- eval(args$tol)
    else tol <- 1e-03

    argsFPTL <- as.list(attr(attr(obj, "summary.fptl"), "FPTLCall"))
    if (is.element("env", names(argsFPTL))){
	if (is.call(argsFPTL$env)) argsFPTL <- c(argsFPTL[3:6], unlist(as.list(argsFPTL$env)[-1]))
	else if (length(argsFPTL$env) > 0) argsFPTL <- c(argsFPTL[3:6], unlist(argsFPTL$env)) else argsFPTL <- argsFPTL[3:6]
    }
    else argsFPTL <- argsFPTL[3:6]   
    
    m1 <- nrow(attr(obj, "summary.fptl"))
    m2 <- nrow(attr(obj, "Steps"))

    logic <- c(from.t0 & (eval(argsFPTL$t0) < attr(obj, "summary.fptl")[1, 
        1]), to.T & (attr(obj, "summary.fptl")[m1, 5] < eval(argsFPTL$T)))

    if (tex) {	  	  
      dollar <- "$"
      noindent <- "\\noindent "
      vskip <- "\\vskip 10pt "	  
      ldots <- "\\ldots"
      labelStep <- c("$h_{", "}")        
      t0.label <- "$t_{0}$"	  
        
	Call <- match.call()		
	Call[[1]] <- as.name("report")
	cat("\n\\vskip 20pt \\noindent \\verb$", deparse(Call), "$", sep="")
	if (is.name(match.call()$obj)){
			cat("\n\n", vskip, noindent, "The fpt.density class object {\\it ", match.call()$obj, "} ", sep="")
	}
	else cat("\n\n", vskip, noindent, "This fpt.density class object ", sep="")
	cat("stores an approximation of the first-passage-time density through the boundary")
	cat("\n$S(t) = ", argsFPTL$S, "$ for the diffusion process $\\{X(t) \\thinspace ; \\thinspace t_0 \\leq t \\leq T \\}$", sep="")		
	cat("\nwith infinitesimal moments $A_1(x,t) = ", attr(attr(obj, "summary.fptl"), "dp")$mean, "$ and $A_2(x,t) = ", attr(attr(obj, "summary.fptl"), "dp")$var, "$, in the particular case", sep="")
	cat("\n", paste(c("$t_0$","$T$","$x_0$",paste("$", names(argsFPTL[-(1:4)]), "$", sep="")), lapply(argsFPTL[-4], deparse), sep= " = ", collapse = ", "), " and $P(X(", deparse(argsFPTL$t0), ") = ", deparse(argsFPTL$x0), ") = 1$.", sep="")	
	cat("\n\n", vskip, noindent, "It was created with the R expression", sep="")	  
	cat("\n\\begin{verbatim}\n")	
	cat(deparse(attr(obj, "Call")), fill = 80)	
	cat("\\end{verbatim}")
	cat("\n\n", vskip, noindent, "The approximation process makes use of the First-Passage-Time Location (FPTL) function to locate the first-passage-time variable.\n", sep="")	
    }
    else {
      dollar <- ""
      noindent <- ""
      vskip <- "\n"
      ldots <- "..."
      labelStep <- c("h[", "]")
      t0.label <- "t0"

	if (is.name(match.call()$obj)){
		cat("\nThe fpt.density class object ", shQuote(match.call()$obj), sep="")
	}
	else cat("\nThis fpt.density class object", sep="")
	cat(" stores an approximation of the first-passage-time density through the")		
	cat("\nboundary S(t) = ", argsFPTL$S, " for the diffusion process {X(t); t0 <= t <= T} with", sep="")
	cat("\ninfinitesimal moments A1(x,t) = ", attr(attr(obj, "summary.fptl"), "dp")$mean, " and A2(x,t) = ", attr(attr(obj, "summary.fptl"), "dp")$var, ", in the particular case", sep="")
	cat("\n", paste(names(argsFPTL[-4]), argsFPTL[-4], sep= " = ", collapse = ", "), " and P(X(", deparse(argsFPTL$t0), ") = ", deparse(argsFPTL$x0), ") = 1.", sep="")	  
	cat("\n\nIt was created with the R expression", sep="")
	cat("\n\n\t", deparse(attr(obj, "Call")))
	cat("\n\nThe approximation process makes use of the First-Passage-Time Location (FPTL) function to locate the first-passage-time variable.")
    }    

    if (report.sfptl) report(attr(obj, "summary.fptl"), tex, digits, title="", head=FALSE)
    
    if (from.t0) 
        t0.label <- paste(" starting from ", t0.label, " = ", eval(argsFPTL$t0), sep = "")
    else t0.label <- ""

    if (to.T) 
        T.label <- paste(rep(" and ", from.t0), "until ", dollar, 
            "T", dollar, " = ", eval(argsFPTL$T), sep = "")
    else T.label <- ""
	
    h <- format((attr(obj, "summary.fptl")[, 3] - attr(obj, 
        "summary.fptl")[, 2])/n, digits=digits)
    intervals <- matrix(format(attr(obj, "summary.fptl")[, c(2,5)], digits=digits),ncol=2)
    intervals <- paste("[", intervals[, 1], ", ", intervals[, 2], "]", sep = "")
    
    cat("\n", vskip, noindent, "From the information provided by the FPTL function and stored in the ", sep="")
    if (is.name(args$sfptl)) cat("summary.fptl class object ", switch(1 + tex, shQuote(args$sfptl), paste("{\\it ", args$sfptl, "}", sep="")), ",", sep="")
    else cat("appropriate summary.fptl class object,")

    if (m1 > 1) {        
	  cat("\nwe must use the following integration steps in order to approximate the first-passage-time density:\n")         
        if (tex) {
		cat("\\begin{center}\n")
            cat("\\setlength{\\tabcolsep}{5pt} \n\\begin{tabular}{rcl}\n")
            cat(paste("$h_{", 1:m1, "}$ = ", h, " & in subinterval & $", 
                intervals, "$ \\\\", sep = ""), sep = "\n")
            cat("\\end{tabular}\n")
		cat("\\end{center}\n")
        }
        else cat(paste("\th[", format(1:m1), "] = ", h, "  in subinterval  ", 
            intervals, sep = ""), sep = "\n")
        cat("\n", vskip, noindent, "With the aim of determining these integration steps, we have subdivided the appropriate subintervals into ", 
            n, " parts.", sep = "")
    }
    else {        
	  cat("\nwe must use the integration step ", labelStep[1], "1", labelStep[2], dollar, " = ", h[1], 
            " in subinterval ", dollar, intervals[1], dollar, " in order to", sep = "")
        cat("\napproximate the first-passage-time density. With the aim of determining this integration step, we have divided the")
        cat("\nappropriate subinterval into ", n, " parts.", 
            sep = "")
    }
    if (m1 > 1) 
        pl <- "s"
    else pl <- ""
    index.h <- c(rep(0, logic[1]), rep(1:round(m2/2), each = 2, 
        length.out = m2 - sum(logic)), rep(m1 + 1, logic[2]))
    h.labels <- paste(labelStep[1], index.h, labelStep[2], sep = "")
    lg <- as.logical(attr(obj, "Steps")[, 4])
    if (any(lg)) {
        h.labels[lg] <- paste(h.labels[lg], "^*", sep = "")
        h.labels[!lg] <- paste(h.labels[!lg], "  ", sep = "")
    }
    h.labels <- paste(h.labels, dollar, " = ", format(attr(obj, 
        "Steps")[, 3], digits=digits), sep = "")
    intervals <- matrix(format(attr(obj, "Steps")[, 1:2], digits=digits),ncol=2)
    tstop <- paste(dollar, "t", dollar, " = ", intervals[length(attr(obj, "cumIntegral")),2], sep = "")
    intervals <- paste("[", intervals[, 1], ", ", intervals[, 2], "]", sep = "")
    cat("\n\n", vskip, noindent, "The f.p.t. density has been approximated ", 
        t0.label, T.label, " using a ", switch(1 + variableStep, 
            "fixed", "variable"), " integration step", sep = "")
    if (skip) 
        cat(" and avoiding the application \nof the numerical algorithm in those subintervals in which this is possible.")
    else cat(".")
    if (variableStep) {
        if (m2 > 1) {
            cat("\n\n", vskip, noindent, "For this specific application of the algorithm, we consider the following subintervals and integration steps:\n", 
                sep = "")
            if (tex) {
                cat("\\begin{center}\n")
	          cat("\\setlength{\\tabcolsep}{5pt} \n\\begin{tabular}{rcl}\n")
                cat(paste(h.labels, " & in subinterval & $", 
                  intervals, "$ \\\\", sep = ""), sep = "\n")
                cat("\\end{tabular}\n")
		    cat("\\end{center}\n")
            }
            else cat(paste("\t", h.labels, "  in subinterval  ", 
                intervals, sep = ""), sep = "\n")
            cat("\n", vskip, noindent, "The endlimits of the subintervals have been readjusted according to the integration steps.", 
                sep = "")
        }
        else {
            cat("\n\n", vskip, noindent, "For this specific application of the algorithm, we consider the integration step ", 
                h.labels, " in subinterval ", 
                dollar, intervals, dollar, ".", sep = "")
            cat("\nThe endlimits of the subinterval have been readjusted according to the integration step.", 
                sep = "")
        }        
    }
    else {
        cat("\n", vskip, noindent, "For this specific application of the algorithm, we consider the fixed integration step ", 
            sep = "")
        if (m1 > 1) 
            cat("given by the minimum of the integration steps ")
        cat("\nprovided by the FPTL function, that is, ", dollar, 
            "h = ", format(attr(obj, "Steps")[1, 3], digits=digits), 
            dollar, sep = "")
        if (m2 > 1) {
            cat(", over the following subintervals:", sep = "")
            if (tex) {
                cat(" \\smallskip \n\n\\centerline{\n\\begin{tabular}{l}\n")
                cat(paste("$", intervals, "$ \\\\", sep = ""), 
                  sep = "\n")
                cat("\\end{tabular}}\n")
            }
            else cat("\t", intervals, sep = "\n\t")
            cat("\n", vskip, noindent, "The endlimits of the subintervals have been readjusted according to the fixed integration step.", 
                sep = "")
        }
        else {
            cat(", in subinterval ", dollar, intervals, dollar, 
                ".", sep = "")
            cat("\nThe endlimits of the subinterval have been readjusted according to the fixed integration step.")
        }
    }
    x <- data.frame(matrix(, nrow = length(attr(obj, "cumIntegral")), ncol = 7 - 
        tex))	    
    x[, 3 - tex] <- attr(obj, "Steps")[1:length(attr(obj, "cumIntegral")), 3]    
    x[, 4 - tex] <- attr(obj, "cumIntegral")
    x[, 5 - tex] <- (attr(obj, "Steps")[1:length(attr(obj, "cumIntegral")), 
        2] - attr(obj, "Steps")[1:length(attr(obj, "cumIntegral")), 
        1])/attr(obj, "Steps")[1:length(attr(obj, "cumIntegral")), 
        3]    
    x[, 6:7 - tex] <- attr(obj, "CPUTime")
    x[attr(obj, "skips"), 5 - tex] <- 0
    it <- sum(x[, 5 - tex])
    
    if (tex) {
	  x <- format(x, digits=digits)
        x[attr(obj, "skips"), 3 - tex] <- ""
        x[, 1] <- paste("$", intervals[1:length(attr(obj, "cumIntegral"))], 
            "$", sep = "")
        cat("\n\n\\begin{table}[h]", sep = "")
        cat("\n\\centering")
        cat("\n\\caption{Approximation summary step by step} \\medskip")
        cat("\n\\label{SummaryApproxfpt}")
        cat("\n\\setlength{\\tabcolsep}{5pt} \n\\begin{tabular}{|r|r|c|r|r|r|}")
        cat("\n\\hline \\multicolumn{1}{|c|}{} &", paste("\\multicolumn{1}{c|}{", 
            c("Integration", "Cumulative", " ", "User", "System"), 
            "}", sep = "", collapse = " & "), "\\\\ \n")
        cat("\\multicolumn{1}{|c|}{Subintervals} &", paste("\\multicolumn{1}{c|}{", 
            c("step", "integral", "Iterations", "time", "time"), 
            "}", sep = "", collapse = " & "), "\\\\ ")
        cat("\n\\hline ")
        cat(apply(x, 1, paste, collapse = " & "), sep = " \\\\ \n")
        cat("\\\\ \\hline ")
        cat("\n\\end{tabular}")
        cat("\n\\end{table}")
        cat("\n\n", vskip, noindent, "Table \\ref{SummaryApproxfpt} shows the approximation process step by step.", 
            sep = "")
    }
    else {
        x[, 1:2] <- format(attr(obj, "Steps")[1:length(attr(obj, 
            "cumIntegral")), 1:2], digits=digits)
	  x <- format(x, digits=digits)
        x[attr(obj, "skips"), 3 - tex] <- ""
	  row.names(x) <- paste("Subinterval", 1:length(attr(obj, 
                "cumIntegral")))	        
        names(x) <- c("Lower end", "Upper end", "Integration step", 
            "Cumulative integral", "Iterations", "User time", 
            "System time")
        cat("\n\nThe table below shows the approximation process step by step: \n\n")
        print(x)
    }
    cumIntegral <- format(attr(obj, "cumIntegral")[length(attr(obj, 
        "cumIntegral"))], digits=digits)
    if (to.T) {
        cat("\nThe value of the cumulative integral of the approximation is ", 
            cumIntegral, ".", sep = "")
    }
    else {        
        if ((cumIntegral > (1 - tol)) & (length(attr(obj, "cumIntegral")) < 
            m2)) 
            cat("\nThe algorithm was stopped at ", tstop, ", since the value of the cumulative integral of the approximation is ", 
                cumIntegral, switch(1 + tex, " >= ", " $\\geq$ "), 
                "1 - tol.", sep = "")
        else cat("\nThe algorithm was stopped at ", tstop, " and the value of the cumulative integral of the approximation is ", 
            cumIntegral, ".", sep = "")
    }
    cat("\nThe total number of iterations is ", it,  
        " and the user time employed was ", format(sum(attr(obj, 
            "CPUTime")[, 1]), 2), " (in seconds).", sep = "")
    cat("\n\n")
}
