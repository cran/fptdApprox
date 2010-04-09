report.fpt.density <-
function (obj, report.sfptl = FALSE, tex = FALSE, digits = 8, 
    ...) 
{
    if (!is.fpt.density(obj)) 
        stop(paste(sQuote("obj"), "is not of class", shQuote("fpt.density")))
    if (report.sfptl) 
        report(attr(obj, "summary.fptl"), tex, digits)
    args <- as.list(attr(obj, "Call"))[-1]
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
        n <- args$n
    else n <- 250
    if (is.element("p", names(args))) 
        p <- args$p
    else p <- 0.2
    if (is.element("tol", names(args))) 
        tol <- args$tol
    else tol <- 1e-05
    argsFPTL <- as.list(attr(attr(obj, "summary.fptl"), "FPTLCall"))[-1]
    m1 <- nrow(attr(obj, "summary.fptl"))
    m2 <- nrow(attr(obj, "Steps"))
    logic <- c(from.t0 & (argsFPTL$t0 < attr(obj, "summary.fptl")[1, 
        1]), to.T & (attr(obj, "summary.fptl")[m1, 5] < argsFPTL$T))
    if (tex) {
        dollar <- "$"
        noindent <- "\\noindent "
        vskip <- "\\vskip 10pt "
        ldots <- "\\ldots"
        labelStep <- c("h_{", "}")
        if (from.t0) 
            t0.label <- paste(" starting from $t_{0} = ", argsFPTL$t0, 
                "$ ", sep = "")
        else t0.label <- " "
    }
    else {
        dollar <- ""
        noindent <- ""
        vskip <- ""
        ldots <- "..."
        labelStep <- c("h[", "]")
        if (from.t0) 
            t0.label <- paste("starting from t0 = ", argsFPTL$t0, 
                " ", sep = "")
        else t0.label <- ""
    }
    if (to.T) 
        T.label <- paste(rep("and ", from.t0), "until ", dollar, 
            "T = ", argsFPTL$T, dollar, " ", sep = "")
    else T.label <- ""
    h <- format((attr(obj, "summary.fptl")[, 3] - attr(obj, 
        "summary.fptl")[, 2])/n, digits=digits)
    intervals <- paste("[", format(attr(obj, "summary.fptl")[, 
        2], digits=digits), ", ", format(attr(obj, "summary.fptl")[, 
        5], digits=digits), "]", sep = "")
    if (m1 > 1) {
        cat("\n", vskip, noindent, "From the information provided by the FPTL function, we must use the following integration steps in order to approximate", 
            sep = "")
        cat("\nthe first-passage-time density:\n")
        if (tex) {
            cat("$$\\setlength{\\arraycolsep}{5pt} \n\\begin{array}{rcl}\n")
            cat(paste("h_{", 1:m1, "} = ", h, " & \\mbox{in subinterval} & ", 
                intervals, " \\\\", sep = ""), sep = "\n")
            cat("\\end{array}\n$$")
        }
        else cat(paste("\th[", format(1:m1), "] = ", h, "  in subinterval  ", 
            intervals, sep = ""), sep = "\n")
        cat("\n", vskip, noindent, "With the aim of determining these integration steps, we have divided the appropriate subintervals into ", 
            n, " parts.", sep = "")
    }
    else {
        cat("\n", vskip, noindent, "From the information provided by the FPTL function, we must use the integration step ", 
            sep = "")
        cat(dollar, labelStep[1], "1", labelStep[2], " = ", h[1], 
            dollar, " in subinterval ", dollar, intervals[1], 
            dollar, sep = "")
        cat("\nin order to approximate the first-passage-time density. With the aim of determining this integration step, we have divided the")
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
    h.labels <- paste(h.labels, " = ", format(attr(obj, 
        "Steps")[, 3], digits=digits), sep = "")
    intervals <- paste("[", format(attr(obj, "Steps")[, 
        1], digits=digits), ", ", format(attr(obj, "Steps")[, 
        2], digits=digits), "]", sep = "")
    cat("\n\n", vskip, noindent, "The f.p.t. density has been approximated ", 
        t0.label, T.label, "using a ", switch(1 + variableStep, 
            "fixed", "variable"), " integration step", sep = "")
    if (skip) 
        cat(" and avoiding the application \nof the numerical algorithm in those subintervals in which this is possible.")
    else cat(".")
    if (variableStep) {
        if (m2 > 1) {
            cat("\n\n", vskip, noindent, "For this specific application of the algorithm, we consider the following subintervals and integration steps:\n", 
                sep = "")
            if (tex) {
                cat("$$\\setlength{\\arraycolsep}{5pt} \n\\begin{array}{rcl}\n")
                cat(paste(h.labels, " & \\mbox{in subinterval} & ", 
                  intervals, " \\\\", sep = ""), sep = "\n")
                cat("\\end{array}\n$$")
            }
            else cat(paste("\t", h.labels, "  in subinterval  ", 
                intervals, sep = ""), sep = "\n")
            cat("\n", vskip, noindent, "The endlimits of the subintervals have been readjusted according to the integration steps.", 
                sep = "")
        }
        else {
            cat("\n\n", vskip, noindent, "For this specific application of the algorithm, we consider the integration step ", 
                dollar, h.labels, dollar, " in subinterval ", 
                dollar, intervals, dollar, ".", sep = "")
            cat("\nThe endlimits of the subinterval have been readjusted according to the integration step.", 
                sep = "")
        }
        if ((m1 > 1) | (any(logic))) {
            if ((m1 + sum(logic)) > 2) 
                pl <- "s"
            else pl <- ""
            cat("\n\n", vskip, noindent, "In order to determine the integration step", 
                pl, " ", sep = "")
            if (logic[1]) {
                cat(dollar, labelStep[1], "0", labelStep[2], 
                  dollar, sep = "")
                if (m1 > 1) {
                  if (logic[2]) 
                    cat(", ")
                  else cat(" and ")
                }
            }
            if (m1 > 1) {
                i <- switch(as.character(m1 - 1), `1` = "1", 
                  "i")
                index <- switch(as.character(m1 - 1), `1` = "", 
                  `2` = paste(", ", dollar, "i=1,2", dollar, 
                    ",", sep = ""), paste(", ", dollar, "i=1,", 
                    ldots, ",", m1 - 1, dollar, ",", sep = ""))
                cat(dollar, labelStep[1], i, labelStep[2], "^*", 
                  dollar, index, sep = "")
            }
            if (logic[2]) {
                if ((m1 > 1) | logic[1]) 
                  cat(" and ")
                cat(dollar, labelStep[1], m1 + 1, labelStep[2], 
                  dollar, sep = "")
            }
            cat(" we have divided the corresponding subinterval", 
                pl, " into ", trunc(n * p + 0.5), " parts.", 
                sep = "")
        }
    }
    else {
        cat("\n\n", vskip, noindent, "For this specific application of the algorithm, we consider the fixed integration step ", 
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
    x <- matrix(, nrow = length(attr(obj, "cumIntegral")), ncol = 7 - 
        tex)	
    x[, 3 - tex] <- attr(obj, "Steps")[1:length(attr(obj, 
        "cumIntegral")), 3]
    x[, 4 - tex] <- attr(obj, "cumIntegral")
    x[, 5 - tex] <- (attr(obj, "Steps")[1:length(attr(obj, "cumIntegral")), 
        2] - attr(obj, "Steps")[1:length(attr(obj, "cumIntegral")), 
        1])/attr(obj, "Steps")[1:length(attr(obj, "cumIntegral")), 
        3]
    x[, 6:7 - tex] <- attr(obj, "CPUTime")
    x[attr(obj, "skips"), 5 - tex] <- 0
    it <- sum(x[, 5 - tex])
    
    if (tex) {
	  x <- format(as.data.frame(x), digits=digits)
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
        x[, 1:2] <- attr(obj, "Steps")[1:length(attr(obj, 
            "cumIntegral")), 1:2]
	  x <- format(as.data.frame(x), digits=digits)
        x[attr(obj, "skips"), 3 - tex] <- ""	
        if (length(attr(obj, "cumIntegral")) > 1) 
            row.names(x) <- paste("Subinterval", 1:length(attr(obj, 
                "cumIntegral")))
        else row.names(x) <- "Subinterval"
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
            dollar, cumIntegral, dollar, ",", sep = "")
    }
    else {
        tstop <- paste(dollar, "t = ", format(obj$x[length(obj$x)], 
            digits=digits), dollar, sep = "")
        if ((cumIntegral > (1 - tol)) & (length(attr(obj, "cumIntegral")) < 
            m2)) 
            cat("\nThe algorithm was stopped at ", tstop, ", since the value of the cumulative integral of the approximation is ", 
                dollar, cumIntegral, switch(1 + tex, ">=", " \\geq "), 
                1 - tol, dollar, ",", sep = "")
        else cat("\nThe algorithm was stopped at ", tstop, ", and the value of the cumulative integral of the approximation is ", 
            dollar, cumIntegral, dollar, ",", sep = "")
    }
    cat("\nthe total number of iterations is ", dollar, it, dollar, 
        ", and the user time employed was ", dollar, format(sum(attr(obj, 
            "CPUTime")[, 1]), 2), dollar, " (in seconds).", sep = "")
    cat("\n\n")
}

