report.summary.fptl <-
function (obj, tex = FALSE, digits = 8, ...) 
{
    if (!is.summary.fptl(obj)) 
        stop(paste(sQuote("obj"), "is not of class", shQuote("summary.fptl")))
    args.FPTL <- as.list(attr(obj, "FPTLCall"))[-1]
    args.summary <- as.list(attr(obj, "Call"))[-1]
    if (is.element("n", names(args.FPTL))) 
        n <- args.FPTL$n
    else n <- 10000
    if (is.element("zeroSlope", names(args.summary))) 
        zeroSlope <- args.summary$zeroSlope
    else zeroSlope <- 0.01
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
    }
    else {
        noindent <- ""
        vskip <- ""
        dollar <- ""
        arrayrowsep <- ""
        ldots <- "..."
        arraycolsep <- " = "
        t0.label <- "t0"
        labels1 <- c("t[", "t[", "tmax[", "tmax[", "tmax[")
        labels2 <- c("]", "]*    ", "]^-", "]  ", "]^+")
    }
    cat("\n", noindent, "The FPTL function was evaluated at ", 
        n, " points in ", sep = "")
    cat(dollar, "[", t0.label, ", T] = [", format(args.FPTL$t0, digits=digits), ", ", 
        format(args.FPTL$T, digits=digits), "]", dollar, ".", sep = "")
    cat(" If we consider that it is constant")
    cat("\nin those growth intervals where the approximate slope of the function is less than ", 
        dollar, zeroSlope, dollar, " degrees, it results that the", 
        sep = "")

    x <- format(obj, digits=digits)

    if (nrow(obj) > 1) {
        pl <- "s"		
        cat("\nfunction starts growing from points:")
        index <- 1:nrow(x)
        labels <- matrix(paste("\t", outer(labels1, index, paste, 
            sep = ""), paste(labels2, t(x), sep = arraycolsep), 
            arrayrowsep, sep = ""), nrow = 5)
        if (tex) {
            cat("\n$$\n\\setlength{\\arraycolsep}{2pt}\n\\begin{array}{rcl}\n")
            cat(labels[1, ], sep = "\n")
            cat("\\end{array}\n$$\n")
            labels[1, ] <- "$$\n\\setlength{\\arraycolsep}{2pt}\n\\begin{array}{rcl}"
            labels <- rbind(paste("\n", vskip, noindent, "For interval $I_{", 
                index, "} = [", x[, 1], ", ", c(x[-1, 1], format(args.FPTL$T, digits=digits)), 
                "]$", ":", sep = ""), labels, "\\end{array}\n$$")
        }
        else {
            cat("", labels[1, ], sep = "\n")
            labels[1, ] <- paste("\n", vskip, noindent, "For interval I[", 
                index, "] = [", x[, 1], ", ", c(x[-1, 1], format(args.FPTL$T, digits=digits)), 
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
        cat("\nfunction starts growing from point ", dollar, 
            labels1[1], "1", labels2[1], " = ", x[1, 1], dollar, 
            ".", sep = "")
        if (tex) {
            cat("\n\n", vskip, noindent, "For interval $I_{1} = [", 
                x[1, 1], ", ", format(args.FPTL$T, digits=digits), "]$:\n", sep = "")
            cat("$$\n\\setlength{\\arraycolsep}{2pt}\n\\begin{array}{rcl}\n")
            cat(paste("\t", labels1[-1], "1", labels2[-1], arraycolsep, 
                x[1, -1], arrayrowsep, sep = ""), sep = "\n")
            cat("\\end{array}\n$$\n")
        }
        else {
            cat("\n\nFor interval I[1] = [", x[1, 1], ", ", format(args.FPTL$T, digits=digits), 
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
        "10^{", -attr(obj, "p0.tol"), "}", dollar, " times the increment of the function between ", 
        dollar, labels1[1], i, labels2[1], dollar, " and ", dollar, 
        labels1[4], i, sub("\\s+$", "", labels2[4]), dollar, 
        ".", sep = "")
    cat("\n\n")
}

