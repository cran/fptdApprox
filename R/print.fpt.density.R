print.fpt.density <-
function (x, ...) 
{
    if (!is.fpt.density(x)) 
        stop(paste(sQuote("x"), "is not of class", shQuote("fpt.density")))

    y <- data.frame(matrix(, nrow = length(attr(x, "cumIntegral")), ncol = 7))
    row.names(y) <- paste("Subinterval", 1:length(attr(x,"cumIntegral")), " ")	        
    names(y) <- c("Lower end", "Upper end", "Integration step", 
            "Cumulative integral", "Iterations", "User time", 
            "System time")
    y[, 1:2] <- format(attr(x, "Steps")[1:length(attr(x, "cumIntegral")), 
        1:2], ...)

    y[, 3] <- attr(x, "Steps")[1:length(attr(x, "cumIntegral")), 
        3]    
    y[, 4] <- attr(x, "cumIntegral")

    y[, 5] <- (attr(x, "Steps")[1:length(attr(x, "cumIntegral")), 
        2] - attr(x, "Steps")[1:length(attr(x, "cumIntegral")), 
        1])/attr(x, "Steps")[1:length(attr(x, "cumIntegral")), 
        3]
    it <- format(sum(y[!attr(x, "skips"), 5]), ...)

    y[, 6:7] <- attr(x, "CPUTime")
    ut <- format(sum(y[, 6]), ...)

    y <- format(y, ...)
    y[attr(x, "skips"), 3] <- ""
    y[attr(x, "skips"), 5] <- 0
	    
    cat("\nAn object of class", shQuote("fpt.density"), "containing")
    cat(paste("\n   $ x: a sequence of", length(x$x), "time instants from", 
        y[1,1], "to", y[length(attr(x, "cumIntegral")),2], sep = " "))	
    cat("\n   $ y: the values of the approximate first-passage-time density function on sequence x")
    cat("\n\nCall:\n")
    cat(deparse(attr(x, "Call"), width.cutoff = 500), "\n")

    cat("\n\nSpecified options to apply the numerical algorithm:")
    args <- as.list(attr(x, "Call"))[-1]
    if (is.element("variableStep", names(args))) 
        variableStep <- as.logical(as.character(args$variableStep))
    else variableStep <- TRUE
    cat("\n\nVariable integration step                                                              ", 
        variableStep)
    if (is.element("from.t0", names(args))) 
        from.t0 <- as.logical(as.character(args$from.t0))
    else from.t0 <- FALSE
    cat("\nCalculate the approximation from the lower end of the interval considered              ", 
        from.t0)
    if (is.element("to.T", names(args))) 
        to.T <- as.logical(as.character(args$to.T))
    else to.T <- FALSE
    cat("\nCalculate the approximation to the upper end of the interval considered                ", 
        to.T)
    if (is.element("skip", names(args))) 
        skip <- as.logical(as.character(args$skip))
    else skip <- TRUE
    cat("\nSkip the intervals at which the FPTL function is near zero                             ", 
        skip)
    if (is.element("tol", names(args))) 
        tol <- args$tol
    else tol <- 1e-03
    # cat("\nStop the algorithm if the cumulative integral of the approximation is over    a        ",
    cat("\nStop the algorithm if the cumulative integral at appropriate time instants is over     ", 
        1 - tol)
    if (variableStep) {
        if (is.element("n", names(args))) 
            n <- args$n
        else n <- 250
        cat("\n\nNumber of points to determine optimal integration steps where the f.p.t. density")
        cat("\nfunction increases or decreases significantly according to the FPTL function           ", 
            n)        
    }
    cat("\n\n\nApplication summary of the numerical algorithm:\n\n")    	
    print(y)	
    cat("\n\nTotal number of iterations   ", it)
    cat("\nTotal user time              ", ut)
    cat("\n\n")
}
