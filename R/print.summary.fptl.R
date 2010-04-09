print.summary.fptl <-
function (x, ...) 
{
    if (!is.summary.fptl(x)) 
        stop(paste(sQuote("x"), "is not of class", shQuote("summary.fptl")))
    cat("\nCall:\n")
    print(attr(x, "Call"))
    cat("\nFPTLCall:\n")
    print(attr(x, "FPTLCall"))
    cat("\n\nSlope to consider that a growing function is constant                 ", 
        attr(x, "zeroSlope"))
    cat("\nNumber of time instants from which the FPTL function starts growing   ", 
        nrow(x))
    cat("\nParameter that controls where the FPTL function begins to increase    ")
    cat("\nsignificantly                                                         ", 
        attr(x, "p0.tol"))
    cat("\n\n\nInteresting time instants:\n\n")
    y <- x[1:nrow(x), ]
    if (!is.matrix(y)) 
        y <- matrix(y, nrow = 1)
    labels <- c("t[i]", "t[i]*", "tmax[i]^-", "tmax[i]", "tmax[i]^+")
    dimnames(y) <- list(format(paste("I[", 1:nrow(y), "]", sep = "")), 
        labels)
    print(y, ...)
    cat("\n")
}

