print.summary.fptl <-
function (x, ...) 
{
    if (!is.summary.fptl(x)) 
        stop(paste(sQuote("x"), "is not of class", shQuote("summary.fptl")))

    cat("\nCall:\n")
    print(attr(x, "Call"))
    args <- as.list(attr(x, "Call"))[-1]
    if (is.element("zeroSlope", names(args))) zeroSlope <- eval(args$zeroSlope) else zeroSlope <- 0.01
    if (is.element("p0.tol", names(args))) p0.tol <- eval(args$p0.tol) else p0.tol <- 8
    if (is.element("k", names(args))) k <- eval(args$k) else k <- 3
			
    cat("\nFPTLCall:\n")
    print(attr(x, "FPTLCall"))
    cat("\n\nSlope to consider that a growing function is constant                 ", zeroSlope)
    cat("\nNumber of time instants from which the FPTL function starts growing   ", nrow(x))
    cat("\nParameter that controls where the FPTL function begins to increase    ")
    cat("\nsignificantly                                                         ", p0.tol)
    cat("\nParameter that controls where the FPTL function decrease slowly       ", k)
    cat("\n\n\nInteresting time instants:\n\n")
    y <- x[1:nrow(x), ]
    if (!is.matrix(y)) 
        y <- matrix(y, nrow = 1)
    labels <- c("t[i]", "t[i]*", "tmax[i]^-", "tmax[i]", "tmax[i]^+")
    dimnames(y) <- list(format(paste("Subinterval", 1:nrow(y), " ")), 
        labels)
    print(y, ...)
    cat("\n")
}
