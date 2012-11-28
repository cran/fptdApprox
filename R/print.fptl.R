print.fptl <-
function (x, ...) 
{
    if (!is.fptl(x)) 
        stop(paste(sQuote("x"), "is not of class", shQuote("fptl")))    
    cat("\nAn object of class", shQuote("fptl"), "containing")
    cat("\n   $ x: a sequence of ", length(x$x), " values from ", 
        format(x$x[1], ...), " to ", format(x$x[length(x$x)], ...), sep = "")
    cat("\n   $ y: the values of the FPTL function on x")
    cat("\n\nCall:\n")
    print(attr(x, "Call"))
    cat("\n")
}
