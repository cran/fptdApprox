is.fptl <-
function (obj) 
{
    if (inherits(obj, "fptl") & is.list(obj) & (length(obj) == 
        2)) 
        if (all(unlist(lapply(obj, is.numeric))) & identical(names(obj), 
            c("x", "y")) & identical(names(attributes(obj)), 
            c("names", "dp", "Call", "class"))) 
            if ((length(obj$x) == length(obj$y)) & is.diffproc(attr(obj, 
                "dp")) & is.call(attr(obj, "Call"))) {
                args <- as.list(attr(obj, "Call"))[-1]
                if (all(is.element(c("dp", "t0", "T", "x0", "S"), 
                  names(args)))) {
                  if (is.element("n", names(args))) 
                    n <- args$n
                  else n <- 10000
                  return((length(obj$x) == n) & all(obj$x >= 
                    args$t0) & all(obj$x <= args$T) & all(obj$y >= 
                    0) & all(obj$y <= 1))
                }
            }
    return(FALSE)
}

